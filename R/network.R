#' @include layer.R
NULL

#' leabra network class
#'
#' Class to simulate a biologically realistic network of units organized in
#' layers
#'
#' This class simulates a biologically realistic artifical neuronal network in
#' the lebra framework. It consists of several \code{layer} objects in the
#' variable (field) \code{lays} and some network-specific variables.
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @keywords data
#' @return Object of \code{\link{R6Class}} with methods for calculating changes
#'   of activity in a layer of neurons
#' @format \code{\link{R6Class}} object.
#'
#' @examples
#' # create a network with 3 layers
#' net <- network$new(list(c(5, 5), c(5, 10), c(5, 5))
#'
#' l$g_e_avg # private values cannot be accessed
#' # if you want to see alle variables, you need to use the function
#' l$get_layer_vars(show_dynamics = T, show_constants = T)
#' # if you want to see a summary of all units without constant values
#' l$get_unit_vars(show_dynamics = T, show_constants = F)
#'
#' # let us clamp the activation of the 25 units to some random values between
#' # 0.05 and 0.95
#' acts <- 0.05 + runif(25, 0, .9)
#' l$avg_act
#' l$clamp_cycle(acts)
#' l$avg_act
#' # what happened to the unit acts?
#' l$get_unit_acts()
#' # compare with acts
#' acts
#' # scaled acts are scaled by the average activity of the layer and should be
#' # smaller
#' l$get_scaled_acts()
#'
#' @field lays a list of layer objects
#' @field lrate learning rate
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new(dim_lays, w_init; g_i_gain)}}{Creates an object of this
#'   class with default parameters.
#'
#'   \code{dim_lays} list of number pairs for rows and columns of the layers
#'
#'   \code{cxn} matrix specifying connection strength between layers, if layer j
#'   sends projections to layer i, then cxn[i, j] = c > 0; 0 otherwise, c
#'   specifies the relative strength of that connection with respect to the
#'   other projections to layer i
#'
#'   \code{w_init} matrix of initial weight matrices (like a cell array in
#'   matlab), this is analogous to cxn, i.e. w_init[i, j] contains the initial
#'   weight matrix for the cxn from layer j to i
#'
#'   \code{g_i_gain} vector of inhibitory conductance gain values for every
#'   layer, this comes in handy to control overall level of inhibition of
#'   specific lays}
#'
#'   \item{\code{cycle(ext_inputs, clamp_inp)}}{cycle iterates one time step
#'   with network object
#'
#'     \code{ext_inputs} a list of matrices; ext_inputs[[i]] is a matrix that
#'     for layer i specifies the external input to each of its units. An empty
#'     matrix denotes no input to that layer. You can also use a vector instead
#'     of a matrix, because the matrix is vectorised anyway.
#'
#'     \code{clamp_inp} a binary flag;
#'     1: layers are clamped to their input value
#'     0: external inputs are summed to ecitatory conductance values}
#'
#'   \item{\code{chg_wt()}}{Changes the weights of the entire network with the
#'   XCAL learning equation}
#'
#'   \item{\code{reset(random = F)}}{ sets the activation of all units in all
#'   layers to 0, and sets all activation time averages to that value, used to
#'   begin trials from a random stationary point, the activity values may also
#'   be set to a random value
#'
#'     \code{random} logical variable, if TRUE set activation randomly between
#'     .05 and .95, if FALSE set activation to 0}
#'
#'   \item{\code{set_weights(w)}}{Sets new weights
#'
#'     \code{w} matrix of matrices (like a cell array in matlab) with new weight
#'     values}
#'
#'   \item{\code{get_unit_vars(show_dynamics = T, show_constants = F)}}{Returns
#'   a data frame with with the current state of all unit variables in the
#'   layer. Every row is a unit. You can choose whether you want dynamic values
#'   and / or constant values. This might be useful if you want to analyse what
#'   happens in units of a layer, which would otherwise not be possible, because
#'   most of the variables (fields) are private in unit class.}
#'
#'   \item{\code{get_layer_vars(show_dynamics = T, show_constants = F)}}{Returns
#'   a data frame with 1 row with the current state of the variables in the
#'   layer.  You can choose whether you want dynamic values and / or constant
#'   values. This might be useful if you want to analyse what happens in a
#'   layer, which would otherwise not be possible, because some of the variables
#'   (fields) are private in the layer class.}
#'
#'   \item{\code{get_weights()}}{Returns the complete weight matrix,
#'   \code{w[rcv, snd]} contains the weight matrix for the projections from
#'   layer \code{snd} to layer \code{rcv}. Note that this is a matrix of
#'   matrices (equivalent to a matlab cell array.)}
#'
#' }
network <-  R6::R6Class("network",
  public = list(
    initialize = function(dim_lays, cxn,
                          lrate = 0.1,
                          g_i_gain = rep(2, length(dim_lays)),
                          w_init_fun = function(x) 0.3 + 0.4 * runif(x),
                          w_init = NULL
    ){
      private$dim_lays <- dim_lays

      private$cxn <- cxn
      private$normalize_cxn()
      private$cxn_greater_zero <- apply(private$cxn, c(1, 2),
                                        function(x) ifelse(x > 0, 1, 0))

      private$g_i_gain <- g_i_gain
      self$lrate <- lrate

      private$n_lays <- length(dim_lays)

      private$create_lays()
      private$set_all_unit_numbers()

      ifelse(is.null(w_init),
             private$create_rnd_wt_matrix(w_init_fun),
             private$w_init <- w_init)

      private$set_all_w_vars()

      private$test_argument_dimensions()


      private$second_test_argument_dimensions()

      private$create_wt_for_lays()
      invisible(self)
      },

    cycle = function(ext_inputs, clamp_inp){
      private$ext_inputs <- lapply(ext_inputs, c)
      private$is_ext_inputs_valid()
      private$has_layer_ext_input <- !sapply(private$ext_inputs, is.null)
      lay_has_no_ext_input_or_is_unclamped <-
        !private$has_layer_ext_input | !clamp_inp

      if (clamp_inp == T) private$clamp_ext_inputs()

      intern_inputs <- private$find_intern_inputs()
      private$cycle_with_unclamped_lays(intern_inputs, ext_inputs,
                                        lay_has_no_ext_input_or_is_unclamped)
      invisible(self)
    },

    chg_wt = function(){
      private$updt_recip_avg_act_n()
      lapply(self$lays, function(x) x$updt_unit_avg_l())

      # Extracting the averages for all layers
      avgs <- lapply(self$lays,
                     function(x) x$get_unit_act_avgs(private$m_in_s,
                                                     private$avg_l_lrn_min,
                                                     private$avg_l_lrn_max))
      avg_l <- lapply(avgs, function(x) x$avg_l)
      avg_m <- lapply(avgs, function(x) x$avg_m)
      avg_s_with_m <- lapply(avgs, function(x) x$avg_s_with_m)
      avg_l_lrn <- lapply(avgs, function(x) x$avg_l_lrn)

      # For each connection matrix, calculate the intermediate vars
      is_cxn_null <- apply(private$cxn_greater_zero, c(1, 2),
                          function(x) ifelse(x == 0, return(NULL), return(x)))

      avg_s_with_m_snd <- private$get_snd_avg(avg_s_with_m, is_cxn_null)
      avg_s_with_m_rcv <- private$get_rcv_avg(avg_s_with_m, is_cxn_null)
      avg_m_snd <- private$get_snd_avg(avg_m, is_cxn_null)
      avg_m_rcv <- private$get_rcv_avg(avg_m, is_cxn_null)

      avg_l_rcv <- private$get_rcv_avg(avg_l, is_cxn_null)
      avg_l_lrn_rcv <- private$get_rcv_avg(avg_l_lrn, is_cxn_null)

      s_hebb <- private$m_mapply(function(x, y) x %o% y,
                         avg_s_with_m_rcv, avg_s_with_m_snd)
      m_hebb <- private$m_mapply(function(x, y) x %o% y, avg_m_rcv, avg_m_snd)

      # this is independent of sending; only the receiving neuron's long term
      # average is important for self-organized learning
      l_rcv <- private$m_mapply(function(x, y) matrix(rep(x, y), ncol = y),
                                avg_l_rcv, private$n_units_in_rcv_lays)

      avg_l_lrn_rcv <- private$m_mapply(function(x, y) matrix(rep(x, y),
                                                              ncol = y),
                                avg_l_lrn_rcv,
                                private$n_units_in_rcv_lays)
      # wt changes
      dwt_m <- private$m_mapply(function(x, y) private$get_dwt(x, y),
                                s_hebb,
                                m_hebb)
      dwt_l <- private$m_mapply(function(x, y) private$get_dwt(x, y),
                                s_hebb,
                                l_rcv)

      # multiply longerm average by individual unit scaling factor
      dwt_l <- private$m_mapply(function(x, y) x * y, dwt_l, avg_l_lrn_rcv)

      # multiply medium average by gain for medium (versus long)
      dwt_m <- apply(dwt_m, c(1, 2), function(x) x[[1]] * private$m_lrn)

      # combine both dwts and multiply by learning rate
      dwt <- private$m_mapply(function(x, y) x + y, dwt_m, dwt_l)
      dwt <- apply(dwt, c(1, 2), function(x) x[[1]] * self$lrate)

      # replace empty cells with null, then use cbind to generate the dwt for
      # every layer in a list
      dwt_null <- apply(dwt, c(1, 2), function(x)
        ifelse(isempty(x[[1]]), return(NULL), return(x[[1]])))
      dwt_list <- apply(dwt_null, 1, function(x) Reduce(cbind, x))

      private$set_new_exp_bounded_wts(dwt_list)
      lapply(self$lays, function(x) x$set_ce_weights(private$off, private$gain))
      invisible(self)
    },

    reset = function(random = F){
      lapply(self$lays, function(x) x$reset(random))
      invisible(self)
    },

    set_weights = function(w){
      # This function receives a cell array w, which is like the cell array
      # w_init in the constructor: w{i,j} is the weight matrix with the initial
      # weights for the cxn from layer j to layer i. The weights are set to the
      # values of w. This whole function is a slightly modified copypasta of the
      # constructor.

      private$w_init <- w

      w_empty <- matrix(sapply(w, isempty), nrow = nrow(w))

      private$is_dim_w_init_equal_to_dim_cxn(w)
      private$connected_lays_have_weight_matrix()
      private$non_connected_lays_do_not_have_weight_matrix()

      private$all_layer_dims_correspond_width_weight_matrix_dims()

      ## Now we set the weights
      # first find how many units project to the layer in all the network
      lay_inp_n <- apply(private$n_units_in_rcv_lays, 1, sum)

      # make one weight matrix for every layer by collapsing receiving layer
      # weights columnwise, only do this for layers that receive something at
      # all (isempty)
      wts <- apply(w, 1, function(x) Reduce("cbind", x))
      self$lays <- Map(function(x, y) {if (!isempty(y)) x$wt <- y; x},
                       self$lays, wts)

      # set the contrast-enhanced version of the weights
      lapply(self$lays, function(x) x$set_ce_weights(private$off, private$gain))
      invisible(self)
    },

    get_weights = function(){
      w <- matrix(Map(private$get_weight_for_one_connection, private$cxn,
                      private$w_index_low, private$w_index_up, self$lays),
                  ncol = ncol(private$w_index_low))
    },

    get_layer_and_unit_vars = function(){

    },

    create_inputs = function(which_layers, n_inputs, prop_active = 0.3){
      lapply(1:n_inputs, function(x)
        private$create_one_input(private$dim_lays, which_layers, prop_active))
    },

    learn_error_driven = function(inputs_minus, inputs_plus, lrate = 0.1,
                                  n_cycles_minus = 50, n_cycles_plus = 25){

      outs <- mapply(private$learn_one_pattern_error_driven,
                     inputs_minus,
                     inputs_plus,
                     pattern_number = seq(length(inputs_minus)),
                     MoreArgs = list(n_cycles_minus = n_cycles_minus,
                                     n_cycles_plus = n_cycles_plus,
                                     lrate = lrate,
                                     number_of_patterns = length(inputs_minus)))
      return(outs)
    },

    # fields -------------------------------------------------------------------
    lrate = 0.1,  # learning rate for XCAL
    lays = list()
  ),
  # private --------------------------------------------------------------------
  private = list(
    # trial
    #
    # runs several cycles with one input (for all layers)
    #
    # input is a list with the length being the number of layers
    #
    # lrate is the learning rate
    #
    # clamp_inp whether input should be clamped
    #
    # reset whether the network should be reset to a stationary point before
    # cycling (this can be extended by using the the random parameter in reset
    # if one wants random activations)
    #
    # returns invisible self
    run_trial = function(input, n_cycles, lrate = 0.1, reset = F){
      self$lrate <- lrate
      if (reset == T) self$reset()
      for(i in seq(n_cycles)) self$cycle(input, clamp_inp = T)
      invisible(self)
    },

    # learn_one_pattern_error_driven
    #
    # it does what it says, learning a pattern in error driven style
    #
    # input_minus one input for all layers, without the correct output
    #
    # input_plus one input for all layers, with correct output
    #
    # n_cycles_minus number of cycles for minus phase
    #
    # n_cycles_plus number of cycles for plus phase
    #
    # lrate learning rate
    #
    # clamp_inp whether input should be clamped
    #
    # reset whether the network should be reset to a stationary point before
    # cycling (this can be extended by using the the random parameter in reset
    # if one wants random activations)
    #
    # returns activties of all units after minus phase (before learning)
    learn_one_pattern_error_driven = function(input_minus, input_plus,
                                              pattern_number,
                                              number_of_patterns,
                                              n_cycles_minus = 50,
                                              n_cycles_plus = 25,
                                              lrate = 0.1,
                                              reset = T){
      # minus phase
      private$run_trial(input_minus, n_cycles_minus, lrate = lrate,
                        reset = reset)

      output <- lapply(self$lays, function(x) x$get_unit_acts())

      # plus phase, do not reset the network!
      private$run_trial(input_plus, n_cycles_plus, lrate = lrate, reset = F)

      # change weights
      self$chg_wt()

      # show progress
      if (pattern_number == number_of_patterns) {
        cat(".\n")
      } else {
          cat(".")
      }
      return(output)
    },

    # functions to create random inputs
    create_one_input = function(dim_lays, which_layers, prop_active = 0.3){
      prob <- c(prop_active, 1 - prop_active)
      #n_units <- lapply(dim_lays, function(x) x[1] * x[2])
      inputs <- mapply(private$create_one_input_for_one_layer,
                       private$n_units_in_lays,
                       private$inv_which(which_layers, private$n_lays),
                       MoreArgs = list(prob = prob))
      inputs
    },

    inv_which = function(indices, totlength){
      is.element(seq_len(totlength), indices)
    },

    create_one_input_for_one_layer = function(n, has_layer_input, prob){
      result <- ifelse(has_layer_input,
                       return(sample(c(.96, 0.01), n, replace = T,
                                     prob = prob)),
                       return(NULL))
      return(result)
    },

    # functions to create random weights
    create_rnd_wt_matrix = function(fun = runif){
      private$w_init <- matrix(mapply(private$create_wt_matrix_for_one_cxn,
                                      private$cxn,
                                      private$n_units_in_rcv_lays,
                                      private$n_units_in_send_lays,
                                      MoreArgs = list(fun = fun)),
                               nrow = nrow(private$cxn))
    },

    create_wt_matrix_for_one_cxn = function(cxn, n_rcv_units, n_snd_units, fun){
      ifelse(cxn > 0,
             return(matrix(fun(n_rcv_units * n_snd_units), nrow = n_snd_units)),
             return(NULL))
    },

    # first argument test in constructor ---------------------------------------
    test_argument_dimensions = function(){
      private$has_every_layer_g_i_gain()
      private$is_cxn_quadratic()
      private$is_nrow_cxn_equal_to_n_lays()
      private$is_dim_w_init_equal_to_dim_cxn()
      private$are_there_negative_cxn()
    },

    has_every_layer_g_i_gain = function(){
      if (length(private$g_i_gain) != length(private$dim_lays)){
        error <- paste("You have to specify g_i_gain for every layer, your
                       g_i_gain vector has length", length(private$g_i_gain),
                       "but you have ", private$n_lays, " layers.")
        stop(error)
      }

    },

    is_cxn_quadratic = function(){
      if (nrow(private$cxn) != ncol(private$cxn)){
        error <- (paste("Cannot create network. Connection matrix must have the
                        same number of rows and columns. It has ",
                        nrow(private$cxn), " rows and ", ncol(private$cxn), "
                        columns.", sep = ""))
        stop(error)
      }
    },

    is_nrow_cxn_equal_to_n_lays = function(){
      if (nrow(private$cxn) != length(private$dim_lays)){
        error <- (paste("Cannot create network. You have ",
                        length(private$dim_lays), " layer(s) and ",
                        nrow(private$cxn), " rows in the ", "connection matrix.
                        You need to specify cxn for every layer, so the number
                        of rows and columns in the connection matrix must equal
                        the number of layers.", sep = ""))
        stop(error)
      }
    },

    is_dim_w_init_equal_to_dim_cxn = function(){
      if (sum(dim(private$w_init) == dim(private$cxn)) < 2){
        error <- paste("Cannot create network. You have an initial matrix of
                       weight matrices with dimensions of ",
                       paste(dim(private$w_init), collapse = ", "), ", and a
                       connection matrix with dimensions of ",
                       paste(dim(private$cxn), collapse = ", "), ". They have to
                       be identical.", sep = "")
        stop(error)
      }
    },

    are_there_negative_cxn = function(){
      if (min(private$cxn) < 0){
        error <- paste("Cannot create network. Negative projection strengths
                       between layers are not allowed in cxn matrix.")
        stop(error)
      }
    },

    # other constructor functions ----------------------------------------------
    # the receiving layer's strengths are normalized, such that every receiving
    # layer has a sum of strengths from all other layers of 1, if there are any
    # connections at all.
    normalize_cxn = function(){
      private$cxn <- plyr::aaply(private$cxn,1,
                                 function(x) if(sum(x) > 0) x / sum(x) else x)
    },

    create_lays = function(){
      self$lays <- mapply(function(x, y) layer$new(x, y), private$dim_lays,
                          private$g_i_gain)
      Map(function(x, y) x$layer_number <- y, self$lays, 1:length(self$lays))
    },

    # there are a couple of variables that count how many units are in the
    # network, in the layers and for every cxn
    set_all_unit_numbers = function(){
      private$n_units_in_lays <- sapply(self$lays, function(x) x$n)
      private$n_units_in_net <- Reduce("+", private$n_units_in_lays)
      private$set_n_units_in_rcv_lays()
      private$set_n_units_in_send_lays()
    },

    set_n_units_in_rcv_lays = function(){
      result <- matrix(rep(private$n_units_in_lays,
                 length(private$n_units_in_lays)
        ), ncol = length(private$n_units_in_lays), byrow = T
      )
      result <- result * private$cxn_greater_zero
      private$n_units_in_rcv_lays <- result
    },

    set_n_units_in_send_lays = function(){
      result <- matrix(rep(private$n_units_in_lays,
                           length(private$n_units_in_lays)),
           ncol = length(private$n_units_in_lays))
      result <- result * private$cxn_greater_zero
      private$n_units_in_send_lays <- result
    },

    set_all_w_vars = function(){
      private$w_init_empty = matrix(sapply(private$w_init, isempty),
                                    nrow = nrow(private$w_init))
      private$w_index_low <- plyr::aaply(private$n_units_in_rcv_lays, 1,
                                         function(x) head(c(1, cumsum(x) + 1),
                                                          -1))
      private$w_index_up <- plyr::aaply(private$n_units_in_rcv_lays, 1,
                                        function(x) cumsum(x))
    },

    create_wt_for_lays = function(){
      wts <- apply(private$w_init, 1, function(x) Reduce("cbind", x))
      Map(function(x, y) if(!isempty(y)) x$wt <- y, self$lays, wts)
      lapply(self$lays, function(x) x$set_ce_weights(private$off, private$gain))
    },

    # second argument test in constructor---------------------------------------
    # translates matrix into character version for error output
    matrix_to_character = function(x){
      apply(which(x == T, arr.ind = T), 1,
            function(x) paste("[", paste(x, collapse = ", "), "] ", sep = ""))
    },

    second_test_argument_dimensions = function(){
      private$connected_lays_have_weigth_matrix()
      private$non_connected_lays_do_not_have_weight_matrix()
      private$all_layer_dims_correspond_with_weight_matrix_dims()
    },

    connected_lays_have_weigth_matrix = function(){
      cxn_no_w <- private$cxn > 0 & private$w_init_empty
      if (sum(cxn_no_w) > 0) {
        stop(paste("Connected layers have no weight matrix. Check thefollowing
                   row(s) and column(s) ([row, column]) in initial weight matrix
                   and connection matrix: \n"),
                   private$matrix_to_character(cxn_no_w))
      }
    },

    non_connected_lays_do_not_have_weight_matrix = function(){
      w_no_cxn <- private$cxn == 0 & !private$w_init_empty
      if (sum(w_no_cxn) > 0) {
        stop(paste("Non-Connected layers have weight matrix. Check the following
                   row(s) and column(s) ([row, column]) in initial connection
                   and weight matrix: \n"),
                   private$matrix_to_character(w_no_cxn))
      }
    },

    all_layer_dims_correspond_with_weight_matrix_dims = function(){
      w_dim_lay_dim <- mapply(
        private$one_layer_dim_corresponds_with_weight_matrix_dim,
        private$w_init,
        private$n_units_in_send_lays,
        private$n_units_in_rcv_lays)
      w_dim_lay_dim <- matrix(w_dim_lay_dim, nrow = nrow(private$w_init))

      if (sum(w_dim_lay_dim) < length(w_dim_lay_dim))
        stop(paste("Dimensions of weights and layers are inconsistent. Check the
                   following row(s) and column(s) ([row, column]) in inital
                   weight matrix and number of units in the corresponding
                   layers.", "\n"), private$matrix_to_character(w_dim_lay_dim))
    },

    one_layer_dim_corresponds_with_weight_matrix_dim =
      function(w_init,lay_n_rcv, lay_n_send){
        if (isempty(w_init) & lay_n_rcv == 0 & lay_n_send == 0){
          return(T)
        }
        sum(dim(w_init) == c(lay_n_rcv, lay_n_send)) == 2
      },

    # cycle subfunctions -------------------------------------------------------
    is_ext_inputs_valid = function(){
      private$is_ext_inputs_a_list()
      private$is_length_ext_inputs_equal_to_n_lays()
    },

    is_ext_inputs_a_list = function(){
      if (!is.list(private$ext_inputs)){
        stop(paste("First argument to cycle function should be a list of
                   matrices or vectors."))
      }
    },

    is_length_ext_inputs_equal_to_n_lays = function(){
      if (length(private$ext_inputs) != private$n_lays){
        stop(paste("Number of layers inconsistent with number of inputs in
                   network cycle. If a layer should not have an input, just use
                   NULL for that layer."))
      }
    },

    clamp_ext_inputs = function(){
      Map(function(x, y) if (!is.null(y)) x$clamp_cycle(y), self$lays,
          private$ext_inputs)
    },

    find_intern_inputs = function(){
      activeness_scaled_acts <- lapply(self$lays,
                                       function(x) x$get_unit_scaled_acts())
      rcv_lays_cxn <- unlist(apply(private$cxn, 1, list), recursive = F)
      lapply(rcv_lays_cxn, private$cxn_scale_acts, activeness_scaled_acts)
    },

    cxn_scale_acts = function(rcv_lay_cxn, activeness_scaled_acts){
      unlist(c(private$mult_list_vec(activeness_scaled_acts[rcv_lay_cxn > 0],
                                     rcv_lay_cxn[rcv_lay_cxn > 0])))
    },

    cycle_with_one_unclamped_lay = function(lay, intern_inputs, ext_inputs,
                                            lay_unclamped){
      if (lay_unclamped) lay$cycle(intern_inputs, ext_inputs)
    },

    cycle_with_unclamped_lays = function(intern_inputs, ext_inputs,
                                         lays_unclamped){
      Map(private$cycle_with_one_unclamped_lay, self$lays, intern_inputs,
          ext_inputs, lays_unclamped)
    },

    # chg_wt subfuctions -------------------------------------------------------
    #
    # get_dwt
    #
    # get weight cahnges with the "check mark" function in XCAL
    #
    # x vector of abscissa values, this is the hebbian short term activity of
    # receiving and sending unit
    #
    # th vector of threshold values with the same size as x, this is either the
    # hebbian medium term activity of receiving and sending unit or the long
    # term receiving unit activity
    #
    get_dwt = function(x, th){
      f <- rep(0, length(x))
      idx1 <- (x > private$d_thr) & (x < private$d_rev * th)
      idx2 <- x >= private$d_rev * th
      f[idx1] <- private$get_m1() * x[idx1]
      f[idx2] <- x[idx2] - th[idx2]
      return(f)
    },

    # get_weight_for_one_connection
    #
    # extracts weights for one connection
    #
    # cxn: cxn matrix (receiving is rows, sending is columns, value is strength)
    #
    # lower: lower index to extract column, which corresponds to the sending
    # layer's first unit
    #
    # upper: upper index to extract column, which corresponds to the sending
    # layer's last unit
    #
    # lay: the receiving layer (to actually get the wt)
    get_weight_for_one_connection = function(cxn, lower, upper, lay){
      ifelse(cxn > 0, return(lay$wt[, lower:upper]), return(NULL))
    },

    # updt_recip_avg_act_n
    #
    # Updates the avg_act_inert and recip_avg_act_n variables, these variables
    # update at the end of plus phases instead of cycle by cycle. This version
    # assumes full connectivity when updating recip_avg_act_n. If partial
    # connectivity were to be used, this should have the calculation in
    # WtScaleSpec::SLayActScale, in LeabraConSpec.cpp
    #
    # recip_avg_act_n is a scaling factor: more active layers are balanced, such
    # that every layer sends approximately the same average scaled activity to
    # other layers.
    #
    updt_recip_avg_act_n = function(){

      lapply(self$lays, function(x) x$updt_recip_avg_act_n())
      invisible(self)
    },

    # get_m1
    #
    # obtains the m1 factor: the slope of the left-hand line in the "check mark"
    # XCAL function.
    #
    get_m1 = function(){
      m <- (private$d_rev - 1) / private$d_rev
      return(m)
    },

    set_new_exp_bounded_wts = function(dwt_list){
      mapply(private$set_new_exp_bounded_wts_for_one_lay, dwt_list, self$lays)
    },

    set_new_exp_bounded_wts_for_one_lay = function(dwt, lay){
      if (!isempty(lay$wt)){
        lay$wt <- ifelse(dwt > 0,
                         lay$wt + (1 - lay$wt) * dwt,
                         lay$wt + lay$wt * dwt)
      }
    },

    get_snd_avg = function(avg, is_cxn_null){
      avg <- matrix(rep(avg, length(avg)), ncol = length(avg), byrow = T)
      private$m_mapply(function(x, y) x * y, is_cxn_null, avg)
    },

    get_rcv_avg = function(avg, is_cxn_null){
      avg <- matrix(rep(avg, length(avg)), ncol = length(avg), byrow = F)
      private$m_mapply(function(x, y) x * y, is_cxn_null, avg)
    },

    # general functions --------------------------------------------------------
    # mult_list_vec
    #
    # multiplies a list of vectors with constants from a vector with the length
    # of
    # the list
    mult_list_vec = function(x, y){
      mapply("*", x, y)
    },

    # m_mapply
    #
    # mapply version for matrix of matrices (matlab cell arrays), returns a
    # matrix
    m_mapply = function(fun, x, y){
      matrix(mapply(fun, x, y), ncol = ncol(x))
    },

    # fields -------------------------------------------------------------------
    dim_lays = NULL,
    cxn = NULL,
    w_init = NULL,
    n_lays = NULL,  # number of lays (number of objs in "lays")
    n_units_in_net = NULL,
    n_units_in_lays = NULL,

    # constants
    avg_l_lrn_max = 0.01, # max amount of "BCM" learning in XCAL
    avg_l_lrn_min = 0.0, # min amount of "BCM" learning in XCAL
    m_in_s = 0.1, # proportion of medium to short term avgs. in XCAL
    m_lrn = 1, # proportion of error-driven learning in XCAL
    d_thr = 0.0001, # threshold for XCAL "check mark" function
    d_rev = 0.1, # reversal value for XCAL "check mark" function
    off = 1, # "offset" in the SIG function for contrast enhancement
    gain = 6, # gain in the SIG function for contrast enhancement

    g_i_gain = 2, # g_i_gain for layers, to control overall inhibition in a specific layer
    avg_l_lrn = list(),
    #dependent
    m1 = NULL, # the slope in the left part of XCAL's "check mark"
    cxn_greater_zero = matrix(), # binary version of cxn
    n_units_in_rcv_lays = matrix(), # number of units of receiving layer in cxn matrix
    n_units_in_send_lays = matrix(), # number of units of sending layer in cxn matrix
    w_index_low = matrix(), # lower index to extract weights from layer matrix
    # in cxn format
    w_index_up = matrix(), # upper index to extract weights from layer matrix in
    # cxn format
    w_init_empty = NULL,
    ext_inputs = NULL,
    has_layer_ext_input = NULL
  )
)
