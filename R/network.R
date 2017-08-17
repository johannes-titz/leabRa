#' @include layer.R


# leabra network class----------------------------------------------------------
#' leabra network class
#'

#' Class to simulate a biologically realistic layer of neurons
#'
#' This class simulates a biologically realistic layer of neurons in the lebra
#' framework. It consists of several \code{unit} objects in the variable (field)
#' \code{units} and some layer-specific variables.
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @keywords data
#' @return Object of \code{\link{R6Class}} with methods for calculating changes of activity in a layer of neurons
#' @format \code{\link{R6Class}} object.
#'
#' @examples
#' l <- layer$new(c(5, 5)) # create a 5 x 5 layer with default leabra values
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
# TODO-----------------------------------------
#' # let us run 10 cycles with unclamped activation and output the activation
#' # produced because of changes in conductance
#' l <- unit$new()
#' cycle_number <- 1:10
#' result <- lapply(cycle_number, function(x)
#'                  l$cycle(intern_input = )$get_vars())
#' # make a data frame out of the list
#' result <- plyr::ldply(result)
#' # plot act
#' plot(result$act, type = "b", xlab = "cycle", ylab = "act")
#' # add conductance g_e to plot, should approach g_e_raw
#' lines(result$g_e, type = "b", col = "blue")
#'
#' @field lays a list of layer objects
#' @field lrate learning rate
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new()}}{Creates an object of this class with default
#'   parameters.
#'
#'   \code{dim_lays} list of number pairs for rows and columns of the layers
#'
#'   \code{cxn} matrix specifying connection strength between layers, if layer j
#' sends projections to layer i, then cxn[i, j] = c > 0; 0 otherwise, c
#' specifies the relative strength of that connection with respect to the other
#' projections to layer i
#'
#'   \code{w_init} matrix of initial weight matrices (3d array!), this is analogous to cxn, i.e. w_init[i, j] contains the initial weight matrix for the cxn from
#' layer j to i
#'
#'   \code{gi} vector of inhibitory conductance gain values for every layer, this comes in handy to control
#' overall level of inhibition of specific lays}
#'
#'   \item{\code{cycle(ext_inputs, clamp_inp)}}{cycle iterates one time step with network object
#'
#'   \code{ext_inputs} a list of matrices; ext_inputs[[i]] is a matrix that for layer i specifies the external input to each of its units. An empty matrix denotes no input to that layer. You can also use a vector instead of a matrix, because the matrix is vectorised anyway.
#'   \code{clamp_inp} a binary flag; 1: layers are clamped to their input value
#'  0: external inputs are summed to ecitatory conductance values.}
#'
#'   \item{\code{chg_wt()}}{Changes the weights of the entire network with the XCAL learning equation}
#'
#'   \item{\code{cycle(intern_input, ext_input)}}{Iterates one time step with
#'   layer object.
#'
#'   \code{intern_input} Vector with inputs from all other layers. Each input
#'   has already been scaled by a reciprocal function of the number of active
#'   units (recip_avg_act_n) of the sending layer and by the connection strength
#'   between the receiving and sending layer. The weight matrix is multiplied
#'   with this input vector to get the excitatory conductance for each unit in
#'   the layer.
#'
#'   \code{ext_input} Vector with inputs not coming from another layer, with
#'   length equalling the number of units in this layer. If empty, no external
#'   inputs are processed. If the external inputs are not clamped, this is
#'   actually an excitatory conductance value, which is added to the conductance
#'   produced by the internal input and weight matrix.
#'   }
#'
#'   \item{\code{clamp_cycle(acts)}}{Iterates one time step with layer
#'   object with clamped acts, meaning that acts are
#'   instantenously set without time integration)
#'
#'   \code{acts} acts you want to clamp to the units in the layer}
#'
#' \item{\code{get_unit_act_avgs()}}{Returns a list with the short, medium and
#' long term activation averages of all units in the layer as vectors. The super
#' short term average is not returned, and the long term average is not updated
#' before being returned. These averages are used by the network class to
#' calculate weight changes.}
#'
#' \item{\code{updt_unit_avg_l()}}{Updates the long-term average (avg_l) of all
#' units in the layer, usually done after a plus phase.}
#'
#' \item{\code{updt_recip_avg_act_n()}}{Updates the avg_act_inert and
#' recip_avg_act_n variables, these variables update at the end of plus phases
#' instead of cycle by cycle. This version assumes full connectivity when
#' updating recip_avg_act_n.}
#'
#' \item{\code{reset(random = F)}}{Sets the activation and activation averages
#' of all units to 0. Used to begin trials from a stationary point.
#'
#' \code{random} Logical variable, if TRUE the activation ist set randomly
#' between .05 and .95 for every unit instead of 0.}
#'
#' \item{\code{set_ce_weights()}}{Sets contrast enhanced weight values}
#'
#' \item{\code{get_unit_vars(show_dynamics = T, show_constants = F)}}{Returns a data frame with with the current state of all unit variables in the layer. Every row is a unit. You can choose whether you want dynamic values and / or constant values. This might be useful if you want to analyse what happens in units of a layer, which would otherwise not be possible, because most of the variables (fields) are private in unit class.}
#' \item{\code{get_layer_vars(show_dynamics = T, show_constants = F)}}{Returns a
#' data frame with 1 row with the current state of the variables in the layer.
#' You can choose whether you want dynamic values and / or constant values. This
#' might be useful if you want to analyse what happens in a layer, which would
#' otherwise not be possible, because some of the variables (fields) are private
#' in the layer class.}

network <-  R6::R6Class("network",
  public = list(
    initialize = function(dim_lays, cxn, w_init, gi = rep(2, length(dim_lays)),
                          lrate = 0.1){
      private$cxn <- cxn
      private$has_every_layer_gi(gi, dim_lays)
      private$test_argument_dimensions(dim_lays, w_init)
      private$normalize_rcv_cxn()

      private$n_lays <- length(dim_lays)

      private$create_layers(dim_lays, gi)

      private$is_cxn_greater_zero <- apply(private$cxn, c(1, 2),
                                   function(x) ifelse(x > 0, 1, 0))

      private$set_all_unit_numbers()
      #private$set_all_w()

      private$w_init_empty <- matrix(sapply(w_init, isempty), nrow = nrow(w_init))

      w_index_low <- t(apply(private$number_of_units_in_receiving_layers, 1,
                             function(x) head(c(1, cumsum(x) + 1), -1)))
      w_index_up <- t(apply(private$number_of_units_in_receiving_layers, 1,
                            function(x) cumsum(x)))

      #test_argument_dimensions2

      ## Second test of argument dimensions
      private$connected_layers_have_weigth_matrix()
      private$non_connected_layers_do_not_have_weight_matrix()

      # test 3
      # checks whether receiving and sending number of units in weight matrix
      # is correct for every layer
      check_w_lay_dim <- function(w_init, lay_n_recv, lay_n_send){
        if (isempty(w_init) & lay_n_recv == 0 & lay_n_send == 0){
          return(F)
        }
        sum(dim(w_init) == c(lay_n_recv, lay_n_send)) != 2
      }

      w_dim_lay_dim <- mapply(check_w_lay_dim, w_init, private$number_of_units_in_sending_layers, private$number_of_units_in_receiving_layers)
      w_dim_lay_dim <- matrix(w_dim_lay_dim, nrow = nrow(w_init))

      if (sum(w_dim_lay_dim) > 0)
        stop(paste("Dimensions of weights and layers are inconsistent. Check the",
                   "following row(s) and column(s) ([row, column]) in inital",
                   "weight matrix and number of units in the corresponding layers.",
                   "\n"),
             private$matrix_to_character(w_dim_lay_dim))

      ## Setting the inital weights for each layer
      lay_inp_n <- apply(private$cxn, 1, function(x) sum(c(0, private$number_of_units_in_layers[x > 0])))

      # cbind weights row-wise (row receives inputs from columns) to have only
      # one weight matrix for every layer, (of course only for non-empty w_init
      # elements)
      wts <- apply(w_init, 1, function(x) Reduce("cbind", x))
      self$lays <- Map(function(x, y) {if(!isempty(y)) x$wt <- y; x}, self$lays, wts)

      # set the contrast-enhanced version of the weights
      lapply(self$lays, function(x) x$set_ce_weights(private$off, private$gain))

      private$n_lays <- n_lays
      private$gi <- gi
      private$w_index_low <- w_index_low
      private$w_index_up <- w_index_up
      self$lrate <- lrate
      invisible(self)
      },
    #' cycle iterates one time step with network object
    #'
    #' @param ext_inputs a list of matrices; ext_inputs[[i]] is a matrix that for layer i
    #' specifies the external input to each of its units. An empty matrix denotes
    #' no input to that layer. You can also use a vector instead of a matrix,
    #' because the matrix is vectorised anyway.
    #' @param clamp_inp a binary flag; 1: lays are clamped to their input value
    #' 0: inputs summed to netins.
    #' @rdname network
    #'  this will be cycle with clamped input and we will write another one
    #'  that has no clamped input later
    cycle = function(ext_inputs = "list", clamp_inp){

      private$ext_inputs <- lapply(ext_inputs, c)
      private$is_ext_inputs_valid()
      private$has_layer_ext_input <- !sapply(private$ext_inputs, is.null)
      has_layer_ext_input_and_is_clamped <- private$has_layer_ext_input & clamp_inp

      if (clamp_inp == T) private$set_all_layers_to_ext_inputs()

      scaled_acts <- lapply(self$lays, function(x) x$get_unit_scaled_acts())

      ## For each unclamped layer, we put all its scaled inputs in one
      #  column vector, and call its cycle function with that vector

      # find internal input, the list recv has the receiving connections for every
      # layer; then simply multiply the scaled_acts with the recv connections and
      # unlist into a vector, now you have a vectorised internal input that is
      # already scaled with the connection strength
      #
      # here we have to remove the connections that are 0!
      # vectorise the matrices!

      recv <- unlist(apply(private$cxn, 1, list), recursive = F)

      get_intern_input <- function(recv, scaled_acts){
        unlist(c(mult_list_vec(scaled_acts[recv > 0], recv[recv > 0])))
      }

      intern_input <- lapply(recv, get_intern_input, scaled_acts)

      # if the layer is not clamped, run cycle with inputs, otherwise return layer
      cycle_non_clamp <- function(lay, intern_input, ext_input, lay_clamped){
        if (lay_clamped) lay$cycle(intern_input, ext_input)
      }
      Map(cycle_non_clamp, self$lays, intern_input, ext_inputs,
          !has_layer_ext_input_and_is_clamped)
      invisible(self)
    },

    #' chg_wt changes the weights of the entire network with the XCAL learning
    #' equation
    #'
    #' @param obj
    #' @rdname network
    chg_wt = function(){
      private$updt_recip_avg_act_n() # update averages used for net input scaling

      # updating the long-term averages
      lapply(self$lays, function(x) x$updt_unit_avg_l())

      ## Extracting the averages for all layers
      avgs <- lapply(self$lays, function(x) x$get_unit_act_avgs(private$m_in_s, private$avg_l_lrn_min, private$avg_l_lrn_max))
      avg_l <- lapply(avgs, function(x) x$avg_l)
      avg_m <- lapply(avgs, function(x) x$avg_m)
      avg_s_with_m <- lapply(avgs, function(x) x$avg_s_with_m)
      avg_l_lrn <- lapply(avgs, function(x) x$avg_l_lrn)

      ## For each connection matrix, calculate the intermediate vars.
      # make cxn-format versions of above averages
      # we need is_cxn_greater_zero with no connection specified with NULL
      cxn_b_null <- apply(private$is_cxn_greater_zero, c(1, 2),
                          function(x) ifelse(x == 0, return(NULL), return(x)))

      # just repeat the avgs row or column-wise and then set those to NULL
      # where there is no connection
      get_snd <- function(avg, cxn_b_null){
        avg <- matrix(rep(avg, length(avg)), ncol = length(avg), byrow = T)
        m_mapply(function(x, y) x * y, cxn_b_null, avg)
      }

      get_rcv <- function(avg, cxn_b_null){
        avg <- matrix(rep(avg, length(avg)), ncol = length(avg), byrow = F)
        m_mapply(function(x, y) x * y, cxn_b_null, avg)
      }

      avg_s_with_m_snd <- get_snd(avg_s_with_m, cxn_b_null)
      avg_s_with_m_rcv <- get_rcv(avg_s_with_m, cxn_b_null)
      avg_m_snd <- get_snd(avg_m, cxn_b_null)
      avg_m_rcv <- get_rcv(avg_m, cxn_b_null)

      avg_l_snd <- get_snd(avg_l, cxn_b_null)
      avg_l_rcv <- get_rcv(avg_l, cxn_b_null)
      avg_l_lrn_rcv <- get_rcv(avg_l_lrn, cxn_b_null)

      s_hebb <- m_mapply(function(x, y) x %o% y, avg_s_with_m_rcv, avg_s_with_m_snd)
      m_hebb <- m_mapply(function(x, y) x %o% y, avg_m_rcv, avg_m_snd)

      # this is independent of sending; only the receiving neuron's long term
      # average is important; so just repeat the vector for every incoming
      # neuron
      l_recv <- m_mapply(function(x, y) matrix(rep(x, y), ncol = y), avg_l_rcv,
                         private$number_of_units_in_receiving_layers)
      avg_l_lrn_rcv <- m_mapply(function(x, y) matrix(rep(x, y), ncol = y),
                                avg_l_lrn_rcv, private$number_of_units_in_receiving_layers)

      dwt_m <- m_mapply(function(x, y) self$get_dwt(x, y), s_hebb, m_hebb)
      dwt_l <- m_mapply(function(x, y) self$get_dwt(x, y), s_hebb, l_recv)

      # multiply longerm average by individual unit scaling factor
      dwt_l <- m_mapply(function(x, y) x * y, dwt_l, avg_l_lrn_rcv)
      # mulitply medium average by gain for medium (versus long)
      dwt_m <- apply(dwt_m, c(1, 2), function(x) x[[1]] * private$m_lrn)
      # combine both dwts and multiply by learning rate
      dwt <- m_mapply(function(x, y) x+y, dwt_m, dwt_l)
      dwt <- apply(dwt, c(1,2), function(x) x[[1]] * self$lrate)

      # replace empty cells with null, then use cbind to generate the dwt for
      # every layer in a list
      dwt_null <- apply(dwt, c(1, 2), function(x)
        ifelse(isempty(x[[1]]), return(NULL), return(x[[1]])))
      dwt_list <- apply(dwt_null, 1, function(x) Reduce(cbind, x))

      # todo 2017-04-26 jt: probably this is wrong in matlab code, write the
      # author, it should be dwt > 0 and not wt > 0
      # weights have to be bounded exponentially
      exp_weights <- function(dwt, lay){
        wt <- lay$wt
        if (!isempty(wt)){
          new_wt <- ifelse(dwt > 0, wt + (1 - wt) * dwt, wt + wt * dwt)
          lay$wt <- new_wt
        }
        return(lay)
      }

      self$lays <- mapply(exp_weights, dwt_list, self$lays)

      # now set the contrast-enhanced version of the weights
      lapply(self$lays, function(x) x$set_ce_weights(private$off, private$gain))
      invisible(self)
    },

    #' get_dwt get weight changes with the "check mark" function in XCAL
    #'
    #' @param x vector of abscissa values
    #' @param th vector of threshold values with the same size as x
    #' @rdname network
    get_dwt = function(x, th){
      f <- rep(0, length(x))
      idx1 <- (x > private$d_thr) & (x < private$d_rev * th)
      idx2 <- x >= private$d_rev * th
      f[idx1] <- private$get_m1() * x[idx1]
      f[idx2] <- x[idx2] - th[idx2]
      return(f)
    },

    #' reset sets the activation of all units in all layers to a random value, and
    #' sets all activation time averages to that value, used to begin trials from a
    #' random stationary point, the activity values may also be set to zero
    #'
    #' @param random logical variable, if TRUE set activation randomly between .05
    #' and .95, if FALSE set activation to 0
    #' @rdname network
    reset = function(random = F){
      lapply(self$lays, function(x) x$reset(random))
      invisible(self)
    },

    #' set_weights sets new weights
    #'
    #' @param w matrix of matrices with new weight values
    #' @rdname network
    set_weights = function(w){
      # This function receives a cell array w, which is like the cell
      # array w_init in the constructor: w{i,j} is the weight matrix with
      # the initial weights for the cxn from layer j to layer
      # i. The weights are set to the values of w.
      # This whole function is a slightly modified copypasta of the
      # constructor.

      w_empty <- matrix(sapply(w, isempty), nrow = nrow(w))

      # translates matrix into character version for error output
      matrix_to_character <- function(x){
        apply(which(x == T, arr.ind = T), 1,
              function(x) paste("[", paste(x, collapse = ", "), "] ", sep = ""))
      }

      ## First we test the dimensions of w
      if (sum(dim(w) == dim(private$cxn)) < 2){
        stop(paste("Cannot set new weights. You have a new matrix of weight ",
                   "matrices with dimensions of ", paste(dim(w), collapse = ", "),
                   " and a connection matrix with dimensions of ",
                   paste(dim(objcxn), collapse = ", "), ". They have ",
                   "to be identical.", sep = ""))
      }

      # test 1
      cxn_no_w <- private$cxn > 0 & w_empty
      if (sum(cxn_no_w) > 0) {
        stop(paste("Connected layers have no weight matrix. Check the",
                   "following row(s) and column(s) ([row, column]) in new ",
                   "weight matrix and connection matrix: \n"),
             private$matrix_to_character(cxn_no_w))
      }
      # test 2
      w_no_cxn <- private$cxn == 0 & !w_empty
      if (sum(w_no_cxn) > 0)
        stop(paste("Non-Connected layers have weight matrix. Check the",
                   "following row(s) and column(s) ([row, column]) in initial ",
                   "connection and new weight matrix: \n"),
             private$matrix_to_character(w_no_cxn))

      # test 3
      # checks whether receiving and sending number of units in weight matrix
      # is correct for every layer
      check_w_lay_dim <- function(w_init, lay_n_recv, lay_n_send){
        if (isempty(w_init) & lay_n_recv == 0 & lay_n_send == 0){
          return(F)
        }
        sum(dim(w_init) == c(lay_n_recv, lay_n_send)) != 2
      }

      w_dim_lay_dim <- mapply(check_w_lay_dim, w_init, number_of_units_in_sending_layers, number_of_units_in_receiving_layers)
      w_dim_lay_dim <- matrix(w_dim_lay_dim, nrow = nrow(w_init))

      if (sum(w_dim_lay_dim) > 0)
        stop(paste("Dimensions of weights and layers are inconsistent. Check the",
                   "following row(s) and column(s) ([row, column]) in new",
                   "weight matrix and number of units in the corresponding layers.",
                   "\n"),
             private$matrix_to_character(w_dim_lay_dim))

      ## Now we set the weights
      # first find how many units project to the layer in all the network
      lay_inp_n <- apply(private$number_of_units_in_receiving_layers, 1, sum)

      # make one weight matrix for every layer by collapsing receiving layer
      # weights columnwise, only do this for layers that receive something at all
      # (isempty)
      wts <- apply(w, 1, function(x) Reduce("cbind", x))
      self$lays <- Map(function(x, y) {if (!isempty(y)) x$wt <- y; x}, self$lays, wts)

      # set the contrast-enhanced version of the weights
      lapply(self$lays, function(x) x$set_ce_weights(private$off, private$gain))
      invisible(self)
    },

    #' get_weights returns a matrix of weight matrices, w[rcv, snd] contains the
    #' weight matrix for the projections from layer snd to layer rcv
    #'
    #' @rdname network
    get_weights = function(){
      # extract for connected layers the weights with the index
      helper <- function(cxn, lower, upper, lay){
        ifelse(cxn > 0, return(lay$wt[, lower:upper]), return(NULL))
      }
      w <- matrix(Map(helper, private$cxn, private$w_index_low, private$w_index_up, self$lays),
                  ncol = ncol(private$w_index_low))
    },
    lrate = 0.1,  # learning rate for XCAL
    lays = list()
    ),
 # private ---------------------------------------------------------------------
  private = list(
    # updates the acts_p_avg and pct_act_scale variables for all lays.
    # These variables update at the end of plus phases instead of
    # cycle by cycle. The avg_l values are not updated here.
    # This version assumes full connectivity when updating
    # pct_act_scale. If partial connectivity were to be used, this
    # should have the calculation in WtScaleSpec::SLayActScale, in
    # LeabraConSpec.cpp
    updt_recip_avg_act_n = function(){

      lapply(self$lays, function(x) x$updt_recip_avg_act_n())
      invisible(self)
    },

    # get_m1 obtains the m1 factor: the slope of the left-hand line in the
    # "check mark" XCAL function.
    #
    get_m1 = function(){
      m <- (private$d_rev - 1) / private$d_rev
      return(m)
    },

    has_every_layer_gi = function(gi, dim_lays){
      if (length(gi) != length(dim_lays)){
        error <- paste("You have to specify gi for every layer, your gi vector
                       has length", length(gi), "but you have ", n_lays,
                       " layers.")
        stop(error)
      }

    },

    test_argument_dimensions = function(dim_lays, w_init){
        private$is_cxn_quadratic()
        private$is_nrow_cxn_equal_to_n_lays(dim_lays)
        private$is_dim_w_init_equal_to_dim_cxn(w_init)
        private$are_there_negative_cxn()
    },

    is_cxn_quadratic = function(){
      if (nrow(private$cxn) != ncol(private$cxn)){
        error <- (paste("Cannot create network. Connection matrix must have the
                        same number of rows and columns. It has ", nrow(private$cxn), "
                        rows and ", ncol(private$cxn), " columns.", sep = ""))
        stop(error)
      }
    },

    is_nrow_cxn_equal_to_n_lays = function(dim_lays){
      if (nrow(private$cxn) != length(dim_lays)){
        error <- (paste("Cannot create network. You have ", length(dim_lays), " layer(s) and
                        ", nrow(private$cxn), " rows in the ", "connection matrix. You
                        need to specify cxn for every layer, so the number of
                        rows and columns in the connection matrix must equal the
                        number of layers.", sep = ""))
        stop(error)
      }
    },

    is_dim_w_init_equal_to_dim_cxn = function(w_init){
      if (sum(dim(w_init) == dim(private$cxn)) < 2){
        error <- paste("Cannot create network. You have an initial matrix of
                        weight matrices with dimensions of ",
                        paste(dim(w_init), collapse = ", "),
                        ", and a connection matrix with dimensions of ",
                        paste(dim(private$cxn), collapse = ", "), ". They have to be
                        identical.", sep = "")
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

    normalize_rcv_cxn = function(){
      private$cxn <- t(apply(private$cxn, 1, function(x) if(sum(x) > 0) x / sum(x) else x))
    },

    create_layers = function(dim_lays, gi){
      self$lays <- mapply(function(dim_lays, gi) layer$new(dim_lays, gi), dim_lays, gi)
      Map(function(x, y) x$layer_number <- y, self$lays, 1:length(self$lays))
    },

    set_all_unit_numbers = function(){
      private$number_of_units_in_layers <- sapply(self$lays, function(x) x$n)
      private$number_of_units_in_net <- Reduce("+", private$number_of_units_in_layers)
      private$get_number_of_units_in_receiving_layers()
      private$get_number_of_units_in_sending_layers()
    },

    get_number_of_units_in_receiving_layers = function(){
      result <- matrix(rep(private$number_of_units_in_layers,
                 length(private$number_of_units_in_layers)
        ), ncol = length(private$number_of_units_in_layers), byrow = T
      )
      result <- result * private$is_cxn_greater_zero
      private$number_of_units_in_receiving_layers <- result
    },

    get_number_of_units_in_sending_layers = function(){
      result <- matrix(rep(private$number_of_units_in_layers,
                           length(private$number_of_units_in_layers)),
           ncol = length(private$number_of_units_in_layers))
      result <- result * private$is_cxn_greater_zero
      private$number_of_units_in_sending_layers <- result
    },

    # translates matrix into character version for error output
    matrix_to_character = function(x){
      apply(which(x == T, arr.ind = T), 1,
            function(x) paste("[", paste(x, collapse = ", "), "] ", sep = ""))
    },

    connected_layers_have_weigth_matrix = function(){
      cxn_no_w <- private$cxn > 0 & private$w_init_empty
      if (sum(cxn_no_w) > 0) {
        stop(paste("Connected layers have no weight matrix. Check the",
                   "following row(s) and column(s) ([row, column]) in initial ",
                   "weight matrix and connection matrix: \n"),
             private$matrix_to_character(cxn_no_w))
      }
    },

    non_connected_layers_do_not_have_weight_matrix = function(){
      w_no_cxn <- private$cxn == 0 & !private$w_init_empty
      if (sum(w_no_cxn) > 0)
        stop(paste("Non-Connected layers have weight matrix. Check the",
                   "following row(s) and column(s) ([row, column]) in initial ",
                   "connection and weight matrix: \n"),
             private$matrix_to_character(w_no_cxn))
    },

    is_ext_inputs_valid = function(){
      private$is_ext_inputs_a_list()
      private$is_length_ext_inputs_equal_to_number_of_layers()
    },

    is_ext_inputs_a_list = function(){
      if (!is.list(private$ext_inputs)){
        stop(paste("First argument to cycle function should be a list of",
                   "matrices or vectors."))
      }
    },

    is_length_ext_inputs_equal_to_number_of_layers = function(){
      if (length(private$ext_inputs) != private$n_lays){
        stop(paste("Number of layers inconsistent with number of inputs in",
                   "network cycle. If a layer should not have an input, just use NULL",
                   "for that layer."))
      }
    },

    set_all_layers_to_ext_inputs = function(){
      Map(function(x, y) if (!is.null(y)) x$clamp_cycle(y), self$lays,
          private$ext_inputs)
    },


    # fields --------------------------
    dim_lays = list(),
    cxn = matrix(),
    w_init = matrix(),
    n_lays = NULL,  # number of lays (number of objs in "lays")
    number_of_units_in_net = NULL,
    number_of_units_in_layers = NULL,

    # constants
    avg_l_lrn_max = 0.01, # max amount of "BCM" learning in XCAL
    avg_l_lrn_min = 0.0, # min amount of "BCM" learning in XCAL
    m_in_s = 0.1, # proportion of medium to short term avgs. in XCAL
    m_lrn = 1, # proportion of error-driven learning in XCAL
    d_thr = 0.0001, # threshold for XCAL "check mark" function
    d_rev = 0.1, # reversal value for XCAL "check mark" function
    off = 1, # "offset" in the SIG function for contrast enhancement
    gain = 6, # gain in the SIG function for contrast enhancement

    gi = 2, # gi for layers, to control overall inhibition in a specific layer
    avg_l_lrn = list(),
    #dependent
    m1 = NULL, # the slope in the left part of XCAL's "check mark"
    is_cxn_greater_zero = matrix(), # binary version of cxn
    number_of_units_in_receiving_layers = matrix(), # number of units of receiving layer in cxn matrix
    number_of_units_in_sending_layers = matrix(), # number of units of sending layer in cxn matrix
    w_index_low = matrix(), # lower index to extract weights from layer matrix
    # in cxn format
    w_index_up = matrix(), # upper index to extract weights from layer matrix in
    # cxn format
    w_init_empty = NULL,
    ext_inputs = NULL,
    has_layer_ext_input = NULL
  )
)
