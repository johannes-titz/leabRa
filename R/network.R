#' @include layer.R

# leabra network class----------------------------------------------------------
#' leabra network class
#'
#' @slot dim_lays list of 2d vectors that contains number of rows and columns
#' of lays
#' @slot cxn matrix specifying connection strength between lays, if
#' layer j sends projections to layer i, then cxn[i, j] = c > 0; 0
#' otherwise, c specifies the relative strength of that connection with respect
#' to the other projections to layer i
#' @slot w_init matrix of initial weight matrices, this is analogous to cxn,
#' i.e. w_init[i, j] contains the initial weight matrix for the cxn from
#' layer j to i
#' @slot gi vector of gi values for every layer, this comes in handy to control
#' overall level of inhibition of specific lays
#' @slot off "offset" in the SIG function for contrast enhancement
#' @slot gain gain in the SIG function for contrast enhancement
#' @slot lrate learning rate
network <-  R6::R6Class("network",
  public = list(
    initialize = function(dim_lays = "list", cxn = "matrix", w_init = "matrix",
                          gi = rep(2, length(dim_lays)), off = 1, gain = 6,
                          lrate = 0.1){
      # constructor
      # set standard value for gi of 2 and check whether length of gi is equal
      # to length of layers
      # if (isempty(gi)) gi <- rep(2, n_lays)
      if (length(gi) != n_lays){
        stop(paste("You have to specify gi for every layer, your gi vector has",
                   "length", length(gi), "but you have", n_lays, "layers"))
      }

      ## Initial test of argument dimensions
      if (nrow(cxn) != ncol(cxn)){
        stop(paste("Cannot create network. Connection matrix must have the same ",
                   "number of rows and columns. It has ", nrow(cxn), " rows and ",
                   ncol(cxn), " columns.", sep = ""))
      }
      if (nrow(cxn) != n_lays){
        stop(paste("Cannot create network. You have ", n_lays, " lays and ",
                   nrow(cxn), " rows in the ", "connection matrix. You need to specify ",
                   "cxn for every layer, so the number of rows in the ",
                   "connection matrix must equal the number of lays",
                   sep = ""))
      }
      if (sum(dim(w_init) == dim(cxn)) < 2){
        stop(paste("Cannot create network. You have an initial matrix of weight ",
                   "matrices with dimensions of ", paste(dim(w_init), collapse = ", "),
                   " and a connection matrix with dimensions of ",
                   paste(dim(cxn), collapse = ", "), ". They have ",
                   "to be identical.", sep = ""))
      }
      if (min(cxn) < 0){
        stop(paste("Cannot create network. Negative projection strengths between",
                   "lays are not allowed in cxn matrix."))
      }

      ## Normalizing the rows of "cxn" so they add to 1
      cxn <- t(apply(cxn, 1, function(x) if(sum(x) > 0) x / sum(x) else x))

      ##
      # list that contains layer objects
      net_lays <- mapply(function(dim_lays, gi) layer$new(dim_lays, gi), dim_lays, gi)
      # binary cxn
      cxn_b <- apply(cxn, c(1, 2), function(x) ifelse(x > 0, 1, 0))
      # som other useful stuff
      lays_n <- sapply(net_lays, function(x) x$n)
      net_n_units <- Reduce("+", lays_n)
      n_lays <- length(dim_lays)

      ## Second test of argument dimensions
      # these are matrices like w_init and cxn to ease testing
      lays_n_recv <- matrix(rep(lays_n, length(lays_n)), ncol = length(lays_n),
                            byrow = T) * cxn_b
      lays_n_send <- matrix(rep(lays_n, length(lays_n)),
                            ncol = length(lays_n)) * cxn_b
      w_init_empty <- matrix(sapply(w0, isempty), nrow = nrow(w0))

      # upper and lower index limits to extract specific weights from specific
      # layers; usually you just have one whole weight matrix for every layer,
      # which is more convenient to operate on; but if you want specific weights
      # you need to extract it from this layer weight matrix, we use shifted head
      # to accomplish this, basically you just look how many incoming
      # connections there are for every row (every receiving layer), then you
      # sum them up (cumsum) add 1; and always begin with index 1
      # for lower limit you obviously do not need the last value, but for upper
      # limit (where you do not need the leading 1)
      w_index_low <- t(apply(lays_n_recv, 1, function(x)
        head(c(1, cumsum(x) + 1), -1)))
      w_index_up <- t(apply(lays_n_recv, 1, function(x) cumsum(x)))

      # translates matrix into character version for error output
      which_matrix <- function(x){
        apply(which(x == T, arr.ind = T), 1,
              function(x) paste("[", paste(x, collapse = ", "), "] ", sep = ""))
      }

      # test 1
      cxn_no_w <- cxn > 0 & w_init_empty
      if (sum(cxn_no_w) > 0) {
        stop(paste("Connected layers have no weight matrix. Check the",
                   "following row(s) and column(s) ([row, column]) in initial ",
                   "weight matrix and connection matrix: \n"),
             which_matrix(cxn_no_w))
      }

      # test 2
      w_no_cxn <- cxn == 0 & !w_init_empty
      if (sum(w_no_cxn) > 0)
        stop(paste("Non-Connected layers have weight matrix. Check the",
                   "following row(s) and column(s) ([row, column]) in initial ",
                   "connection and weight matrix: \n"),
             which_matrix(w_no_cxn))

      # test 3
      # checks whether receiving and sending number of units in weight matrix
      # is correct for every layer
      check_w_lay_dim <- function(w_init, lay_n_recv, lay_n_send){
        if (isempty(w_init) & lay_n_recv == 0 & lay_n_send == 0){
          return(F)
        }
        sum(dim(w_init) == c(lay_n_recv, lay_n_send)) != 2
      }

      w_dim_lay_dim <- mapply(check_w_lay_dim, w_init, lays_n_send, lays_n_recv)
      w_dim_lay_dim <- matrix(w_dim_lay_dim, nrow = nrow(w_init))

      if (sum(w_dim_lay_dim) > 0)
        stop(paste("Dimensions of weights and layers are inconsistent. Check the",
                   "following row(s) and column(s) ([row, column]) in inital",
                   "weight matrix and number of units in the corresponding layers.",
                   "\n"),
             which_matrix(w_dim_lay_dim))

      ## Setting the inital weights for each layer
      lay_inp_n <- apply(cxn, 1, function(x) sum(c(0, lays_n[x > 0])))

      # cbind weights row-wise (row receives inputs from columns) to have only
      # one weight matrix for every layer, (of course only for non-empty w_init
      # elements)
      wts <- apply(w_init, 1, function(x) Reduce("cbind", x))
      self$lays <- Map(function(x, y) {if(!isempty(y)) x$wt <- y; x}, net_lays, wts)

      # set the contrast-enhanced version of the weights
      self$lays <- lapply(self$lays, function(x) x$set_ce_weights(off, gain))
      private$cxn <- cxn
      private$n_lays <- n_lays
      private$n_units <- net_n_units
      private$off <- off
      private$gain <- gain
      private$gi <- gi
      private$cxn_b <- cxn_b
      private$lays_n_send <- lays_n_send
      private$lays_n_recv <- lays_n_recv
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
    cycle = function(ext_inputs = "list", clamp_inp){
      # Testing the arguments and reshaping the input
      if (!is.list(ext_inputs)){
        stop(paste("First argument to cycle function should be a list of",
                   "matrices or vectors."))
      }
      if (length(ext_inputs) != private$n_lays){
        stop(paste("Number of layers inconsistent with number of inputs in",
                   "network cycle. If a layer should not have an input, just use NULL",
                   "for that layer."))
      }

      # reshaping ext_inputs into vectors
      ext_inputs <- lapply(ext_inputs, c)

      # First we set all clamped layers to their input values
      if (clamp_inp){
        self$lays <- Map(function(x, y) ifelse(is.null(y), return(x),
                                                  return(x$clamp_cycle(y))),
                            self$lays, ext_inputs)
        lays_clamped <- !sapply(ext_inputs, is.null)
      }

      # We make a copy of the scaled activity for all layers
      scaled_acts <- lapply(self$lays, function(x) x$get_scaled_acts())

      ## For each unclamped layer, we put all its scaled ext_inputs in one
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
        ifelse(lay_clamped, return(lay), return(lay$cycle(intern_input,
                                                       ext_input)))
      }
      self$lays <- Map(cycle_non_clamp, self$lays, intern_input, ext_inputs,
                          lays_clamped)
      invisible(self)
    },

    #' chg_wt changes the weights of the entire network with the XCAL learning
    #' equation
    #'
    #' @param obj
    #' @rdname network
    chg_wt = function(){
      # updating the long-term averages
      self$lays <- lapply(self$lays, function(x) x$updt_avg_l())

      ## Extracting the averages for all layers
      avgs <- lapply(self$lays, function(x) x$get_act_avgs(private$m_in_s, private$avg_l_lrn_min, private$avg_l_lrn_max))
      avg_l <- lapply(avgs, function(x) x$avg_l)
      avg_m <- lapply(avgs, function(x) x$avg_m)
      avg_s_eff <- lapply(avgs, function(x) x$avg_s_eff)
      avg_l_lrn <- lapply(avgs, function(x) x$avg_l_lrn)

      ## For each connection matrix, calculate the intermediate vars.
      # make cxn-format versions of above averages
      # we need cxn_b with no connection specified with NULL
      cxn_b_null <- apply(private$cxn_b, c(1, 2),
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

      avg_s_eff_snd <- get_snd(avg_s_eff, cxn_b_null)
      avg_s_eff_rcv <- get_rcv(avg_s_eff, cxn_b_null)
      avg_m_snd <- get_snd(avg_m, cxn_b_null)
      avg_m_rcv <- get_rcv(avg_m, cxn_b_null)

      avg_l_snd <- get_snd(avg_l, cxn_b_null)
      avg_l_rcv <- get_rcv(avg_l, cxn_b_null)
      avg_l_lrn_rcv <- get_rcv(avg_l_lrn, cxn_b_null)

      s_hebb <- m_mapply(function(x, y) x %o% y, avg_s_eff_rcv, avg_s_eff_snd)
      m_hebb <- m_mapply(function(x, y) x %o% y, avg_m_rcv, avg_m_snd)

      # this is independent of sending; only the receiving neuron's long term
      # average is important; so just repeat the vector for every incoming
      # neuron
      l_recv <- m_mapply(function(x, y) matrix(rep(x, y), ncol = y), avg_l_rcv,
                         private$lays_n_recv)
      avg_l_lrn_rcv <- m_mapply(function(x, y) matrix(rep(x, y), ncol = y),
                                avg_l_lrn_rcv, private$lays_n_recv)

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
      self$lays <- lapply(self$lays, function(x) x$set_ce_weights(private$off, private$gain))
      invisible(self)
    },

    #' updt_pct_act updates the acts_p_avg and pct_act_scale variables for all
    #' layers. These variables update at the end of plus phases instead of cycle by
    #' cycle.
    #'
    #' @rdname network
    updt_pct_act = function(){
      # updates the acts_p_avg and pct_act_scale variables for all lays.
      # These variables update at the end of plus phases instead of
      # cycle by cycle. The avg_l values are not updated here.
      # This version assumes full connectivity when updating
      # pct_act_scale. If partial connectivity were to be used, this
      # should have the calculation in WtScaleSpec::SLayActScale, in
      # LeabraConSpec.cpp
      self$lays <- lapply(self$lays, function(x) x$updt_pct_act())
      invisible(self)
    },

    #' get_m1 obtains the m1 factor: the slope of the left-hand line in the
    #' "check mark" XCAL function.
    #'
    #' @rdname network
    get_m1 = function(){
      # obtains the m1 factor: the slope of the left-hand line in the
      # "check mark" XCAL function. Notice it includes the negative
      # sign.
      m <- (private$d_rev - 1) / private$d_rev
      return(m)
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
      f[idx1] <- self$get_m1() * x[idx1]
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
      self$lays <- lapply(self$lays, function(x) x$reset(random))
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
      which_matrix <- function(x){
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
             which_matrix(cxn_no_w))
      }
      # test 2
      w_no_cxn <- private$cxn == 0 & !w_empty
      if (sum(w_no_cxn) > 0)
        stop(paste("Non-Connected layers have weight matrix. Check the",
                   "following row(s) and column(s) ([row, column]) in initial ",
                   "connection and new weight matrix: \n"),
             which_matrix(w_no_cxn))

      # test 3
      # checks whether receiving and sending number of units in weight matrix
      # is correct for every layer
      check_w_lay_dim <- function(w_init, lay_n_recv, lay_n_send){
        if (isempty(w_init) & lay_n_recv == 0 & lay_n_send == 0){
          return(F)
        }
        sum(dim(w_init) == c(lay_n_recv, lay_n_send)) != 2
      }

      w_dim_lay_dim <- mapply(check_w_lay_dim, w_init, lays_n_send, lays_n_recv)
      w_dim_lay_dim <- matrix(w_dim_lay_dim, nrow = nrow(w_init))

      if (sum(w_dim_lay_dim) > 0)
        stop(paste("Dimensions of weights and layers are inconsistent. Check the",
                   "following row(s) and column(s) ([row, column]) in new",
                   "weight matrix and number of units in the corresponding layers.",
                   "\n"),
             which_matrix(w_dim_lay_dim))

      ## Now we set the weights
      # first find how many units project to the layer in all the network
      lay_inp_n <- apply(private$lays_n_recv, 1, sum)

      # make one weight matrix for every layer by collapsing receiving layer
      # weights columnwise, only do this for layers that receive something at all
      # (isempty)
      wts <- apply(w, 1, function(x) Reduce("cbind", x))
      self$lays <- Map(function(x, y) {if (!isempty(y)) x$wt <- y; x},
                          self$lays, wts)

      # set the contrast-enhanced version of the weights
      self$lays <- lapply(self$lays, set_ce_weights, private$off, private$gain)
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

  private = list(

    dim_lays = list(),
    cxn = matrix(),
    w_init = matrix(),
    n_lays = NULL,  # number of lays (number of objs in "lays")
    n_units = NULL, # total number of units in the network

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
    cxn_b = matrix(), # binary version of cxn
    lays_n_recv = matrix(), # number of units of receiving layer in cxn matrix
    lays_n_send = matrix(), # number of units of sending layer in cxn matrix
    w_index_low = matrix(), # lower index to extract weights from layer matrix
    # in cxn format
    w_index_up = matrix() # upper index to extract weights from layer matrix in
    # cxn format
  )
)
