#' @include layer.R
#' @include cxn.R

# leabra network class----------------------------------------------------------
#' leabra network class
#'
#' @field dim_lays list of 2d vectors that contains number of rows and columns
#' of lays
#' @field cxn matrix specifying connection strength between lays, if
#' layer j sends projections to layer i, then cxn[i, j] = c > 0; 0
#' otherwise, c specifies the relative strength of that connection with respect
#' to the other projections to layer i
#' @field w_init matrix of initial weight matrices, this is analogous to cxn,
#' i.e. w_init[i, j] contains the initial weight matrix for the cxn from
#' layer j to i
#' @field gi vector of gi values for every layer, this comes in handy to control
#' overall level of inhibition of specific lays
#' @field off "offset" in the SIG function for contrast enhancement
#' @field gain gain in the SIG function for contrast enhancement
#' @field lrate learning rate
network <-  R6::R6Class("network",
  public = list(
    #' initialize
    #'
    #' constructor for network class
    #' @rdname network
    initialize = function(lays_tbl, cxn_tbl, off = 1, gain = 6,
                          lrate = 0.1){

      #if (isempty(lays_tbl$g_i)) cat("now g_i specified in layer")
      # todo: test if cxn and dim_lays are valid objects-----
      # cxn$strength > 0
      # number of units per layer
      lays_tbl <- dplyr::mutate(lays_tbl, n = n_rows * n_cols)

      cxn_tbl <- dplyr::group_by(cxn_tbl, lay_recv)
        # Normalizing the rows of "cxn" so they add to 1
      cxn_tbl <- dplyr::mutate(cxn_tbl,
                        strength = strength / sum(strength))
        # combine with lays_tbl
      cxn_tbl <- dplyr::left_join(cxn_tbl, lays_tbl,
                                  by = c("lay_send" = "layer"))
      cxn_tbl <- dplyr::left_join(cxn_tbl, lays_tbl,
                                  by = c("lay_recv" = "layer"),
                           suffix = c("_s", "_r"))

      # list that contains layer objects
      self$lays <- plyr::dlply(lays_tbl, 1, function(x)
        layer$new(c(x$n_rows, x$n_cols), x$g_i, x$layer))
      # set the contrast-enhanced version of the weights
      lapply(self$lays, function(x) x$set_ce_weights(off, gain))

      # Setting the inital weights for each layer
      self$cxn <- cxn$new(cxn_tbl)
      self$lays_tbl <- lays_tbl
      self$cxn_tbl <- cxn_tbl
      private$n_lays <- nrow(lays_tbl)
      private$n_units <- sum(lays_tbl$n)
      self$lrate <- lrate
      invisible(self)
      },
    #' cycle iterates one time step with network object
    #'
    #' @param ext_inputs a data frame with five variables: input, layer, type
    #' unit_id and activation
    #' @param clamp_inp a binary flag; 1: lays are clamped to their input value
    #' 0: inputs summed to netins.
    #' @rdname network
    cycle = function(ext_input, clamp_inp){
      # todo: tests -------------
      # Testing the arguments and reshaping the input
      # test whether any layers or units are used that do not exist
      # test whether activation is between 0 and 1
      # test whether any hidden unit does have activations or any input / output
      # unit does not have activations
      # test whether number of inputs per layer correspond with number of units
      # in layer
      ext_inputs <- split(ext_inputs$activation, ext_inputs$layer)

      # set all layers to unclamped; if clamp_inp is true, then set all layers
      # to unclamped that have an intern_input
      lays_clamped = dplyr::transmute(self$lays_tbl, layer = layer, clamped = F)

      # First we set all clamped layers to their input values
      if (clamp_inp){
        Map(function(x, y) if (!all(is.na(y))) x$clamp_cycle(y),
            self$lays, ext_inputs)
        # if mean activation is NA, then there was no input, so the layer is
        # not clamped and intern_input has to be calculated
        lays_clamped <- ext_input %>%
          group_by(layer) %>%
          summarize(clamped = !is.na(mean(activation)))
      }

      # We make a copy of the scaled activity for all layers
      scaled_acts <- plyr::ldply(self$lays, function(x) x$get_unit_scaled_acts())

      # find intern_input for unclamped layers; you will need scaled_acts and
      # the cxn_strength
      #
      # actually we could multiply it with the weight already, right?
      wt_table <- dplyr::left_join(self$cxn$wt_table, scaled_acts,
                       by = c("lay_send" = "layer", "u_id_send" = "u_id"))
      # scale activity by strength as well
      wt_table <- dplyr::mutate(wt_table, g_e = act_scaled * strength * ce_wt)
      wt_table <- dplyr::group_by(wt_table, lay_recv, u_id_recv)
      wt_table <- dplyr::summarize(wt_table, g_e_intern = sum(g_e))
      g_e_intern <- split(wt_table$g_e_intern, wt_table$lay_recv)

      # for unclamped layers
      # if the layer is not clamped, run cycle with inputs
      cycle_non_clamp <- function(lay, g_e_intern, ext_input, lay_clamped){
        if (lays_clamped$clamped) lay$cycle(g_e_intern, ext_input)
      }
      self$lays <- Map(cycle_non_clamp, self$lays, g_e_intern, ext_inputs,
                          lays_clamped)
      invisible(self)
    },

    #' chg_wt changes the weights of the entire network with the XCAL learning
    #' equation
    #'
    #' @rdname network
    chg_wt = function(){
      # updating the long-term averages
      self$lays <- lapply(self$lays, function(x) x$updt_unit_avg_l())
      # this has to stay, but then we can call cxn...

      ## Extracting the averages for all layers
      avgs <- lapply(self$lays, function(x) x$get_act_avgs(private$m_avg_prc_in_s_avg, private$avg_l_lrn_min, private$avg_l_lrn_max))
      avgs <- plyr::ldply(avgs)
      # avg_l <- lapply(avgs, function(x) x$avg_l)
      # avg_m <- lapply(avgs, function(x) x$avg_m)
      # avg_s_with_m <- lapply(avgs, function(x) x$avg_s_with_m)
      # avg_l_lrn <- lapply(avgs, function(x) x$avg_l_lrn)

      # now can we call the cxn functions?
      self$cxn$chg_wt(avgs, ...)

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
      self$lays <- lapply(self$lays, function(x) x$set_ce_weights(self$off, self$gain))
      invisible(self)
    },

    #' updt_pct_act updates the avg_act_inert and recip_avg_act_n variables for
    #' all layers. These variables update at the end of plus phases instead of
    #' cycle by cycle.
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

    #' get_dwt get weight changes with the "check mark" function in XCAL
    #'
    #' @param x vector of abscissa values
    #' @param th vector of threshold values with the same size as x
    #' @rdname network
    get_dwt = function(x, th){
      f <- rep(0, length(x))
      idx1 <- (x > private$d_thr) & (x < private$xcal_rev * th)
      idx2 <- x >= private$xcal_rev * th
      f[idx1] <- self$m1 * x[idx1]
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
      self$lays <- lapply(self$lays, set_ce_weights, self$off, self$gain)
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
    # fields -------------------------------------------------------------------
    lrate = 0.1,  # learning rate for XCAL
    lays = list(),
    cxn = NULL,
    lays_tbl = NULL,
    cxn_tbl = NULL,
    off = 1, # "offset" in the SIG function for contrast enhancement
    gain = 6 # gain in the SIG function for contrast enhancement
    ),
  # active ---------------------------------------------------------------------
  active = list(
    #' get_m1 obtains the m1 factor: the slope of the left-hand line in the
    #' "check mark" XCAL function.
    #'
    #' @rdname network
    m1 = function(){
      # obtains the m1 factor: the slope of the left-hand line in the "check
      # mark" XCAL function. Notice it includes the negative sign.
      (private$xcal_rev - 1) / private$xcal_rev
    }
  ),
# private -------------
# fields ----------
  private = list(
    dim_lays = list(),
    #cxn = matrix(),
    w_init = matrix(),
    n_lays = NULL,  # number of lays (number of objs in "lays")
    n_units = NULL, # total number of units in the network

    # constants
    avg_l_lrn_max = 0.01,# max amount of "BCM" learning in XCAL
    avg_l_lrn_min = 0.0, # min amount of "BCM" learning in XCAL
    m_avg_prc_in_s_avg = 0.1, # proportion of medium to short term avgs. in XCAL
    m_lrn = 1, # proportion of error-driven learning in XCAL
    d_thr = 0.0001, # threshold for XCAL "check mark" function
    xcal_rev = 0.1, # reversal value for XCAL "check mark" function

    gi = 2, # gi for layers, to control overall inhibition in a specific layer
    avg_l_lrn = list(),
    #dependent
    #m1 = NULL, # the slope in the left part of XCAL's "check mark"
    cxn_b = matrix(), # binary version of cxn
    lays_n_recv = matrix(), # number of units of receiving layer in cxn matrix
    lays_n_send = matrix(), # number of units of sending layer in cxn matrix
    w_index_low = matrix(), # lower index to extract weights from layer matrix
    # in cxn format
    w_index_up = matrix() # upper index to extract weights from layer matrix in
    # cxn format
  )
)
