#' @include unit.R

# leabra "cxn" class----------------------------------------------------------
#' Leabra connection class
#'
#' @field units a list with all the unit objects of the layer
cxn <- R6::R6Class("cxn",
  #public ----------------------------------------------------------------------
  public = list(
    # constructor
    initialize = function(cxn_table){
        self$wt_table <- plyr::ddply(cxn_table, 1, function(x) self$create_wt(x))
      },

    create_wt = function(cxn_table_row){
      # todo: add random -----
      d <- data.frame(strength = cxn_table_row$strength,
                      lay_send = cxn_table_row$lay_send,
                      lay_recv = cxn_table_row$lay_recv,
                      expand.grid(u_id_send = seq(cxn_table_row$n_s),
                                  u_id_recv = seq(cxn_table_row$n_r)))

      wt_table <- dplyr::mutate(d, wt = 0.3 + 0.4 * runif(nrow(d)),
                                ce_wt = private$set_ce_wt(wt,
                                                    off = private$off,
                                                    gain = private$gain))
      return(wt_table)
    },

#' chg_wt changes the weights of the entire network with the XCAL learning
#' equation
#'
#' @rdname cxn
chg_wt = function(avgs){
  # updating the long-term averages
  #lapply(self$lays, function(x) x$updt_unit_avg_l()) # this has to be done by
  #network

  ## Extracting the averages for all layers
  # avg_l <- lapply(avgs, function(x) x$avg_l)
  # avg_m <- lapply(avgs, function(x) x$avg_m)
  # avg_s_with_m <- lapply(avgs, function(x) x$avg_s_with_m)
  # avg_l_lrn <- lapply(avgs, function(x) x$avg_l_lrn)

  # combine the avgs with the cxn_table?
  # once with layer rcv, once with layer snd

  hebb_s <- avg_s_with_m_rcv * avg_s_with_m_snd
  hebb_m <- avg_m_rcv * avg_m_snd

  dwt_m <- m_mapply(function(x, y) self$get_dwt(x, y), hebb_s, hebb_m)
  dwt_l <- m_mapply(function(x, y) self$get_dwt(x, y), hebb_s, l_recv)

  #  todo: rename m_lrn ---------
  dwt <- self$lrate * (dwt_m * private$m_lrn + avg_l_lrn * dwt_l)

  self$bound_weights_exp(dwt)
  self$set_ce_wt()
  invisible(self)
},

# fields -----------------------------------------------------------------------
  cxn = NULL, #ce_wt =  # contrast enhanced version of weights
  wt_table = NULL # contains all wt values with enough information to calculate weight
# changes
),
# private ----------------------------------------------------------------------
private = list(
  #' get_dwt
  #'
  #' get_dwt gets the weight changes with the "check mark" function in XCAL
  #'
  #' @param hebb_s activation vector of sending and receiving units (hebbian
  #'   co-product)
  #'
  #' @param thres vector of threshold values with the same size as hebb_s, this
  #'   is either the hebbian medium product, or the receiving unit's longterm
  #'   average; it is the point at which weight changes start to get negative in
  #'   the xcal function
  #'
  #' @rdname cxn
  get_dwt = function(hebb_s, thres){
    f <- rep(0, length(hebb_s))
    xcal_neg_slope_thres <- private$xcal_neg_slope_constant * thres
    xcal_neg_slope <- (private$xcal_rev - 1) / private$xcal_rev
    idx1 <- (hebb_s > private$xcal_min_hebb_s) & (hebb_s < xcal_neg_slope_thres)
    f[idx1] <- hebb_s[idx1] * xcal_neg_slope
    idx2 <- hebb_s >= xcal_neg_slope_thres
    f[idx2] <- hebb_s[idx2] - thres[idx2]
    return(f)
  },
  #' set_ce_weights
  #'
  #' set_ce_weights sets contrast enhanced weight values
  #'
  #' @rdname cxn
             set_ce_wt = function(wt, off, gain){
               1 / (1 + (off * (1 - wt) / wt) ^ gain)
             },

  #' bound_weights_exp bounds weights exponentially; weight increases
  #' exponentially toward upper bound of 1, and decreases toward lower bound of
  #' 0. based on linear, non-contrast enhanced weights.
  #'
  #' @rdname cxn
             bound_weights_exp = function(dwt){
               self$wt <- ifelse(dwt > 0, wt + (1 - wt) * dwt, wt + wt * dwt)
               invisible(self)
             },
  # fields -----------#
  xcal_neg_slope_constant = 0.1,
  xcal_min_hebb_s = 0.0001,
  off = 1,
  gain = 6
)
)
