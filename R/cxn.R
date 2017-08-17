
cxn <- R6::R6Class("cxn",
  #public ----------------------------------------------------------------------
  public = list(
    # constructor
    initialize = function(cxn_table, dim_lays){
        dim_lays <- mutate(dim_lays, n_units = n_rows * n_cols)
        private$cxn_table <- cxn_table
        private$dim_lays <- dim_lays
        self$weight_matrix_per_layer <- private$create_wt_tbl()
      },
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
    #
    create_wt_tbl = function(){
      private$cxn_table <- dplyr::group_by(private$cxn_table, lay_recv)
      private$normalize_rows()
      private$combine_cxn_table_with_dim_lays()

      weights <- plyr::ddply(private$cxn_table, 1,
                             function(x) private$create_wt_for_one_cxn(x))
      # weight_matrix_of_net <- reshape2::acast(weights,
      #                                         unique_id_r~unique_id_s,
      #                                         value.var = "wt")
      private$weight_matrix_per_layer <- plyr::dlply(weights, "layer_r",
                                             private$create_weights_from_ids)
    },

    #
    normalize_rows = function(){
      private$cxn_table <- dplyr::mutate(private$cxn_table,
                                       strength = strength / sum(strength))
    },

    #
    combine_cxn_table_with_dim_lays = function(){
      private$cxn_table <- dplyr::left_join(private$cxn_table,
                                          private$dim_lays,
                                          by = c("lay_send" = "layer"))
      private$cxn_table <- dplyr::left_join(private$cxn_table,
                                          private$dim_lays,
                                          by = c("lay_recv" = "layer"),
                                          suffix = c("_s", "_r"))
      invisible(self)
    },

    create_wt_for_one_cxn = function(cxn_table_row){
      # todo: add random -----
      d <- data.frame(strength = cxn_table_row$strength,
                      layer_s = cxn_table_row$lay_send,
                      layer_r = cxn_table_row$lay_recv,
                      expand.grid(unit_s = seq(cxn_table_row$n_units_s),
                                  unit_r = seq(cxn_table_row$n_units_r)))
      # create unique names for sending and receiving
      wt_table <- dplyr::mutate(d, wt = 0.3 + 0.4 * runif(nrow(d)),
                                unique_id_s = paste("layer", layer_s, "unit", unit_s, sep=""),
                                unique_id_r = paste("layer", layer_r, "unit", unit_r, sep=""))#,
      #ce_wt = set_ce_weights(wt)),
      #
      wt_table <- mutate(wt_table,
                         unique_id_s = private$char_to_factor_keep_order(unique_id_s),
                         unique_id_r = private$char_to_factor_keep_order(unique_id_r))
      return(wt_table)
    },

    # todo -- the other way araound----------
    create_weights_from_ids = function(weights){
      reshape2::acast(weights, unique_id_r~unique_id_s, value.var = "wt", fill = 0)
    },

    char_to_factor_keep_order = function(char){
      factor(char, levels = unique(char))
    },

    # get_dwt
    #
    # get_dwt gets the weight changes with the "check mark" function in XCAL
    #
    # @param hebb_s activation vector of sending and receiving units (hebbian
    #   co-product)
    #
    # @param thres vector of threshold values with the same size as hebb_s, this
    #   is either the hebbian medium product, or the receiving unit's longterm
    #   average; it is the point at which weight changes start to get negative in
    #   the xcal function
    #
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
    # set_ce_weights
    #
    # set_ce_weights sets contrast enhanced weight values
    #
    # @rdname cxn
               set_ce_wt = function(wt, off, gain){
                 1 / (1 + (off * (1 - wt) / wt) ^ gain)
               },

    # bound_weights_exp bounds weights exponentially; weight increases
    # exponentially toward upper bound of 1, and decreases toward lower bound of
    # 0. based on linear, non-contrast enhanced weights.
    #
    # @rdname cxn
               bound_weights_exp = function(dwt){
                 self$wt <- ifelse(dwt > 0, wt + (1 - wt) * dwt, wt + wt * dwt)
                 invisible(self)
               },


    # fields -----------#
    xcal_neg_slope_constant = 0.1,
    xcal_min_hebb_s = 0.0001,
    off = 1,
    gain = 6,
    cxn_table = NULL,
    dim_lays = NULL,
    weight_matrix_per_layer = NULL
  )
)
