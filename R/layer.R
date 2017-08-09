#' @include unit.R

# leabra "layer" class----------------------------------------------------------
#' Leabra layer class
#'
#' @field units a list with all the unit objects of the layer
#' @field recip_avg_act_n scaling factor for the outputs coming out from
#' this layer, notice this is different from C++ version, only updated with
#' updt_pct_act
#' @field avg_act_inert time-averaged version of avg_act, it is
#' updated only at the end of plus phases, so it is not a proper time
#' average, see updt_pct_act.
#' @field g_e_avg average excitatory conductance (net input) for all
#' units during last cycle
#' @field wt an NxI weights matrix, where the n-th row has the current wt values
#' for all inputs coming to the n-th unit
#' @field ce_wt contrast-enhanced version of wt matrix, ce_wt = SIG(wt)
#' @field n number of units
#' @field g_fbi feedback inhibition
#' @field g_i_gain overall gain on inhibition
layer <- R6::R6Class("layer",
  #public ----------------------------------------------------------------------
  public = list(
    # constructor
    initialize = function(dims, g_i_gain = 2){
      if (length(c(dims)) == 2){
        self$n <- prod(dims)
        unit1 <- unit$new()
        # use cloning
        self$units <- lapply(seq(self$n), function(x) unit1$clone(deep = T))
      } else{
        stop("dims argument should be of the type c(rows, columns)")
      }
      # recip_avg_act_n is initialized with number of units that are > 0.4
      # (n_act_lrgr_forty)
      private$avg_act_inert <- self$avg_act
      n_act_lrgr_forty <- max(c(sum(self$get_unit_acts() > 0.4), 1))
      private$recip_avg_act_n <- 1 / (n_act_lrgr_forty + 2)
      private$g_fbi <- private$g_fbi_gain * self$avg_act
      invisible(self)
    },
    #' get_unit_acts returns a vector with the activities of all units of a
    #' layer
    #'
    #' @rdname layer
    get_unit_acts = function(){
      sapply(self$units, function(x) x$act)
    },
    #' get_unit_scaled_acts returns a vector with the scaled activities of all
    #' units of a layer, scaling is done with recip_avg_act_n, a reciprocal
    #' function of the number of active units
    #'
    #' @rdname layer
    get_unit_scaled_acts = function(){
      data.frame(act_scaled = private$recip_avg_act_n * self$get_unit_acts(),
                 u_id = 1:self$n)
    },
    #' cycle iterates one time step with layer object
    #'
    #' @param g_e_intern single combined vector with g_e from all layers.
    #' @param ext_input vector with inputs not coming from another layer, with
    #'   length equalling the number of units in this layer. If empty, no
    #'   external inputs are processed. Note that this is basically an
    #'   excitatory conductance value (g_e), although it is also used as an
    #'   activation value when layers are clamped
    #' @rdname layer
    cycle = function(g_e_intern, ext_input = NULL){
      ## obtaining the g_e_per_unit
      g_e_per_unit <- g_e_intern # contrast enhanced weights are used
      if (!isempty(ext_input)) g_e_per_unit <- g_e_per_unit + ext_input
      ## obtaining inhibition
      private$g_e_avg <- mean(g_e_per_unit)
      g_ffi <- private$g_ffi_gain * max(c(private$g_e_avg - private$g_ffi_thres, 0))
      private$g_fbi <- private$g_fbi + private$g_fbi_dt *
        (private$g_fbi_gain * self$avg_act - private$g_fbi)
      g_i <- private$g_i_gain * (g_ffi + private$g_fbi)
      ## calling the cycle method for all units
      # todo: split cycle function, so that nxx1 call can be done with a
      # functional, much faster!
      Map(function(x, y, z) x$cycle(y, z), self$units, g_e_per_unit, g_i)
      invisible(self)
    },

    #' clamp_cycle iterates one time step with layer object with clamped
    #' activations (i.e activations are instantenously set without dt)
    #'
    #' @param activations vector of activations that you want to clamp
    #'
    #' @rdname layer
    clamp_cycle = function(activations){
      self$units <- Map(function(x, y) x$clamp_cycle(y), self$units, activations)
      # updating inhibition for the next cycle
      private$g_fbi <- private$g_fbi_gain * self$avg_act
      invisible(self)
    },

    #' get_unit_act_avgs returns a list with the s,m,l averages in the layer as
    #' vectors, notice that the ss average is not returned, and avg_l is not
    #' updated before being returned
    #'
    #' @rdname layer
    get_unit_act_avgs = function(m_avg_prc_in_s_avg, avg_l_lrn_min,
                                 avg_l_lrn_max){
      avg_s <- sapply(self$units, function(x) x$avg_s)
      avg_m <- sapply(self$units, function(x) x$avg_m)
      avg_l <- sapply(self$units, function(x) x$avg_l)

      # obatining avg_s_with_m
      avg_s_with_m <- m_avg_prc_in_s_avg * avg_m +
        (1 - m_avg_prc_in_s_avg) * avg_s

      # obtaining avg_l_lrn
      # avg_l_lrn will be a percentage value between max and min
      avg_l_lrn <- avg_l_lrn_min + self$get_unit_avg_l_prc() *
        (avg_l_lrn_max - avg_l_lrn_min)

      data.frame(unit_id = 1:length(avg_s), "avg_s" = avg_s, "avg_m" = avg_m,
                 "avg_l" = avg_l, "avg_s_with_m" = avg_s_with_m,
                 "avg_l_lrn" = avg_l_lrn)
    },

    #' get_unit_avg_l_prc returns the percentage value of avg_l in the
    #' range of the minimum and maximum values for avg_l specified with the
    #' parameters
    #'
    #' @rdname layer
    #'
    # maybe move to private? ---------------------------------------------------
    get_unit_avg_l_prc = function(){
      sapply(self$units, function(x) x$avg_l_prc)
    },

    #' updt_unit_avg_l updates the long-term average (avg_l) of all the units in
    #' the layer, usually done after a plus phase
    #'
    #' @rdname layer
    updt_unit_avg_l = function(){
      self$units <- lapply(self$units, function(x) x$updt_avg_l())
      invisible(self)
    },

    #' updt_pct_act updates the avg_act_inert and recip_avg_act_n variables,
    #' these variables update at the end of plus phases instead of cycle by
    #' cycle. This version assumes full connectivity when updating
    #' recip_avg_act_n. If partial connectivity were to be used, this should
    #' have the calculation in WtScaleSpec::SLayActScale, in LeabraConSpec.cpp
    #'
    #' @rdname layer
    updt_pct_act = function(){
      private$avg_act_inert <- private$avg_act_inert + private$avg_act_inert_dt *
        (self$avg_act - private$avg_act_inert)
      n_units_avg_act <- max(round(private$avg_act_inert * private$n), 1)
      private$recip_avg_act_n <- 1 / (n_units_avg_act + 2)
      invisible(self)
    },

    #' reset sets the activation and activation averages of all units to 0, used
    #' to begin trials from a random stationary point, the activity values may
    #' also be set to random values between .05 and .95
    #'
    #' @param random logical variable, if TRUE set activation randomly between
    #'   .05 and .95, if FALSE set activation to 0
    #' @rdname layer
    reset = function(random = F){
      self$units <- lapply(self$units, function(x) x$reset(random = random))
      invisible(self)
    },

    #' set_ce_weights sets contrast enhanced weight values
    #'
    #' @rdname layer
    set_ce_weights = function(off, gain){
      self$ce_wt <- 1 / (1 + (off * (1 - self$wt) / self$wt) ^ gain)
      invisible(self)
    },
    # fields ----------------------
    n = NULL,          # number of units
    wt = NULL,         # An NxI weights matrix, where the n-th row has the
    # current wt values for all inputs coming to the n-th unit
    ce_wt = NULL,      # contrast-enhanced version of wt. ce_wt = SIG(wt).
    units = NULL        # an array with all the unit objs of the layer
  ),
  # private --------------------------------------------------------------------
  private = list(
    # dynamic

    # recip_avg_act_n is the scaling factor for the outputs coming out from THIS
    # layer. Notice this is different from C++ version. Only updated with
    # updt_pct_act.
    recip_avg_act_n = NULL,
    avg_act_inert = NULL,

    g_e_avg = 0,       # average g_e for all units during last cycle
    g_fbi = NULL,      # feedback inhibition

    # constant
    g_ffi_gain = 1,          # gain for feedforward inhibition
    g_ffi_thres = 0.1,       # threshold for feedforward inhibition
    g_fbi_gain = 0.5,        # gain for feedback inhibition
    g_fbi_dt = 1 / 1.4,      # time step for fb inhibition (fb_tau = 1.4)
    g_i_gain = 2,            # overall gain on inhibition
    avg_act_inert_dt = 0.01  # time step constant for updating avg_act_inert
  ),

  active = list(
    # dependent
    #
    #' avg_act returns the average activity of all units
    #'
    #' @rdname layer
    avg_act = function(){
      mean(self$get_unit_acts())
    })
)
