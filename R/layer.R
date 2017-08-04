#' @include unit.R

# leabra "layer" class----------------------------------------------------------
#' Leabra layer class
#'
#' @field units a list with all the unit objects of the layer
#' @field recip_avg_act_n scaling factor for the outputs coming out from
#' this layer, notice this is different from C++ version, only updated with
#' updt_pct_act
#' @field acts_p_avg time-averaged version of acts_avg, it is
#' updated only at the end of plus phases, so it is not a proper time
#' average, see updt_pct_act.
#' @field g_e_avg average excitatory conductance (net input) for all
#' units during last cycle
#' @field wt an NxI weights matrix, where the n-th row has the current wt values
#' for all inputs coming to the n-th unit
#' @field ce_wt contrast-enhanced version of wt matrix, ce_wt = SIG(wt)
#' @field n number of units
#' @field fbi feedback inhibition
#' @field gi overall gain on inhibition
layer <- R6::R6Class("layer",
  public = list(
    # constructor
    initialize = function(dims, gi = 2){
      if (length(c(dims)) == 2){
        self$n <- prod(dims)
        unit1 <- unit$new()
        # use cloning
        self$units <- lapply(seq(self$n), function(x) unit1$clone(deep = T))
        # notice how the unit array is 1-D. The 2-D structure of the
        # layer doesn't have meaning in this part of the code
      } else{
        stop("dims argument should be of the type c(rows, columns)")
      }
      # recip_avg_act_n is initialized as a reciprocal function of number of
      # units that are > 0.4 (avg_act_n)
      private$acts_p_avg <- self$get_acts()
      avg_act_n <- max(c(sum(self$get_acts() > 0.4), 1))
      private$recip_avg_act_n <- 1 / (avg_act_n + 2)
      private$fbi <- private$fb * self$acts_avg
      invisible(self)
    },
    #' get_acts returns a vector with the activities of all units of a layer
    #'
    #' @rdname layer
    get_acts = function(){
      sapply(self$units, function(x) x$act)
    },
    #' get_scaled_act returns a vector with the scaled activities of all units
    #' of a layer, scaling is done with recip_avg_act_n, it changes only
    #'
    #' @rdname layer
    get_scaled_acts = function(){
      private$recip_avg_act_n * self$get_acts()
    },
    #' cycle iterates one time step with layer object
    #'
    #' @param intern_input An Ix1 matrix, where I is the total number of inputs from
    #' all layers. Each input has already been scaled by the recip_avg_act_n of its
    #' layer of origin and by the wt_scale_rel factor.
    #' @param ext_input An Nx1 matrix denoting inputs that don't come from another
    #' layer, where N is the number of units in this layer. An empty matrix
    #' indicates that there are no external inputs.
    #' @rdname layer
    cycle = function(intern_input, ext_input = 0){
      ## obtaining the g_es
      g_es <- self$ce_wt %*% intern_input  # you use contrast-enhanced weights
      if (!isempty(ext_input)) g_es <- g_es + ext_input
      ## obtaining inhibition
      private$g_e_avg <- mean(g_es)
      ffi <- private$ff * max(c(private$g_e_avg - private$ff0, 0))
      private$fbi <- private$fbi + private$fb_dt * (private$fb * self$get_avg_act() - private$fbi)
      gc_i <- private$gi * (ffi + private$fbi)
      ## calling the cycle method for all units

      Map(function(x, y, z) x$cycle(y, z), self$units, g_es, gc_i)
      invisible(self)
    },

    #' clamp_cycle iterates one time step with layer object with clamped input (i.e
    #' activation is instantenously set to input)
    #'
    #' @param input vector of activations that you want to clamp
    #' @rdname layer
    clamp_cycle = function(input){
      self$units <- Map(function(x, y) x$clamp_cycle(y), self$units, input)
      # updating inhibition for the next cycle
      private$fbi <- private$fb * self$get_avg_act()
      invisible(self)
    },

    #' get_act_avgs returns a list with the s,m,l averages in the layer as vectors,
    #' notice that the ss average is not returned, and avg_l is not updated before
    #' being returned
    #'
    #' @rdname layer
    get_act_avgs = function(m_in_s, avg_l_lrn_min, avg_l_lrn_max){
      avg_s <- sapply(self$units, function(x) x$avg_s)
      avg_m <- sapply(self$units, function(x) x$avg_m)
      avg_l <- sapply(self$units, function(x) x$avg_l)

      # obatining avg_s_eff
      avg_s_eff <- m_in_s * avg_m + (1 - m_in_s) * avg_s

      # obtaining avg_l_lrn
      avg_l_lrn <- avg_l_lrn_min + self$get_rel_avg_l() *
        (avg_l_lrn_max - avg_l_lrn_min)

      list("avg_s" =  avg_s, "avg_m" = avg_m, "avg_l" = avg_l,
           "avg_s_eff" = avg_s_eff, "avg_l_lrn" = avg_l_lrn)
    },

    #' get_rel_avg_l returns the relative values of avg_l, these are the dependent
    #' variables rel_avg_l in all units used in latest XCAL
    #'
    #' @rdname layer
    get_rel_avg_l = function(){
      l_avg_rel <- sapply(self$units, function(x) x$rel_avg_l)
    },

    #' updt_avg_l updates the long-term average (avg_l) of all the units in the
    #' layer, usually done after a plus phase
    #'
    #' @rdname layer
    updt_avg_l = function(){

      self$units <- lapply(self$units, function(x) x$updt_avg_l())
      invisible(self)
    },

    #' updt_pct_act updates the acts_p_avg and recip_avg_act_n variables, these
    #' variables update at the end of plus phases instead of cycle by cycle. This
    #' version assumes full connectivity when updating recip_avg_act_n. If partial
    #' connectivity were to be used, this should have the calculation in
    #' WtScaleSpec::SLayActScale, in LeabraConSpec.cpp
    #'
    #' @rdname layer
    updt_pct_act = function(){
      private$acts_p_avg <- private$acts_p_avg + private$acts_p_avg_dt *
        (self$get_avg_act() - private$acts_p_avg)
      r_avg_act_n <- max(round(private$acts_p_avg * private$n), 1)
      private$recip_avg_act_n <- 1 / (r_avg_act_n + 2)
      invisible(self)
    },

    #' reset sets the activation of all units to a random value, and sets all
    #' activation time averages to that value, used to begin trials from a random
    #' stationary point, the activity values may also be set to zero
    #'
    #' @param random logical variable, if TRUE set activation randomly between .05
    #' and .95, if FALSE set activation to 0
    #' @rdname layer
    reset = function(random = F){
      self$units <- lapply(self$units, function(x) x$reset(random = random))
      invisible(self)
    },

    #' get_avg_act gets the value of acts_avg, the mean of unit activities for
    #' the layer
    #'
    #' @rdname layer
    get_avg_act = function(){
      # get the value of acts_avg, the mean of unit activities
      avg_acts <- mean(self$get_acts())
    },
    #' set_ce_weights sets contrast enhanced weight values
    #'
    #' @rdname layer
    set_ce_weights = function(off, gain){
      self$ce_wt <- 1 / (1 + (off * (1 - self$wt) / self$wt) ^ gain)
      invisible(self)
    },
    n = NULL,          # number of units
    wt = NULL,         # An NxI weights matrix, where the n-th row has the
    # current wt values for all inputs coming to the n-th unit
    ce_wt = NULL,      # contrast-enhanced version of wt. ce_wt = SIG(wt).
    units = NULL        # an array with all the unit objs of the layer
  ),
  private = list(
    # dynamic

    recip_avg_act_n = NULL,  # The scaling factor for the outputs coming out
    # from THIS layer. Notice this is different from C++
    # version. Only updated with updt_pct_act.
    acts_p_avg = NULL,

    g_e_avg = 0,       # average net input for all units during last cycle
    fbi = NULL,        # feedback inhibition

    # constant
    ff = 1,            # gain for feedforward inhibition
    ff0 = 0.1,         # threshold for feedforward inhibition
    fb = 0.5,          # gain for feedback inhibition
    fb_dt = 1/1.4,     # time step for fb inhibition (fb_tau = 1.4)
    gi = 2,            # overall gain on inhibition
    acts_p_avg_dt = 0.01  # time step constant for updating acts_p_avg
  ),

  active = list(
    # dependent
    # average activity of all units
    #' acts_avg returns a vector with the activities of all units of a layer
    #'
    #' @rdname layer
    acts_avg = function(){
      mean(self$get_acts())
    })
)
