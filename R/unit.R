#' Leabra unit class
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @return Object of \code{\link{R6Class}}
#' @format \code{\link{R6Class}} object.
#' @field act activation ("firing rate") of the unit (sent to other units)
#' @field g_e time-averaged excitatory conductance, asymptotically approaches
#' g_e_raw (see cycle method), aka net input
#' @field v_m membrane potential
#' @field v_m_eq equilibrium membrane potential -- not reset by spikes -- just
#' keeps integrating
#' @field spike a flag that indicates spiking threshold was crossed
#' @field adapt adaptation current
#' @field avg_ss super short-term running average activation
#' @field avg_s short-term running average activation, integrates over avg_ss,
#' represents plus phase learning signal
#' @field avg_m medium-term running average activation, integrates over avg_s,
#' represents minus phase learning signal
#' @field avg_l long-term running average activation, integrates over avg_m,
#' drives long-term floating average for Hebbian learning
#'
#' @section Methods:
#' \describe{
#'   \item{\code{example_method(parameter_1 = 3)}}{This method uses \code{parameter_1} to ...}
#'
#'   \item{\code{nxx1(x = 3)}}{nxx1 calculates noisy x/(x+1) function
#' (convolution of x/(x+1) with a Gaussian function),
#'
#' x is a numeric vector to calculate nxx1 (x is either v_m_eq - v_m_thr or g_e - g_e_thr)}
#' }
unit <- R6::R6Class("unit",
  # public ---------------------------------------------------------------------
  public = list(
    initialize = function(){

    },
    cycle = function(g_e_raw, gc_i){
      ## updating g_e input
      private$g_e <- private$g_e + private$cyc_dt * private$g_e_dt *
        (g_e_raw - private$g_e)

      ## Finding membrane potential
      # excitatory, inhibitory and leak current
      i_e <- private$g_e * (private$e_rev_e - private$v_m)
      i_i <- gc_i * (private$e_rev_i - private$v_m)
      i_l <- private$gc_l * (private$e_rev_l - private$v_m)
      i_net <- i_e + i_i + i_l

      # almost half-step method for updating v_m (adapt doesn't half step)
      v_m_h <- private$v_m + 0.5 * private$cyc_dt * private$v_m_dt *
        (i_net - private$adapt)
      i_e_h <- private$g_e * (private$e_rev_e - v_m_h)
      i_i_h <- gc_i * (private$e_rev_i - v_m_h)
      i_l_h <- private$gc_l * (private$e_rev_l - v_m_h)
      i_net_h <- i_e_h + i_i_h + i_l_h

      private$v_m <- private$v_m + private$cyc_dt * private$v_m_dt *
        (i_net_h - private$adapt)
      private$v_m_eq <- private$v_m_eq + private$cyc_dt * private$v_m_dt *
        (i_net_h - private$adapt)

      ## Finding activation
      # finding threshold excitatory conductance
      g_e_thr <- (gc_i * (private$e_rev_i - private$v_m_thr) +
                    private$gc_l * (private$e_rev_l - private$v_m_thr) -
                    private$adapt) / (private$v_m_thr - private$e_rev_e)


      # finding whether there's an action potential
      cat("private$v_m = ", private$v_m, "\n")
      if (private$v_m > private$spk_thr){
        private$spike <- 1
        private$v_m <- private$v_m_r
      }  else {
        private$spike <- 0
      }

      # finding instantaneous rate due to input
      if (private$v_m_eq <= private$v_m_thr){
        new_act <- private$nxx1(private$v_m_eq - private$v_m_thr)
      } else {
        new_act <- private$nxx1(private$g_e - g_e_thr)
      }

      # update activity
      self$act <- self$act + private$cyc_dt * private$v_m_dt *
        (new_act - self$act)

      ## Updating adaptation
      private$adapt <- private$adapt + private$cyc_dt *
        (private$adpt_dt * (private$v_m_gain * (private$v_m - private$e_rev_l)
                            - private$adapt) + private$spike * private$spike_gain)

      ## updating averages
      private$avg_ss <- private$avg_ss + private$cyc_dt * private$ss_dt *
        (self$act - private$avg_ss)
      self$avg_s <- self$avg_s + private$cyc_dt * private$s_dt *
        (private$avg_ss - self$avg_s)
      self$avg_m <- self$avg_m + private$cyc_dt * private$m_dt *
        (self$avg_s - self$avg_m)
      invisible(self)
    },
    clamp_cycle = function(activation){
      ## Clamping the activty to the activation
      self$act <- activation

      ## updating averages
      private$avg_ss <- private$avg_ss + private$cyc_dt * private$ss_dt *
        (self$act - private$avg_ss)
      self$avg_s <- self$avg_s + private$cyc_dt * private$s_dt *
        (private$avg_ss - self$avg_s)
      self$avg_m <- self$avg_m + private$cyc_dt * private$m_dt *
        (self$avg_s - self$avg_m)
      invisible(self)
    },
    updt_avg_l = function(){
        # This fuction updates the long-term average "avg_l"
        # u = this unit
        # Based on the description in:
        # https://grey.colorado.edu/ccnlab/index.php/#Leabra_Hog_Prob_Fix#Adaptive_Contrast_Impl
        # it tends to move towards the minium avg_l (0.1) if avg_m is smaller
        # 0.2 if it is larger, than it will tend to got to avg_m
        if (self$avg_m > 0.2){
          self$avg_l <- self$avg_l + private$avg_l_dt *
            (private$avg_l_max - self$avg_m)
        } else{
          self$avg_l <- self$avg_l + private$avg_l_dt *
            (private$avg_l_min - self$avg_m)
        }
        invisible(self)
      },
    reset = function(random = F){
      ifelse(random == T, self$act <- 0.05 + 0.9 * runif(1),
             self$act <- 0)
      # does this influence pct_act_rel?
      private$avg_ss <- self$act
      self$avg_s <- self$act
      self$avg_m <- self$act
      self$avg_l <- self$act
      private$g_e <- 0
      private$v_m <- 0.3
      private$v_m_eq <- 0.3
      private$adapt <- 0
      private$spike <- 0
      invisible(self)
    },
    # fields -------------------------------------------------------------------
    act = 0.2,
    avg_s = 0.2,
    avg_m = 0.2,
    avg_l = 0.2
  ),
  # active ---------------------------------------------------------------------
  active = list(
               avg_l_prc = function(){
                 (self$avg_l - private$avg_l_min) /
                   (private$avg_l_max - private$avg_l_min)
               }),
  # private --------------------------------------------------------------------
  private = list(

    #' @rdname unit
              nxx1 = function(x){
                # nxx1_df is a df that is used as a lookup-table, it is stored internally
                # but you can generate the data with the create_nxx1 function
                approx(nxx1_df$nxx1_dom, nxx1_df$nxoxp1, x, method = "constant",
                       rule = 2)$y
              },
  # fields ---------------------------------------------------------------------
    # dynamic values
    g_e = 0,
    avg_ss = 0.2,
    v_m = 0.3,
    v_m_eq = 0.3,
    adapt = 0,
    spike = 0,
    # constant values
    g_e_dt = 1 / 1.4,   # time step constant for update of "g_e"
    cyc_dt = 1,       # time step constant for integration of cycle dynamics
    v_m_dt = 1 / 3.3,   # time step constant for membrane potential
    l_dn_dt = 1 / 2.5,  # time step constant for avg_l decrease
    adpt_dt = 1 / 144, # time step constant for adaptation
    ss_dt = 0.5,        # time step for super-short average
    s_dt = 0.5,         # time step for short average
    m_dt = 0.1,         # time step for medium-term average
    avg_l_dt = 0.1,     # time step for long-term average
    avg_l_max = 1.5,    # max value of avg_l; why can this be larger than 1?----
    avg_l_min = 0.1,    # min value of avg_l
    e_rev_e = 1,        # excitatory reversal potential
    e_rev_i = .25,      # inhibitory reversal potential
    e_rev_l = 0.3,      # leak reversal potential
    gc_l = 0.1,         # leak conductance
    v_m_thr = 0.5,      # normalized "rate threshold", -50mV (0: -100mV,
    # 2: 100mV)
    spk_thr = 1.2,      # normalized spike threshold
    v_m_r = 0.3,        # reset membrane potential after spike
    v_m_gain = 0.04,    # gain that voltage produces on adaptation
    spike_gain = 0.00805, # effect of spikes on adaptation
    l_up_inc = 0.2      # increase in avg_l if avg_m has been "large"
  )
)
