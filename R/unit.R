#' @include misc.R
NULL

#' Class to simulate a biologically realistic neuron
#'
#' This class simulates a biologically realistic neuron in the lebra framework. A layer has the field "units", which is a list of unit objects.
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @keywords data
#' @return Object of \code{\link{R6Class}} with methods for calculating neuron activity changes
#' @format \code{\link{R6Class}} object.
#'
#' @examples
#' unit$new()
#' unit$cycle()
#'
#' @field act percentage activation ("firing rate") of the unit, which is sent
#'   to other units, think of it as a percentage of how many neurons are active
#'   in a microcolumn of 100 neurons
#' @field g_e excitatory conductance, asymptotically approaches g_e_raw (see
#'   cycle method), aka net input
#' @field v membrane potential in the range between 0 and 2, 0 is -100mV, 2 is
#'   100mV, thus 1 is 0mV (v should usually be between 0 and 1)
#' @field v_eq equilibrium membrane potential, which is not reset by spikes, it
#'   just keeps integrating
#' @field spike a flag that indicates spiking threshold was crossed
#' @field i_adapt adaptation current like in the AdEx model
#' @field avg_ss super short-term running average activation, integrates over
#'   act
#' @field avg_s short-term running average activation, integrates over avg_ss,
#'   represents plus phase learning signal
#' @field avg_m medium-term running average activation, integrates over avg_s,
#'   represents minus phase learning signal
#' @field avg_l long-term running average activation, integrates over avg_m,
#'   drives long-term floating average for self-organized learning learning
#'@field avg_l_prc this is avg_l in percentage of minimum (default is 0.1) and maximum (default is 1.5)
#'
#' @section Methods:
#' \describe{
#'   \item{Documentation}{For full documentation of each method go to https://github.com/lightning-viz/lightining-r/}
#'   \item{\code{new()}}{Creates an object of this class with default parameters.}
#'
#'   \item{\code{cycle(g_e_raw, g_i)}}{Cycles 1 ms with given excitatory conductance  \code{g_e_raw} and inhibitory conductance \code{g_i}. Excitatory conductance depends on the weights to other units and the activity of those other units. Inhibitory conductance depends on feedforward and feedback inhibition. See layer cycle method.}
#'
#'   \item{\code{clamp_cycle(activation)}}{Clamps the value of \code{activation} to the \code{act} variable of the unit without any time integration. Then updates averages. This is usually done when presenting external input.}
#'
#'   \item{\code{updt_avg_l()}}{This method updates the variable \code{avg_l}. This usually happens before the weights are changed in the network (after the plus phase), and not every cycle. It tends to move towards a constant minium avg_l (0.1) if avg_m is smaller than 0.2; if it is larger, than it will tend to go to avg_m}
#'
#'   \item{\code{nxx1(x)}}{Calculates the activation of a unit as a function of v or g_e and their thresholds}}

unit <- R6::R6Class("unit",
  # public ---------------------------------------------------------------------
  public = list(
    cycle = function(g_e_raw, g_i){
      ## updating g_e input
      self$g_e <- self$g_e + private$cyc_dt * private$g_e_dt *
        (g_e_raw - self$g_e)

      ## Finding membrane potential
      # excitatory, inhibitory and leak current
      i_e <- self$g_e * (private$v_rev_e - self$v)
      i_i <- g_i * (private$v_rev_i - self$v)
      i_l <- private$g_l * (private$v_rev_l - self$v)
      i_net <- i_e + i_i + i_l

      # almost half-step method for updating v (i_adapt doesn't half step)
      v_h <- self$v + 0.5 * private$cyc_dt * private$v_dt *
        (i_net - self$i_adapt)
      i_e_h <- self$g_e * (private$v_rev_e - v_h)
      i_i_h <- g_i * (private$v_rev_i - v_h)
      i_l_h <- private$g_l * (private$v_rev_l - v_h)
      i_net_h <- i_e_h + i_i_h + i_l_h

      self$v <- self$v + private$cyc_dt * private$v_dt *
        (i_net_h - self$i_adapt)
      self$v_eq <- self$v_eq + private$cyc_dt * private$v_dt *
        (i_net_h - self$i_adapt)

      ## Finding activation
      # finding threshold excitatory conductance
      g_e_thr <- (g_i * (private$v_rev_i - private$v_thr) +
                    private$g_l * (private$v_rev_l - private$v_thr) -
                    self$i_adapt) / (private$v_thr - private$v_rev_e)

      # finding whether there's an action potential
      if (self$v > private$spk_thr){
        self$spike <- 1
        self$v <- private$v_reset
      }  else {
        self$spike <- 0
      }

      # finding instantaneous rate due to input
      if (self$v_eq <= private$v_thr){
        new_act <- private$nxx1(self$v_eq - private$v_thr)
      } else {
        new_act <- private$nxx1(self$g_e - g_e_thr)
      }

      # update activity
      self$act <- self$act + private$cyc_dt * private$v_dt *
        (new_act - self$act)

      ## Updating adaptation current
      self$i_adapt <- self$i_adapt + private$cyc_dt *
        (private$i_adapt_dt * (private$v_gain * (self$v - private$v_rev_l)
                            - self$i_adapt) + self$spike * private$spike_gain_adapt)

      private$update_averages()
      invisible(self)
    },

    clamp_cycle = function(activation){
      ## Clamping the activty to the activation
      self$act <- activation
      private$update_averages()
      invisible(self)
    },

    updt_avg_l = function(){
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
      self$avg_ss <- self$act
      self$avg_s <- self$act
      self$avg_m <- self$act
      self$avg_l <- self$act
      self$g_e <- 0
      self$v <- 0.3
      self$v_eq <- 0.3
      self$i_adapt <- 0
      self$spike <- 0
      invisible(self)
    },
    # fields -------------------------------------------------------------------
    act = 0.2,
    avg_s = 0.2,
    avg_m = 0.2,
    avg_l = 0.2,
    avg_ss = 0.2,
    g_e = 0,
    v = 0.3,
    v_eq = 0.3,
    i_adapt = 0,
    spike = 0
  ),
  # active ---------------------------------------------------------------------
  active = list(
    avg_l_prc = function() {
      (self$avg_l - private$avg_l_min) / (private$avg_l_max - private$avg_l_min)
    }
  ),
  # private --------------------------------------------------------------------
  private = list(
    nxx1 = function(x){
      # nxx1_df is a df that is used as a lookup-table, it is stored internally
      # but you can generate the data with the create_nxx1 function
      approx(nxx1_df$nxx1_dom, nxx1_df$nxoxp1, x, method = "constant",
             rule = 2)$y
    },

    update_averages = function(){
    self$avg_ss <- self$avg_ss + private$cyc_dt * private$ss_dt *
      (self$act - self$avg_ss)
    self$avg_s <- self$avg_s + private$cyc_dt * private$s_dt *
      (self$avg_ss - self$avg_s)
    self$avg_m <- self$avg_m + private$cyc_dt * private$m_dt *
      (self$avg_s - self$avg_m)
    invisible(self)
    },
  # fields ---------------------------------------------------------------------
    # dynamic values
  # constant values
  # time step constants
  g_e_dt = 1 / 1.4,     # time step constant for update of "g_e"
  cyc_dt = 1,           # time step constant for integration of cycle dynamics
  v_dt = 1 / 3.3,       # time step constant for membrane potential
  i_adapt_dt = 1 / 144, # time step constant for adaptation
  ss_dt = 0.5,          # time step constant for super-short average
  s_dt = 0.5,           # time step constant for short average
  m_dt = 0.1,           # time step constant for medium-term average
  avg_l_dt = 0.1,       # time step constant for long-term average
  # other
  avg_l_max = 1.5,      # max value of avg_l
  avg_l_min = 0.1,      # min value of avg_l
  v_rev_e = 1,          # excitatory reversal potential
  v_rev_i = .25,        # inhibitory reversal potential
  v_rev_l = 0.3,        # leak reversal potential
  g_l = 0.1,            # leak conductance
  v_thr = 0.5,          # normalized "rate threshold", corresponds with -50mV
  spk_thr = 1.2,        # normalized spike threshold
  v_reset= 0.3,         # reset membrane potential after spike
  v_gain = 0.04,        # gain that voltage produces on adaptation
  spike_gain_adapt = 0.00805 # effect of spikes on adaptation
  )
)
