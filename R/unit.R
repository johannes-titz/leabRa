#' @include misc.R
NULL

#' Class to simulate a biologically realistic neuron
#'
#' This class simulates a biologically realistic neuron in the lebra framework.
#' When you use a layer class, you will see that a layer object has a variable
#' (field) \code{units}, which is a list of unit objects.
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @keywords data
#' @return Object of \code{\link{R6Class}} with methods for calculating neuron
#'   activity changes
#' @format \code{\link{R6Class}} object.
#'
#' @examples
#' u <- unit$new() # creates a new unit with default leabra values
#'
#' print(u) # a lot of private values
#' u$v # private values cannot be accessed
#' # if you want to see alle variables, you need to use the function
#' u$get_vars(show_dynamics = T, show_constants = T)
#'
#' # let us clamp the activation to 0.7
#' u$act
#' u$clamp_cycle(0.7)
#' c(u$act, u$avg_s, u$avg_m, u$avg_l, u$avg_l_prc)
#' # act is indeed 0.7, but avg_l was not updated, this only happens before the
#' # weights are changed, let us update it now
#' u$updt_avg_l()
#' c(u$act, u$avg_s, u$avg_m, u$avg_l, u$avg_l_prc)
#' # seems to work
#'
#' # let us run 10 cycles with unclamped activation and output the activation
#' # produced because of changes in conductance
#' u <- unit$new()
#' cycle_number <- 1:10
#' result <- lapply(cycle_number, function(x)
#'                  u$cycle(g_e_raw = 0.5, g_i = 0.5)$get_vars())
#' # make a data frame out of the list
#' result <- plyr::ldply(result)
#' # plot act
#' plot(result$act, type = "b", xlab = "cycle", ylab = "act")
#' # add conductance g_e to plot, should approach g_e_raw
#' lines(result$g_e, type = "b", col = "blue")
#'
#' @field act percentage activation ("firing rate") of the unit, which is sent
#'   to other units, think of it as a percentage of how many neurons are active
#'   in a microcolumn of 100 neurons
#' @field avg_s short-term running average activation, integrates over avg_ss (a
#'   private variable, which integrates over act), represents plus phase
#'   learning signal
#' @field avg_m medium-term running average activation, integrates over avg_s,
#'   represents minus phase learning signal
#' @field avg_l long-term running average activation, integrates over avg_m,
#'   drives long-term floating average for self-organized learning
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new()}}{Creates an object of this class with default
#'   parameters.}
#'
#'   \item{\code{cycle(g_e_raw, g_i)}}{Cycles 1 ms with given excitatory
#'   conductance  \code{g_e_raw} and inhibitory conductance \code{g_i}.
#'   Excitatory conductance depends on the weights to other units and the
#'   activity of those other units. Inhibitory conductance depends on
#'   feedforward and feedback inhibition. See layer cycle method.}
#'
#'   \item{\code{clamp_cycle(activation)}}{Clamps the value of \code{activation}
#'   to the \code{act} variable of the unit without any time integration. Then
#'   updates averages. This is usually done when presenting external input.}
#'
#'   \item{\code{updt_avg_l()}}{Updates the variable \code{avg_l}. This usually
#'   happens before the weights are changed in the network (after the plus
#'   phase), and not every cycle. If avg_m is smaller than 0.2 (or equal) avg_l
#'   tends to move towards the constant minium avg_l value (0.1). If avg_m is
#'   larger than 0.2, than avg_l will tend to go to the constant maximum avg_l
#'   value (1.5)}
#'
#'   \item{\code{get_vars(show_dynamics = T, show_constants = F)}}{Returns a
#'   data frame with 1 row with the current state of all the variables of the
#'   unit. You can choose whether you want dynamic values and / or constant
#'   values. This might be useful if you want to analyse what happens in a unit,
#'   which would otherwise not be possible, because most of the variables
#'   (fields) are private in this class.}}
#'
unit <- R6::R6Class("unit",
  # public ---------------------------------------------------------------------
  public = list(
    initialize = function(){
      private$nxx1_df <- create_nxx1()
    },
    cycle = function(g_e_raw, g_i){
      ## updating g_e input
      private$g_e <- private$g_e + private$cyc_dt * private$g_e_dt *
        (g_e_raw - private$g_e)

      ## Finding membrane potential
      # excitatory, inhibitory and leak current
      i_e <- private$g_e * (private$v_rev_e - private$v)
      i_i <- g_i * (private$v_rev_i - private$v)
      i_l <- private$g_l * (private$v_rev_l - private$v)
      i_net <- i_e + i_i + i_l

      # almost half-step method for updating v (i_adapt doesn't half step)
      v_h <- private$v + 0.5 * private$cyc_dt * private$v_dt *
        (i_net - private$i_adapt)
      i_e_h <- private$g_e * (private$v_rev_e - v_h)
      i_i_h <- g_i * (private$v_rev_i - v_h)
      i_l_h <- private$g_l * (private$v_rev_l - v_h)
      i_net_h <- i_e_h + i_i_h + i_l_h

      private$v <- private$v + private$cyc_dt * private$v_dt *
        (i_net_h - private$i_adapt)
      private$v_eq <- private$v_eq + private$cyc_dt * private$v_dt *
        (i_net_h - private$i_adapt)

      ## Finding activation
      # finding threshold excitatory conductance
      g_e_thr <- (g_i * (private$v_rev_i - private$v_thr) +
                    private$g_l * (private$v_rev_l - private$v_thr) -
                    private$i_adapt) / (private$v_thr - private$v_rev_e)

      # finding whether there's an action potential
      if (private$v > private$spk_thr){
        private$spike <- 1
        private$v <- private$v_reset
      }  else {
        private$spike <- 0
      }

      # finding instantaneous rate due to input
      if (private$v_eq <= private$v_thr){
        new_act <- private$nxx1(private$v_eq - private$v_thr)
      } else {
        new_act <- private$nxx1(private$g_e - g_e_thr)
      }

      # update activity
      self$act <- self$act + private$cyc_dt * private$v_dt *
        (new_act - self$act)

      ## Updating adaptation current
      private$i_adapt <- private$i_adapt + private$cyc_dt *
        (private$i_adapt_dt * (private$v_gain * (private$v - private$v_rev_l)
                               - private$i_adapt)
         + private$spike * private$spike_gain_adapt)

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
          self$avg_l <- self$avg_l + private$l_dt *
            (private$avg_l_max - self$avg_m)
        } else{
          self$avg_l <- self$avg_l + private$l_dt *
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
      private$v <- 0.3
      private$v_eq <- 0.3
      private$i_adapt <- 0
      private$spike <- 0
      invisible(self)
    },
    # return a data frame with all variables
    get_vars = function(show_dynamics = T, show_constants = F){
      df <- data.frame(unit = self$unit_number)
      dynamic_vars <- data.frame(
        act = self$act,
        avg_ss = private$avg_ss,
        avg_s = self$avg_s,
        avg_m = self$avg_m,
        avg_l = self$avg_l,
        avg_l_prc = self$avg_l_prc,
        g_e = private$g_e,
        v = private$v,
        v_eq = private$v_eq,
        i_adapt = private$i_adapt,
        spike = private$spike
      )
      constant_vars <-
        data.frame(
          cyc_dt           = private$cyc_dt,
          g_e_dt           = private$g_e_dt,
          ss_dt            = private$ss_dt,
          s_dt             = private$s_dt,
          l_dt             = private$l_dt,
          m_dt             = private$m_dt,
          v_dt             = private$v_dt,
          i_adapt_dt       = private$i_adapt_dt,
          avg_l_max        = private$avg_l_max,
          avg_l_min        = private$avg_l_min,
          v_rev_e          = private$v_rev_e,
          v_rev_l          = private$v_rev_l,
          g_l              = private$g_l,
          v_thr            = private$v_thr,
          spk_thr          = private$spk_thr,
          v_reset          = private$v_reset,
          v_gain           = private$v_gain,
          spike_gain_adapt = private$spike_gain_adapt
        )
      if (show_dynamics == T) df <- cbind(df, dynamic_vars)
      if (show_constants == T) df <- cbind(df, constant_vars)
      return(df)
    },
    # fields -------------------------------------------------------------------
    act = 0.2,
    avg_s = 0.2,
    avg_m = 0.2,
    avg_l = 0.2,
    # number of unit in the layer, if you create a single unit object, this is
    # one, otherwise the layer will set this value
    unit_number = 1
  ),
  # active ---------------------------------------------------------------------
  active = list(
    avg_l_prc = function() {
      (self$avg_l - private$avg_l_min) / (private$avg_l_max - private$avg_l_min)
    }
  ),
  # private --------------------------------------------------------------------
  # \item{\code{nxx1(x)}}{Calculates the activation of a unit as a function of v or g_e and their thresholds}}
  private = list(
    nxx1 = function(x){
      # nxx1_df is a df that is used as a lookup-table, it is stored internally
      # but you can generate the data with the create_nxx1 function
      closest_value <- which(abs(private$nxx1_df$nxx1_dom-x) ==
                               min(abs(private$nxx1_df$nxx1_dom-x)))
      private$nxx1_df$nxoxp1[closest_value]
      #approx(private$nxx1_df$nxx1_dom, private$nxx1_df$nxoxp1, x,
      #       method = "constant", rule = 2)$y
    },

    create_nxx1 = function(){
      # we don't have precalculated vectors for interpolation
      n_x <- 2000 # size of the precalculated vectors
      mid <- 2 # mid length of the domain
      domain <- seq(-mid, mid, length.out = n_x) # will be "nxx1_dom"
      # domain of Gaussian
      dom_g <- seq(-2 * mid, 2 * mid, length.out = 2 * n_x)
      values <- rep(0, n_x) # will be "nxoxp1"
      sd <- .005 # standard deviation of the Gaussian
      gaussian <- exp(- (dom_g ^ 2) / (2 * sd ^ 2)) / (sd * sqrt(2 * pi))

      XX1 <- function(x, gain = 100){
        x[x <= 0] <- 0
        x[x > 0] <- gain * x[x > 0] / (gain * x[x > 0] + 1) # gain = 100 default
        return(x)
      }

      for (p in 1:n_x){
        low <- n_x - p + 1
        high <- 2 * n_x - p
        values[p] <- sum(XX1(domain) * gaussian[low:high])
        values[p] <- values[p] / sum(gaussian[low:high])
      }
      nxx1_dom <- domain
      nxoxp1 <- values
      as.data.frame(cbind(nxoxp1, nxx1_dom))
    },

    update_averages = function(){
    private$avg_ss <- private$avg_ss + private$cyc_dt * private$ss_dt *
      (self$act - private$avg_ss)
    self$avg_s <- self$avg_s + private$cyc_dt * private$s_dt *
      (private$avg_ss - self$avg_s)
    self$avg_m <- self$avg_m + private$cyc_dt * private$m_dt *
      (self$avg_s - self$avg_m)
    invisible(self)
    },
  # fields ---------------------------------------------------------------------
  # dynamic values
  avg_ss = 0.2,
  g_e = 0,
  v = 0.3,
  v_eq = 0.3,
  i_adapt = 0,
  spike = 0,
  # constant values
  # time step constants
  g_e_dt = 1 / 1.4,     # time step constant for update of "g_e"
  cyc_dt = 1,           # time step constant for integration of cycle dynamics
  v_dt = 1 / 3.3,       # time step constant for membrane potential
  i_adapt_dt = 1 / 144, # time step constant for adaptation
  ss_dt = 0.5,          # time step constant for super-short average
  s_dt = 0.5,           # time step constant for short average
  m_dt = 0.1,           # time step constant for medium-term average
  l_dt = 0.1,       # time step constant for long-term average
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
  spike_gain_adapt = 0.00805, # effect of spikes on adaptation
  nxx1_df = NULL        # this is the nxx1 activation function as a data frame
  )
)
