#' @include unit.R
NULL

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
#' @return Object of \code{\link{R6Class}} with methods for calculating changes
#'   of activity in a layer of neurons
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
#' l <- layer$new(c(5, 5))
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
#' @field units a list with all the \link{unit} objects of the layer
#' @field wt a receiving x sending weight matrix, where the receiving unit (row)
#'   has the current weight values for all sending units (columns)
#' @field ce_wt sigmoidal contrast-enhanced version of the weight matrix wt
#' @field n number of units in layer
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new(dim, g_i_gain = 2)}}{Creates an object of this class with
#'   default parameters.
#'
#'     \code{dim} A a pair of numbers giving the rows and columns of the layer.
#'
#'     \code{g_i_gain} Gain factor for inhibitory conductance. If you want less
#'     activation in a layer, set this higher.}
#'
#'   \item{\code{get_unit_acts()}}{Returns a vector with the activities of all
#'   units of a layer.}
#'
#'   \item{\code{get_unit_scaled_acts()}}{get_unit_scaled_acts returns a vector
#'   with the scaled activities of all units of a layer. Scaling is done with
#'   recip_avg_act_n, a reciprocal function of the number of active units.}
#'
#'   \item{\code{cycle(intern_input, ext_input)}}{Iterates one time step with
#'   layer object.
#'
#'    \code{intern_input} Vector with inputs from all other layers. Each input
#'   has already been scaled by a reciprocal function of the number of active
#'   units (recip_avg_act_n) of the sending layer and by the connection strength
#'   between the receiving and sending layer. The weight matrix is multiplied
#'   with this input vector to get the excitatory conductance for each unit in
#'   the layer.
#'
#'     \code{ext_input} Vector with inputs not coming from another layer, with
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
#'     \code{acts} acts you want to clamp to the units in the layer}
#'
#'   \item{\code{get_unit_act_avgs()}}{Returns a list with the short, medium and
#' long term activation averages of all units in the layer as vectors. The super
#' short term average is not returned, and the long term average is not updated
#' before being returned. These averages are used by the network class to
#' calculate weight changes.}
#'
#'   \item{\code{updt_unit_avg_l()}}{Updates the long-term average (avg_l) of all
#' units in the layer, usually done after a plus phase.}
#'
#'   \item{\code{updt_recip_avg_act_n()}}{Updates the avg_act_inert and
#' recip_avg_act_n variables, these variables update at the end of plus phases
#' instead of cycle by cycle. This version assumes full connectivity when
#' updating recip_avg_act_n.}
#'
#'   \item{\code{reset(random = F)}}{Sets the activation and activation averages
#' of all units to 0. Used to begin trials from a stationary point.
#'
#'     \code{random} Logical variable, if TRUE the activation ist set randomly
#' between .05 and .95 for every unit instead of 0.}
#'
#'   \item{\code{set_ce_weights()}}{Sets contrast enhanced weight values}
#'
#'   \item{\code{get_unit_vars(show_dynamics = T, show_constants = F)}}{Returns
#'   a data frame with with the current state of all unit variables in the
#'   layer. Every row is a unit. You can choose whether you want dynamic values
#'   and / or constant values. This might be useful if you want to analyse what
#'   happens in units of a layer, which would otherwise not be possible, because
#'   most of the variables (fields) are private in unit class.}
#'
#'   \item{\code{get_layer_vars(show_dynamics = T, show_constants = F)}}{Returns a
#' data frame with 1 row with the current state of the variables in the layer.
#' You can choose whether you want dynamic values and / or constant values. This
#' might be useful if you want to analyse what happens in a layer, which would
#' otherwise not be possible, because some of the variables (fields) are private
#' in the layer class.}}
layer <- R6::R6Class("layer",
  #public ----------------------------------------------------------------------
  public = list(
    # constructor
    initialize = function(dim, g_i_gain = 2){
      if (length(c(dim)) == 2){
        self$n <- prod(dim)
        unit1 <- unit$new()
        # use cloning
        self$units <- lapply(seq(self$n), function(x) unit1$clone(deep = T))
      } else{
        stop("dim argument should be of the type c(rows, columns)")
      }

      Map(function(x, y) x$unit_number <- y, self$units, 1:self$n)
      # recip_avg_act_n is initialized with number of units that are > 0.4
      # (n_act_lrgr_forty)
      private$avg_act_inert <- self$avg_act
      n_act_lrgr_forty <- max(c(sum(self$get_unit_acts() > 0.4), 1))
      private$recip_avg_act_n <- 1 / (n_act_lrgr_forty + 2)
      private$g_fbi <- private$g_fbi_gain * self$avg_act
      private$g_ffi <- private$g_ffi_gain * max(c(private$g_e_avg -
                                                    private$g_ffi_thr, 0))
      invisible(self)
    },

    get_unit_acts = function(){
      sapply(self$units, function(x) x$act)
    },

    get_unit_scaled_acts = function(){
      private$recip_avg_act_n * self$get_unit_acts()
    },

    cycle = function(intern_input, ext_input = NULL){
      # obtaining the excitatory conductance because of input
      # contrast enhanced weights are used
      g_e_per_unit <- self$ce_wt %*% intern_input
      if (!isempty(ext_input)) g_e_per_unit <- g_e_per_unit + ext_input

      # obtaining inhibitory conductance
      private$g_e_avg <- mean(g_e_per_unit)
      private$g_ffi <- private$g_ffi_gain * max(c(private$g_e_avg -
                                                    private$g_ffi_thr, 0))
      private$g_fbi <- private$g_fbi + private$g_fbi_dt *
        (private$g_fbi_gain * self$avg_act - private$g_fbi)
      g_i <- private$g_i_gain * (private$g_ffi + private$g_fbi)

      # calling the cycle method for all units
      Map(function(x, y, z) x$cycle(y, z), self$units, g_e_per_unit, g_i)
      invisible(self)
    },

    clamp_cycle = function(acts){
      Map(function(x, y) x$clamp_cycle(y), self$units, acts)
      # updating inhibition for the next cycle
      private$g_fbi <- private$g_fbi_gain * self$avg_act
      invisible(self)
    },

    get_unit_act_avgs = function(m_avg_prc_in_s_avg, avg_l_lrn_min,
                                 avg_l_lrn_max){
      avg_s <- sapply(self$units, function(x) x$avg_s)
      avg_m <- sapply(self$units, function(x) x$avg_m)
      avg_l <- sapply(self$units, function(x) x$avg_l)

      # obatining avg_s_with_m
      avg_s_with_m <- m_avg_prc_in_s_avg * avg_m +
        (1 - m_avg_prc_in_s_avg) * avg_s

      # obtaining avg_l_lrn
      avg_l_lrn <- avg_l_lrn_min + private$get_unit_avg_l_prc() *
        (avg_l_lrn_max - avg_l_lrn_min)

      list(
        "avg_s" = avg_s,
        "avg_m" = avg_m,
        "avg_l" = avg_l,
        "avg_s_with_m" = avg_s_with_m,
        "avg_l_lrn" = avg_l_lrn
      )
    },

    updt_unit_avg_l = function(){
      lapply(self$units, function(x) x$updt_avg_l())
      invisible(self)
    },

    updt_recip_avg_act_n = function(){
      private$avg_act_inert <- private$avg_act_inert +
        private$avg_act_inert_dt * (self$avg_act - private$avg_act_inert)
      n_units_avg_act <- max(round(private$avg_act_inert * private$n), 1)
      private$recip_avg_act_n <- 1 / (n_units_avg_act + 2)
      invisible(self)
    },

    reset = function(random = F){
      lapply(self$units, function(x) x$reset(random = random))
      invisible(self)
    },

    set_ce_weights = function(off, gain){
      self$ce_wt <- 1 / (1 + (off * (1 - self$wt) / self$wt) ^ gain)
      invisible(self)
    },

    get_unit_vars = function(show_dynamics = T, show_constants = F){
      unit_vars_list <- lapply(self$units, function(x)
        x$get_vars(show_dynamics, show_constants))
      unit_vars_df <- do.call(rbind, unit_vars_list)
    },

    get_layer_vars = function(show_dynamics = T, show_constants = F){
      df <- data.frame(layer = self$layer_number)
      dynamic_vars <- data.frame(
        avg_act = self$avg_act,
        avg_act_inert = private$avg_act_inert,
        g_e_avg = private$g_e_avg,
        g_fbi = private$g_fbi,
        g_ffi = private$g_ffi,
        recip_avg_act_n = private$recip_avg_act_n
      )
      constant_vars <-
        data.frame(
          n = self$n,
          g_ffi_gain           = private$g_ffi_gain,
          g_ffi_thr            = private$g_ffi_thr,
          g_fbi_gain           = private$g_fbi_gain,
          g_fbi_dt             = private$g_fbi_dt,
          g_i_gain             = private$g_i_gain,
          avg_act_inert_dt     = private$avg_act_inert_dt
        )
      if (show_dynamics == T) df <- cbind(df, dynamic_vars)
      if (show_constants == T) df <- cbind(df, constant_vars)
      return(df)
    },

    # fields -------------------------------------------------------------------
    n = NULL,          # number of units
    # An n x i weights matrix, where the n-th row has the
    # current wt values for all inputs coming to the n-th unit
    wt = NULL,
    ce_wt = NULL,      # contrast-enhanced version of wt. ce_wt = SIG(wt).
    units = NULL,      # a list with all the unit objects of the layer
    layer_number = 1   # number of layer in the network
  ),

  # private --------------------------------------------------------------------
  private = list(
    # get_unit_avg_l_prc returns the percentage value of avg_l in the
    # range of the minimum and maximum values for avg_l specified with the
    # parameters
    get_unit_avg_l_prc = function(){
      sapply(self$units, function(x) x$avg_l_prc)
    },

    # fields -------------------------------------------------------------------
    # dynamic ------------------------------------------------------------------

    # recip_avg_act_n is the scaling factor for the outputs coming out from THIS
    # layer. Notice this is different from C++ version. Only updated with
    # updt_recip_avg_act_n.
    recip_avg_act_n = NULL,
    avg_act_inert = NULL,

    g_e_avg = 0,       # average g_e for all units during last cycle
    g_fbi = NULL,      # feedback inhibition
    g_ffi = NULL,      # feedforward inhibition

    # constant -----------------------------------------------------------------
    g_ffi_gain = 1,          # gain for feedforward inhibition
    g_ffi_thr = 0.1,         # threshold for feedforward inhibition
    g_fbi_gain = 0.5,        # gain for feedback inhibition
    g_fbi_dt = 1 / 1.4,      # time step for fb inhibition (fb_tau = 1.4)
    g_i_gain = 2,            # overall gain on inhibition
    avg_act_inert_dt = 0.01  # time step constant for updating avg_act_inert
  ),
  # active ---------------------------------------------------------------------
  active = list(
    # dependent
    #
    # avg_act returns the average activity of all units
    #
    avg_act = function(){
      mean(self$get_unit_acts())
    })
)
