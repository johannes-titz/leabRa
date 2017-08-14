#' @include misc.R
NULL

#' Class to simulate a biologically realistic neuron
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @keywords data
#' @return Object of \code{\link{R6Class}} with methods for calculating neuron activity changes. A layer consists of multiple unit objects.
#' @format \code{\link{R6Class}} object.
#'
#' @examples
#' Lightning$new("http://localhost:3000/")
#' Lightning$new("http://your-lightning.herokuapp.com/")
#'
#' @field act percentage activation ("firing rate") of the unit, which is sent
#'   to other units, think of it as a percentage of how many neurons are active
#'   in a microcolumn of 100 neurons
#' @field g_e excitatory conductance, asymptotically approaches g_e_raw (see cycle method), aka net input
#' @field v_m membrane potential in the range between 0 and 2, 0 is -100mV, 2 is 100mV, thus 1 is 0mV (v_m should usually be between 0 and 1)
#' @field v_m_eq equilibrium membrane potential, which is not reset by spikes, it just
#' keeps integrating
#' @field spike a flag that indicates spiking threshold was crossed
#' @field i_adapt adaptation current
#' @field avg_ss super short-term running average activation, integrates over act
#' @field avg_s short-term running average activation, integrates over avg_ss,
#' represents plus phase learning signal
#' @field avg_m medium-term running average activation, integrates over avg_s,
#' represents minus phase learning signal
#' @field avg_l long-term running average activation, integrates over avg_m,
#' drives long-term floating average for self-organized learning learning
#'
#' @section Methods:
#' \describe{
#'   \item{Documentation}{For full documentation of each method go to https://github.com/lightning-viz/lightining-r/}
#'   \item{\code{new(serveraddress)}}{This method is used to create object of this class with \code{serveraddress} as address of the server object is connecting to.}
#'
#'   \item{\code{sethost(serveraddress)}}{This method changes server that you are contacting with to \code{serveraddress}.}
#'   \item{\code{createsession(sessionname = "")}}{This method creates new session on the server with optionally given name in \code{sessionname}.}
#'   \item{\code{usesession(sessionid)}}{This method changes currently used session on the server to the one with id given in \code{sessionid} parameter.}
#'   \item{\code{openviz(vizid = NA)}}{This method by default opens most recently created by this object visualization. If \code{vizid} parameter is given, it opens a visualization with given id instead.}
#'   \item{\code{enableautoopening()}}{This method enables auto opening of every visualisation that you create since that moment. Disabled by default.}
#'   \item{\code{disableautoopening()}}{This method disables auto opening of every visualisation that you create since that moment. Disabled by default.}
#'   \item{\code{line(series, index = NA, color = NA, label = NA, size = NA, xaxis = NA, yaxis = NA, logScaleX = "false", logScaleY = "false")}}{This method creates a line visualization for vector/matrix with each row representing a line, given in \code{series}.}
#'   \item{\code{scatter(x, y, color = NA, label = NA, size = NA, alpha = NA, xaxis = NA, yaxis = NA)}}{This method creates a scatterplot for points with coordinates given in vectors \code{x, y}.}
#'   \item{\code{linestacked(series, color = NA, label = NA, size = NA)}}{This method creates a plot of multiple lines given in matrix \code{series}, with an ability to hide and show every one of them.}
#'   \item{\code{force(matrix, color = NA, label = NA, size = NA)}}{This method creates a force plot for matrix given in \code{matrix}.}
#'   \item{\code{graph(x, y, matrix, color = NA, label = NA, size = NA)}}{This method creates a graph of points with coordinates given in \code{x, y} vectors, with connection given in \code{matrix} connectivity matrix.}
#'   \item{\code{map(regions, weights, colormap)}}{This method creates a world (or USA) map, marking regions given as a vector of abbreviations (3-char for countries, 2-char for states) in \code{regions} with weights given in \code{weights} vector and with \code{colormap} color (string from colorbrewer).}
#'   \item{\code{graphbundled(x, y, matrix, color = NA, label = NA, size = NA)}}{This method creates a bundled graph of points with coordinates given in \code{x, y} vectors, with connection given in \code{matrix} connectivity matrix. Lines on this graph are stacked a bit more than in the \code{graph} function.}
#'   \item{\code{matrix(matrix, colormap)}}{This method creates a visualization of matrix given in \code{matrix} parameter, with its contents used as weights for the colormap given in \code{colormap} (string from colorbrewer).}
#'   \item{\code{adjacency(matrix, label = NA)}}{This method creates a visualization for adjacency matrix given in \code{matrix} parameter.}
#'   \item{\code{scatterline(x, y, t, color = NA, label = NA, size = NA)}}{This method creates a scatterplot for coordinates in vectors \code{x, y} and assignes a line plot to every point on that plot. Each line is given as a row in \code{t} matrix.}
#'   \item{\code{scatter3(x, y, z, color = NA, label = NA, size = NA, alpha = NA)}}{This method creates a 3D scatterplot for coordinates given in vectors \code{x, y, z}.}
#'   \item{\code{image(imgpath)}}{This method uploads image from file \code{imgpath} to the server and creates a visualisation of it.}
#'   \item{\code{gallery(imgpathvector)}}{This method uploads images from vector of file paths \code{imgpathvector} to the server and creates a gallery of these images.}}

unit <- R6::R6Class("unit",
  # public ---------------------------------------------------------------------
  public = list(
    cycle = function(g_e_raw, g_i){
      ## updating g_e input
      private$g_e <- private$g_e + private$cyc_dt * private$g_e_dt *
        (g_e_raw - private$g_e)

      ## Finding membrane potential
      # excitatory, inhibitory and leak current
      i_e <- private$g_e * (private$e_rev_e - private$v_m)
      i_i <- g_i * (private$e_rev_i - private$v_m)
      i_l <- private$gc_l * (private$e_rev_l - private$v_m)
      i_net <- i_e + i_i + i_l

      # almost half-step method for updating v_m (i_adapt doesn't half step)
      v_m_h <- private$v_m + 0.5 * private$cyc_dt * private$v_m_dt *
        (i_net - private$i_adapt)
      i_e_h <- private$g_e * (private$e_rev_e - v_m_h)
      i_i_h <- g_i * (private$e_rev_i - v_m_h)
      i_l_h <- private$gc_l * (private$e_rev_l - v_m_h)
      i_net_h <- i_e_h + i_i_h + i_l_h

      private$v_m <- private$v_m + private$cyc_dt * private$v_m_dt *
        (i_net_h - private$i_adapt)
      private$v_m_eq <- private$v_m_eq + private$cyc_dt * private$v_m_dt *
        (i_net_h - private$i_adapt)

      ## Finding activation
      # finding threshold excitatory conductance
      g_e_thr <- (g_i * (private$e_rev_i - private$v_m_thr) +
                    private$gc_l * (private$e_rev_l - private$v_m_thr) -
                    private$i_adapt) / (private$v_m_thr - private$e_rev_e)

      # finding whether there's an action potential
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

      ## Updating adaptation current
      private$i_adapt <- private$i_adapt + private$cyc_dt *
        (private$i_adapt_dt * (private$v_m_gain * (private$v_m - private$e_rev_l)
                            - private$i_adapt) + private$spike * private$spike_gain)

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
      private$i_adapt <- 0
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
    g_e = 0,
    avg_ss = 0.2,
    v_m = 0.3,
    v_m_eq = 0.3,
    i_adapt = 0,
    spike = 0,
    # constant values
    g_e_dt = 1 / 1.4,   # time step constant for update of "g_e"
    cyc_dt = 1,       # time step constant for integration of cycle dynamics
    v_m_dt = 1 / 3.3,   # time step constant for membrane potential
    l_dn_dt = 1 / 2.5,  # time step constant for avg_l decrease
    i_adapt_dt = 1 / 144, # time step constant for adaptation
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
