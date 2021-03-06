% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/unit.R
\docType{class}
\name{unit}
\alias{unit}
\title{Leabra unit (neuron) class}
\format{\code{\link{R6Class}} object.}
\usage{
unit
}
\value{
Object of \code{\link{R6Class}} with methods for calculating neuron
  activation changes.
}
\description{
This class simulates a biologically realistic neuron (also called unit) in
the Leabra framework. When you use the layer class, you will see that a
\link{layer} object has a variable (field) \code{units}, which is a list of
unit objects.
}
\section{Fields}{

\describe{
\item{\code{activation}}{Percentage activation ("firing rate") of the unit, which is
sent to other units, think of it as a percentage of how many neurons are
active in a microcolumn of 100 neurons.}

\item{\code{avg_s}}{Short-term running average activation, integrates over avg_ss (a
private variable, which integrates over activation), represents plus phase
learning signal.}

\item{\code{avg_m}}{Medium-term running average activation, integrates over avg_s,
represents minus phase learning signal.}

\item{\code{avg_l}}{Long-term running average activation, integrates over avg_m,
drives long-term floating average for self-organized learning.}

\item{\code{unit_number}}{Number of unit in layer, if the unit is not created within
a layer, this value will be 1.}
}}

\section{Methods}{

\describe{
  \item{\code{new()}}{Creates an object of this class with default
  parameters.}

  \item{\code{cycle(g_e_raw, g_i)}}{Cycles 1 ms with given excitatory
  conductance \code{g_e_raw} and inhibitory conductance \code{g_i}.
  Excitatory conductance depends on the connection weights to other units and
  the activity of those other units. Inhibitory conductance depends on
  feedforward and feedback inhibition. See \link{layer} cycle method.

    \describe{
      \item{\code{g_e_raw}}{Raw excitatory conductance. The actual excitatory
      conductance will incrementally approach this value with every cycle.}

      \item{\code{g_i}}{Inhibitory conductance.}
      }
  }

  \item{\code{clamp_cycle(activation)}}{Clamps the value of \code{activation}
  to the \code{activation} variable of the unit without any time integration.
  Then updates averages (\code{avg_ss}, \code{avg_s}, \code{avg_m}). This is
  usually done when presenting external input.

    \describe{
      \item{\code{activation}}{Activation to clamp.}
      }
  }

  \item{\code{updt_avg_l()}}{Updates the variable \code{avg_l}. This usually
  happens before the weights are changed in the network (after the plus
  phase), and not every cycle.}

  \item{\code{get_vars(show_dynamics = TRUE, show_constants =
  FALSE)}}{Returns a data frame with 1 row with the current state of all the
  variables of the unit. You can choose whether you want dynamic values and /
  or constant values. This might be useful if you want to analyze what
  happens in a unit, which would otherwise not be possible, because most of
  the variables (fields) are private in this class.

    \describe{
      \item{\code{show_dynamics}}{Should dynamic values be shown? Default is
      TRUE
      }

      \item{\code{show_constants}}{Should constant values be shown? Default
      is FALSE
      }
    }
  }
}
}

\examples{
u <- unit$new() # creates a new unit with default leabra values

print(u) # a lot of private values
u$v # private values cannot be accessed
# if you want to see alle variables, you need to use the function
u$get_vars(show_dynamics = TRUE, show_constants = TRUE)

# let us clamp the activation to 0.7
u$activation
u$clamp_cycle(0.7)
c(u$activation, u$avg_s, u$avg_m, u$avg_l)
# activation is indeed 0.7, but avg_l was not updated, this only happens
# before the weights are changed, let us update it now
u$updt_avg_l()
c(u$activation, u$avg_s, u$avg_m, u$avg_l)
# seems to work

# let us run 10 cycles with unclamped activation and output the activation
# produced because of changes in conductance
u <- unit$new()
cycle_number <- 1:10
result <- lapply(cycle_number, function(x)
                 u$cycle(g_e_raw = 0.5, g_i = 0.5)$get_vars())
# make a data frame out of the list
result <- plyr::ldply(result)
# plot activation
plot(result$activation, type = "b", xlab = "cycle", ylab = "activation")
# add conductance g_e to plot, should approach g_e_raw
lines(result$g_e, type = "b", col = "blue")

}
\references{
O'Reilly, R. C., Munakata, Y., Frank, M. J., Hazy, T. E., and
  Contributors (2016). Computational Cognitive Neuroscience. Wiki Book, 3rd
  (partial) Edition. URL: \url{http://ccnbook.colorado.edu}

Have also a look at
  \url{https://grey.colorado.edu/emergent/index.php/Leabra} (especially the
  link to the 'MATLAB' code) and \url{https://en.wikipedia.org/wiki/Leabra}
}
\keyword{data}
