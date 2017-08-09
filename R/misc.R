.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to LeabRa")
}

.onLoad <- function(libname, pkgname) {
  op <- options()
  op.devtools <- list(
    devtools.path = "~/R-dev",
    devtools.install.args = "",
    devtools.name = "Johannes Titz",
    devtools.desc.author = '"Johannes Titz <johannes.titz@gmail.com> [aut, cre]"',
    #devtools.desc.author = '"Sergio Verduzco-Flores <> [aut]"',
    devtools.desc.license = "GPL-3",
    devtools.desc.suggests = NULL,
    devtools.desc = list()
    #devtools::create()
  )
  toset <- names(op.devtools) %in% names(op)
  if(any(toset)) options(op.devtools[toset])

  invisible()
}

# packages needed
devtools::use_package("signal")
devtools::use_package("R.cache")
devtools::use_package("R6")
devtools::use_package("dplyr")
devtools::use_package("plyr")

isempty <- function(x)
{
  length(x) == 0
}

mult_list_vec <- function(x, y){
  mapply("*", x, y)
}

#' m_mapply
#'
#' mapply version for matrix of matrices, returns a matrix
m_mapply <- function(fun, x, y){
  matrix(mapply(fun, x, y), ncol = ncol(x))
}

# f = nxx1(x) calculates the noisy x/(x+1) function for all values in the
# vector "x". The returned values come from the convolution of x/(x+1)
# with a Gaussian function.
# To avoid calculating the convolution every time, the first time the function
# is called a vector "nxoxp1" is created, corresponding to the values of nxx1
# at all the x in the vector "nxx1_dom". Once these vectors are in the
# workspace, subsequent calls to nxx1 use interpolation with these vectors in
# order to calculate their return values.

create_nxx1 <- function(){
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
}
