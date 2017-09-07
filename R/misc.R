.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to LeabRa")
}

.onLoad <- function(libname, pkgname) {
  op <- options()
  op.devtools <- list(
    devtools.path = "~/R-dev",
    devtools.install.args = "",
    devtools.name = "Johannes Titz",
    devtools.desc.author =
      '"Johannes Titz <johannes.titz@gmail.com> [aut, cre]"',
    #devtools.desc.author = '"Sergio Verduzco-Flores <> [aut]"',
    devtools.desc.license = "GPL-3",
    devtools.desc.suggests = NULL,
    devtools.desc = list()
  )
  toset <- names(op.devtools) %in% names(op)
  if (any(toset)) options(op.devtools[toset])

  invisible()
}

isempty <- function(x){
  length(x) == 0
}
