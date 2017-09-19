<!-- README.md is generated from README.Rmd. Please edit that file -->
leabRa
======

This package provides the Leabra artificial neural network algorithm (O’Reilly, 1996) for R. Leabra stands for “local error driven and associative biologically realistic algorithm.” It is the Rolls Royce of artificial neural networks because it combines error driven learning and self organized learning in an elegant way, while focusing on a biologically plausible learning rule. If you have never heard of Leabra, you should read about it first. A good place to start is the computational cognitive neuroscience book (Part I), available at <https://grey.colorado.edu/CompCogNeuro/index.php/CCNBook/Main> (O’Reilly et al., 2016).

Installation
------------

<!--via CRAN: 


```r
#install.pacakges("leabRa")
```
-->
via git (and with building the vignette; this takes a while)

``` r
devtools::install_github("johannes-titz/leabRa", build_vignettes = TRUE)
```

Now you can load the package:

``` r
library(leabRa)
#> Welcome to leabRa
```

To see how it works, have a look at the vignette:

``` r
vignette("leabRa")
#> starting httpd help server ... done
```

And the help files:

``` r
?network
?layer
?unit
```

Citation
--------

Johannes Titz (2017). *leabRa: An R implementation of the artificial neural networks algorithm Leabra*. R package version 0.1.0. <https://github.com/johannes-titz/leabRa>
