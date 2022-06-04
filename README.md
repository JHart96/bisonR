# bisonR - An R package for Bayesian Inference of Social Networks <img src="man/figures/logo.png" width=128 align="right" />

<!-- badges: start -->
  [![R-CMD-check](https://github.com/JHart96/bisonR/workflows/R-CMD-check/badge.svg)](https://github.com/JHart96/bisonR/actions)
  <!-- badges: end -->

:warning: **Do not use this package! (yet): This package is still in very early stages of development and not ready for real use. If you're looking to use BISoN for social network analysis check out the examples repository in Stan and INLA here: [github.com/JHart96/bison_examples](https://github.com/JHart96/bison_examples). **

bisonR is an R package implementing the BISoN framework for conducting Bayesian analysis of social networks. BISoN estimates uncertainty over edge weights in social networks from empirical data (such as observations) and builds networks with uncertainty. The networks can then be visualised *with* uncertainty, and can also be analysed with a fully Bayesian methodology using standard tools such as regression.

## Installation

bisonR isn't currently on CRAN, but it can be installed from GitHub. To do this, make sure you have the `devtools` package installed. Then run the following command:

```
devtools::install_github("JHart96/bisonR")
```
