# bisonR - An R package for Bayesian Inference of Social Networks <img src="man/figures/logo.png" width=128 align="right" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/JHart96/bisonR/workflows/R-CMD-check/badge.svg)](https://github.com/JHart96/bisonR/actions)
[![pkgdown](https://github.com/JHart96/bisonR/workflows/pkgdown/badge.svg)](https://github.com/JHart96/bisonR/actions)
[![Codecov test coverage](https://codecov.io/gh/JHart96/bisonR/branch/main/graph/badge.svg)](https://app.codecov.io/gh/JHart96/bisonR?branch=main)
<!-- badges: end -->

:warning: **This package is still in early stages of development and has not yet been fully tested**

bisonR is an R package implementing the BISoN framework for conducting Bayesian analysis of social networks. BISoN estimates uncertainty over edge weights in social networks from empirical data (such as observations) and builds networks with uncertainty. The networks can then be visualised *with* uncertainty, and can also be analysed with a fully Bayesian methodology using standard tools such as regression.

## Installation

cmdstanR needs to be installed for the bisonR package to work. To install cmdstanR, follow the instructions at https://mc-stan.org/cmdstanr/, then run `cmdstanr::install_cmdstan()`.

bisonR isn't currently on CRAN, but it can be installed from GitHub. To do this, make sure you have the `remotes` package installed. Then run the following command:

```
remotes::install_github("JHart96/bisonR")
```
