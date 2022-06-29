# bisonR - An R package for Bayesian Inference of Social Networks <img src="man/figures/logo.png" width=128 align="right" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/JHart96/bisonR/workflows/R-CMD-check/badge.svg)](https://github.com/JHart96/bisonR/actions)
[![pkgdown](https://github.com/JHart96/bisonR/workflows/pkgdown/badge.svg)](https://github.com/JHart96/bisonR/actions)
[![Codecov test coverage](https://codecov.io/gh/JHart96/bisonR/branch/main/graph/badge.svg)](https://app.codecov.io/gh/JHart96/bisonR?branch=main)
<!-- badges: end -->

:warning: **This package is still in early stages of development and has not yet been fully tested**

bisonR is an R package implementing the BISoN framework for conducting Bayesian analysis of social networks. BISoN estimates uncertainty over edge weights in social networks from empirical data (such as observations) and builds networks with uncertainty. The networks can then be visualised *with* uncertainty, and can also be analysed with a fully Bayesian methodology using standard tools such as regression.

## Installation

### Installing CmdStanR

The bisonR package is written in R, but uses the Stan programming language to fit Bayesian models. Stan is a separate program, and interfaces with bisonR using an R package called cmdstanR. Stan and cmdstanR are installed in a different way to standard R packages, so require a specific series of installation steps. 

On Windows, before proceeding you may need to install the version of Rtools appropriate for your version of R. Rtools can be found here: https://cran.r-project.org/bin/windows/Rtools/.

The full instructions to install cmdstanR can be found at https://mc-stan.org/cmdstanr/. We've found that the following steps often work, but depending on your operating system and version of R, the process may be more involved:

```r
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
install_cmdstan()
```

### Installing bisonR

bisonR isn't currently on CRAN, but it can be installed from GitHub. To do this, make sure you have the `remotes` package installed. Then run the following command:

```r
remotes::install_github("JHart96/bisonR")
```

#### Development version

If you want to install the latest development version, bugs and all, you can install from the development branch. This option is not recommended for standard users.

```r
remotes::install_github("JHart96/bisonR@dev")
```
