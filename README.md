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

The full instructions to install cmdstanR can be found at https://mc-stan.org/cmdstanr/. We've found that the following steps often work, but depending on your operating system and version of R, the process may be more involved.

```r
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
install_cmdstan() # On networked PCs, specify a local directory here with the argument dir=path_to_local_directory
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

## Quick Start

You can check your `bisonR` installation and get started with the package using the code below:

```r
library(bisonR)

sim_data <- simulate_edge_model("binary", aggregated = TRUE)
df <- sim_data$df_sim

priors <- get_default_priors("binary")

fit_edge <- edge_model(
  (event | duration) ~ dyad(node_1_id, node_2_id), 
  data=df, 
  data_type="binary",
  priors=priors
)

summary(fit_edge)
```

This will fit a basic edge model model to simulated data and generate an output like this below:

```
=== Fitted BISoN edge model ===
Data type: binary
Formula: (event | duration) ~ dyad(node_1_id, node_2_id)
Number of nodes: 10
Number of dyads: 45
Directed: FALSE
=== Edge list summary ===
         median    5%   95%
1 <-> 2   0.156 0.014 0.530
1 <-> 3   0.702 0.219 0.973
2 <-> 3   0.793 0.356 0.982
...
```

For a more detailed example check out the *Getting Started* page here: https://jhart96.github.io/bisonR/articles/getting_started.html.
