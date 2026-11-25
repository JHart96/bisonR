# bisonR <img src="man/figures/hex.svg" width="140" align="right" alt="bisonR hex logo" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/JHart96/bisonR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/JHart96/bisonR/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/JHart96/bisonR/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/JHart96/bisonR/actions/workflows/pkgdown.yaml)
[![Codecov test coverage](https://codecov.io/gh/JHart96/bisonR/branch/main/graph/badge.svg)](https://app.codecov.io/gh/JHart96/bisonR?branch=main)
<!-- badges: end -->

**Bayesian Inference of Social Networks in R.**

bisonR implements the BISoN framework for analysing animal social networks. Instead of treating edge weights as fixed numbers, BISoN estimates each edge weight *with uncertainty* from the observation data (how often a pair was seen interacting, and how much they were observed). That uncertainty is then carried through everything downstream: network plots, network metrics, regressions, and comparisons between networks.

Documentation, tutorials and the function reference are at **<https://bisonr.org>**.

## Quick start

```r
library(bisonR)

# Simulated observations: one row per dyad, with the number of sampling
# periods each pair was seen together (event) out of those observed (duration).
sim_data <- simulate_bison_model("binary", aggregated = TRUE)
df <- sim_data$df_sim

fit_edge <- bison_model(
  (event | duration) ~ dyad(node_1_id, node_2_id),
  data = df,
  model_type = "binary",
  priors = get_default_priors("binary")
)

summary(fit_edge)      # Edge weights with credible intervals
plot_network(fit_edge) # Sociogram showing uncertainty in each edge
```

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
...
```

The [Getting started](https://bisonr.org/articles/getting_started.html) article walks through a complete analysis.

## Installation

bisonR uses [Stan](https://mc-stan.org/) through the [cmdstanr](https://mc-stan.org/cmdstanr/) package to fit models, so you need cmdstanr and CmdStan before installing bisonR.

### 1. Install cmdstanr and CmdStan

On Windows you may first need [Rtools](https://cran.r-project.org/bin/windows/Rtools/) for your version of R (see the [Windows installation](https://bisonr.org/articles/windows_installation.html) article).

```r
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
cmdstanr::check_cmdstan_toolchain(fix = TRUE)
cmdstanr::install_cmdstan() # On networked PCs, set dir = "a/local/folder"
```

The full instructions are in the [cmdstanr documentation](https://mc-stan.org/cmdstanr/articles/cmdstanr.html).

### 2. Install bisonR

bisonR isn't on CRAN. Install the latest release from GitHub with:

```r
# install.packages("remotes")
remotes::install_github("JHart96/bisonR")
```

The first time you fit a model, bisonR compiles its Stan models, which takes a minute or two. Compiled models are cached, so this only happens again when bisonR or CmdStan is updated.

### Development versions

The `main` branch always holds the latest release. Work on the next release happens on a `release/<version>` branch (for example `release/0.5.0`) and is merged into `main` when it's ready. To try an upcoming release before it's finished, bugs and all, install from its branch:

```r
remotes::install_github("JHart96/bisonR@release/0.5.0")
```

You can see which release branches exist on the [branches page](https://github.com/JHart96/bisonR/branches). Development versions aren't recommended for published analyses.

## Getting help and contributing

* Questions and ideas: [GitHub Discussions](https://github.com/JHart96/bisonR/discussions)
* Bugs and feature requests: [GitHub Issues](https://github.com/JHart96/bisonR/issues)
* Contributing code or documentation: see the [contributing guide](.github/CONTRIBUTING.md)

## Citing bisonR

If you use bisonR in published work, please cite the BISoN paper:

> Hart, J. D. A., Weiss, M. N., Franks, D. W., & Brent, L. J. N. (2023). BISoN: A Bayesian framework for inference of social networks. *Methods in Ecology and Evolution*, 14(9), 2411–2420. <https://doi.org/10.1111/2041-210X.14171>

Run `citation("bisonR")` in R for the citation in BibTeX format.
