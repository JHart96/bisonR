## bisonR

Bayesian inference of animal social networks in R

bisonR implements the BISoN framework (Hart et al. 2023). It estimates
each edge weight from observation data as a posterior distribution (the
range of plausible values given the data and prior), so the precision of
each estimate reflects how much that dyad was sampled. These posteriors
can then be used to calculate network metrics, fit regression models and
compare networks.

[Getting started](https://bisonr.org/articles/getting_started.md)
[Function reference](https://bisonr.org/reference/index.md)

![A social network plotted with plot_network(), in which the width of
each edge shows the uncertainty in its weight.](hero-network.png)

### Edge weight models

Binary and count models of social events that account for observation
effort, with optional fixed and random effects, partial pooling and zero
inflation. Faster conjugate models, fitted without MCMC (the sampling
method that Stan uses), are available for simple cases.

[Preparing data](https://bisonr.org/articles/preparing_data.md)

### Diagnostics and plots

Trace plots of the MCMC sampler, prior and posterior predictive checks,
and network plots that show the credible interval of each edge.

[Predictive checks](https://bisonr.org/reference/plot_predictions.md)

### Downstream analyses

Posterior distributions of node and network metrics, dyadic and nodal
regression with brms, edge mixture models and non-random association
tests.

[Regression with brms](https://bisonr.org/articles/bison_brms.md)

## A first model

``` r

library(bisonR)

sim_data <- simulate_bison_model("binary", aggregated = TRUE)

fit_edge <- bison_model(
  (event | effort) ~ dyad(node_1_id, node_2_id),
  data = sim_data$df_sim,
  model_type = "binary",
  priors = get_default_priors("binary")
)

summary(fit_edge)                                # Edge weights with credible intervals
plot_network(fit_edge, lwd = 5)                  # Sociogram showing uncertainty
extract_metric(fit_edge, "node_strength")        # Posterior of each node's strength
```

The left-hand side of the formula, `(event | effort)`, gives the number
of social events seen and how much each dyad was observed. The
right-hand side, `dyad(node_1_id, node_2_id)`, estimates an edge weight
for each dyad (pair of individuals). The [Getting
started](https://bisonr.org/articles/getting_started.md) article
explains each step.

## Installation

bisonR fits models with [Stan](https://mc-stan.org/) through the
[cmdstanr](https://mc-stan.org/cmdstanr/) package, so install cmdstanr
and CmdStan first:

``` r

install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
cmdstanr::check_cmdstan_toolchain(fix = TRUE)
cmdstanr::install_cmdstan()
```

Then install bisonR from GitHub:

``` r

# install.packages("remotes")
remotes::install_github("JHart96/bisonR")
```

On Windows you may need to install
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) first. The
[Windows
installation](https://bisonr.org/articles/windows_installation.md)
article covers common problems. Upcoming releases are developed on
`release/<version>` branches. To install one, use, for example,
`remotes::install_github("JHart96/bisonR@release/1.1.0")`.

## Citing bisonR

> Hart, J. D. A., Weiss, M. N., Franks, D. W., & Brent, L. J. N. (2023).
> BISoN: A Bayesian framework for inference of social networks. *Methods
> in Ecology and Evolution*, 14(9), 2411–2420.
> <https://doi.org/10.1111/2041-210X.14171>

## Getting help

Ask questions in
[Discussions](https://github.com/JHart96/bisonR/discussions) and report
bugs in [Issues](https://github.com/JHart96/bisonR/issues).
Contributions are welcome; see the [contributing
guide](https://github.com/JHart96/bisonR/blob/main/.github/CONTRIBUTING.md).
