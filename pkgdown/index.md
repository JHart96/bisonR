<div class="bison-hero">
<div class="bison-hero-text">
<h2 class="bison-title">bisonR</h2>
<p class="bison-subtitle">Bayesian inference of animal social networks in R</p>
<p class="bison-lead">bisonR implements the BISoN framework (Hart et al. 2023). Edge weights are estimated from observation data as posterior distributions, so the precision of each estimate reflects how much the dyad was sampled. These posteriors can then be used in network metrics, regression models and comparisons between networks.</p>
<div class="bison-actions">
<a class="btn btn-primary" href="articles/getting_started.html">Getting started</a>
<a class="btn btn-outline-primary" href="reference/index.html">Function reference</a>
</div>
</div>
<div class="bison-hero-figure">
<img src="hero-network.png" alt="A network plotted with plot_network(). Each edge has a dark line for the lower bound of its credible interval and a lighter band for the upper bound." />
</div>
</div>

<div class="bison-features">
<div class="bison-feature">
<h3>Edge weight models</h3>
<p>Binary and count models of social events with observation effort, optional fixed and random effects, partial pooling and zero inflation. Conjugate models are available for simple cases.</p>
<a href="articles/preparing_data.html">Preparing data</a>
</div>
<div class="bison-feature">
<h3>Diagnostics and plots</h3>
<p>MCMC trace plots, prior and posterior predictive checks, and network plots that display the credible interval of each edge.</p>
<a href="reference/plot_predictions.html">Predictive checks</a>
</div>
<div class="bison-feature">
<h3>Downstream analyses</h3>
<p>Posterior distributions of node and network metrics, dyadic and nodal regression with brms, edge mixture models and non-random association tests.</p>
<a href="articles/bison_brms.html">Regression with brms</a>
</div>
</div>

## A first model

```r
library(bisonR)

sim_data <- simulate_bison_model("binary", aggregated = TRUE)

fit_edge <- bison_model(
  (event | duration) ~ dyad(node_1_id, node_2_id),
  data = sim_data$df_sim,
  model_type = "binary",
  priors = get_default_priors("binary")
)

summary(fit_edge)                                # Edge weights with credible intervals
plot_network(fit_edge, lwd = 5)                  # Sociogram showing uncertainty
extract_metric(fit_edge, "node_strength")        # Posterior of each node's strength
```

The left-hand side of the formula, `(event | duration)`, says how many social events were seen and how much each dyad was observed. The right-hand side, `dyad(node_1_id, node_2_id)`, asks for an edge weight for every pair. The [Getting started](articles/getting_started.html) article explains each step.

## Installation

bisonR fits models with [Stan](https://mc-stan.org/) via the [cmdstanr](https://mc-stan.org/cmdstanr/) package, so install cmdstanr and CmdStan first:

```r
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
cmdstanr::check_cmdstan_toolchain(fix = TRUE)
cmdstanr::install_cmdstan()
```

Then install bisonR from GitHub:

```r
# install.packages("remotes")
remotes::install_github("JHart96/bisonR")
```

On Windows you may need [Rtools](https://cran.r-project.org/bin/windows/Rtools/) first; the [Windows installation](articles/windows_installation.html) article covers common problems. Upcoming releases are developed on `release/<version>` branches and can be installed with, for example, `remotes::install_github("JHart96/bisonR@release/1.1.0")`.

## Citing bisonR

> Hart, J. D. A., Weiss, M. N., Franks, D. W., & Brent, L. J. N. (2023). BISoN: A Bayesian framework for inference of social networks. *Methods in Ecology and Evolution*, 14(9), 2411–2420. <https://doi.org/10.1111/2041-210X.14171>

## Getting help

Ask questions in [Discussions](https://github.com/JHart96/bisonR/discussions) and report bugs in [Issues](https://github.com/JHart96/bisonR/issues). Contributions are very welcome; see the [contributing guide](https://github.com/JHart96/bisonR/blob/main/.github/CONTRIBUTING.md).
