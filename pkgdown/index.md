<div class="bison-hero">
<div class="bison-hero-text">
<p class="bison-eyebrow">bisonR · Bayesian Inference of Social Networks</p>
<h2 class="bison-title">Social networks, with honest uncertainty</h2>
<p class="bison-lead">bisonR estimates animal social networks from observation data, so a pair seen together twice isn't treated as confidently as a pair seen together two hundred times. That uncertainty flows through to everything you do next: plots, network metrics, regressions and comparisons between networks.</p>
<div class="bison-actions">
<a class="btn btn-primary btn-lg" href="articles/getting_started.html">Get started</a>
<a class="btn btn-outline-primary btn-lg" href="#installation">Install</a>
</div>
</div>
<div class="bison-hero-figure">
<img src="hero-network.png" alt="A social network drawn by bisonR. Each edge has a dark line showing the lower bound of its weight and a lighter band showing the upper bound." />
<p class="bison-figure-caption">Dark lines show how weak each edge could be; light bands show how strong.</p>
</div>
</div>

<div class="bison-features">
<div class="bison-feature">
<h3>Estimate edges from raw data</h3>
<p>Fit binary or count edge models to observations, with sampling effort, fixed and random effects, partial pooling and zero inflation. Fast conjugate models give exact answers in seconds.</p>
<a href="articles/preparing_data.html">Preparing your data →</a>
</div>
<div class="bison-feature">
<h3>Check and visualise</h3>
<p>Trace plots, posterior predictive checks and network plots that show the uncertainty in every edge, so you can see what the data do and don't tell you.</p>
<a href="reference/plot_network.html">Plotting networks →</a>
</div>
<div class="bison-feature">
<h3>Analyse with uncertainty</h3>
<p>Calculate posterior distributions of network metrics, and fit dyadic, nodal and network-level regressions with brms, without throwing away the uncertainty in the network.</p>
<a href="articles/bison_brms.html">Regression with brms →</a>
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

On Windows you may need [Rtools](https://cran.r-project.org/bin/windows/Rtools/) first; the [Windows installation](articles/windows_installation.html) article covers common problems. Upcoming releases are developed on `release/<version>` branches and can be installed with, for example, `remotes::install_github("JHart96/bisonR@release/0.5.0")`.

## Citing bisonR

> Hart, J. D. A., Weiss, M. N., Franks, D. W., & Brent, L. J. N. (2023). BISoN: A Bayesian framework for inference of social networks. *Methods in Ecology and Evolution*, 14(9), 2411–2420. <https://doi.org/10.1111/2041-210X.14171>

## Getting help

Ask questions in [Discussions](https://github.com/JHart96/bisonR/discussions) and report bugs in [Issues](https://github.com/JHart96/bisonR/issues). Contributions are very welcome; see the [contributing guide](https://github.com/JHart96/bisonR/blob/main/.github/CONTRIBUTING.md).
