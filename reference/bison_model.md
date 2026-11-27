# Fit an edge weight model

Fits a BISoN edge weight model, which estimates a posterior distribution
for the edge weight of each dyad from observations of social events and
the observation effort for each dyad. An edge weight is a dyad's
underlying tendency to interact or associate, as a probability or a
rate. The posterior distribution describes the values consistent with
the data and priors, so edges observed rarely get wide uncertainty. This
is usually the first step of an analysis: the fitted model is then used
to plot the network, calculate network metrics and fit regressions.

## Usage

``` r
bison_model(
  formula,
  data,
  model_type = c("binary", "count", "binary_conjugate", "count_conjugate"),
  directed = FALSE,
  partial_pooling = FALSE,
  zero_inflated = FALSE,
  priors = NULL,
  refresh = 0,
  mc_cores = 4,
  iter_sampling = 1000,
  iter_warmup = 1000,
  priors_only = FALSE,
  chains = 4,
  ...
)
```

## Arguments

- formula:

  A formula of the form `(event | effort) ~ dyad(node_1, node_2)`. The
  left-hand side gives the column of social events (0/1 or counts) and
  the observation effort (a column, or a single number). Effort is how
  much each row was observed, such as the number of sampling periods or
  the hours of observation. The right-hand side gives the terms that
  explain the events. See *Formula syntax* below.

- data:

  A data frame with one row per observation (disaggregated) or one row
  per dyad (aggregated).

- model_type:

  The type of edge model: `"binary"` (events are present/absent in each
  sampling period, edge weights are probabilities), `"count"` (events
  are counts, edge weights are rates per unit effort), or the fast
  conjugate versions `"binary_conjugate"` and `"count_conjugate"`, which
  only support a `dyad()` term.

- directed:

  `TRUE` or `FALSE`: whether the network is directed. For undirected
  networks, rows for A-B and B-A refer to the same dyad.

- partial_pooling:

  Whether to estimate the spread of edge weights from the data, so that
  information is shared between edges. How far estimates for dyads with
  few observations are pulled towards the mean of the edge prior then
  depends on how much edge weights vary across the network.

- zero_inflated:

  Whether to use a zero-inflated model, for data with more zeros than
  the edge weights alone would predict (see the `"stats"` check in
  [`plot_predictions()`](https://bisonr.org/reference/plot_predictions.md)).
  Zero inflation applies to each row of the data: each row has
  probability `zero_prob` of being a structural zero (a zero regardless
  of the edge weight). For binary data with one sampling period per row
  this can't be separated from the edge weights, so aggregate the data
  first.

- priors:

  A list of priors, in the format returned by
  [`get_default_priors()`](https://bisonr.org/reference/get_default_priors.md).
  A prior describes which values of a parameter are plausible before
  seeing the data. If `NULL`, the defaults are used with a message;
  check them with
  [`prior_check()`](https://bisonr.org/reference/prior_check.md) and
  [`prior_predictive_check()`](https://bisonr.org/reference/prior_predictive_check.md).

- refresh:

  How often the sampler prints progress (0 for silent).

- mc_cores:

  Number of chains to run in parallel.

- iter_sampling:

  Number of posterior draws to keep per chain, after warmup. For
  conjugate models, `iter_sampling * chains` draws are taken.

- iter_warmup:

  Number of warmup iterations per chain, during which the sampler tunes
  itself. These draws are discarded.

- priors_only:

  If `TRUE`, ignore the data and sample from the priors (used for prior
  predictive checks).

- chains:

  Number of MCMC chains.

- ...:

  Further arguments passed to the `$sample()` method of
  [cmdstanr::CmdStanModel](https://mc-stan.org/cmdstanr/reference/CmdStanModel.html),
  such as `seed` or `adapt_delta`.

## Value

A `bison_model` object. Use
[`summary.bison_model()`](https://bisonr.org/reference/summary.bison_model.md)
or [`get_edgelist()`](https://bisonr.org/reference/get_edgelist.md) to
see the edge weights,
[`plot_network()`](https://bisonr.org/reference/plot_network.md) to
visualise them,
[`plot_trace()`](https://bisonr.org/reference/plot_trace.md) and
[`plot_predictions()`](https://bisonr.org/reference/plot_predictions.md)
for diagnostics, and
[`extract_metric()`](https://bisonr.org/reference/extract_metric.md) or
[`bison_brm()`](https://bisonr.org/reference/bison_brm.md) for
downstream analyses.

## Details

Models other than the conjugate ones are fitted by MCMC (Markov chain
Monte Carlo) in Stan, which draws samples from the posterior in several
independent runs called chains. Conjugate models have a posterior that
can be calculated exactly, so they are fitted in seconds without Stan,
but they can't include other effects.

## Formula syntax

- `dyad(node_1, node_2)`: estimate an edge weight for each dyad. The two
  columns hold node names or IDs. If they are both factors with the same
  levels, nodes that never appear in the data are kept.

- `x`: a fixed effect of column `x`. Factors get one parameter per
  level.

- `(1 | g)`: a random intercept for each level of column `g`. Random
  effects are centred on 0, so a model without a `dyad()` term or factor
  fixed effect needs an intercept (`1`) to set the overall rate of
  events.

- `1` or `0`: include or exclude a global intercept (excluded by
  default).

Including other terms changes the interpretation of the edge weights,
which then represent each dyad's tendency to interact after controlling
for those effects. On the natural scale (e.g. in
[`get_edgelist()`](https://bisonr.org/reference/get_edgelist.md) and
[`extract_metric()`](https://bisonr.org/reference/extract_metric.md)),
they are the probability or rate for an observation in which the other
effects (and any intercept) are 0, so continuous predictors may be
easier to interpret if they are centred.

## See also

[`get_default_priors()`](https://bisonr.org/reference/get_default_priors.md),
[`prior_check()`](https://bisonr.org/reference/prior_check.md),
[`summary.bison_model()`](https://bisonr.org/reference/summary.bison_model.md),
[`get_effects()`](https://bisonr.org/reference/get_effects.md),
[`model_comparison()`](https://bisonr.org/reference/model_comparison.md),
[`bison_save()`](https://bisonr.org/reference/bison_save.md)

## Examples

``` r
sim_data <- simulate_bison_model("binary", aggregated = TRUE)
df <- sim_data$df_sim

# Fast conjugate model, no MCMC needed
fit_edge <- bison_model(
  (event | effort) ~ dyad(node_1_id, node_2_id),
  data = df,
  model_type = "binary_conjugate",
  priors = get_default_priors("binary_conjugate")
)
summary(fit_edge)
#> === Fitted BISoN edge model ===
#> Data type: binary
#> Formula: (event | effort) ~ dyad(node_1_id, node_2_id)
#> Number of nodes: 10
#> Number of dyads: 45
#> Directed: FALSE
#> === Edge list summary ===
#>          median    5%   95%
#> 1 <-> 2   0.385 0.091 0.751
#> 1 <-> 3   0.159 0.011 0.530
#> 1 <-> 4   0.212 0.017 0.626
#> 1 <-> 5   0.306 0.075 0.643
#> 1 <-> 6   0.614 0.254 0.902
#> 1 <-> 7   0.498 0.189 0.823
#> 1 <-> 8   0.291 0.031 0.763
#> 1 <-> 9   0.162 0.014 0.520
#> 1 <-> 10  0.503 0.195 0.812
#> 2 <-> 3   0.584 0.283 0.850
#> 2 <-> 4   0.739 0.429 0.937
#> 2 <-> 5   0.610 0.239 0.904
#> 2 <-> 6   0.708 0.229 0.972
#> 2 <-> 7   0.384 0.100 0.757
#> 2 <-> 8   0.792 0.366 0.984
#> 2 <-> 9   0.735 0.421 0.937
#> 2 <-> 10  0.679 0.343 0.925
#> 3 <-> 4   0.611 0.253 0.902
#> 3 <-> 5   0.267 0.063 0.572
#> 3 <-> 6   0.692 0.349 0.922
#> 3 <-> 7   0.690 0.357 0.921
#> 3 <-> 8   0.793 0.369 0.984
#> 3 <-> 9   0.682 0.348 0.917
#> 3 <-> 10  0.871 0.539 0.990
#> 4 <-> 5   0.683 0.346 0.920
#> 4 <-> 6   0.842 0.483 0.988
#> 4 <-> 7   0.502 0.132 0.863
#> 4 <-> 8   0.709 0.223 0.974
#> 4 <-> 9   0.202 0.015 0.617
#> 4 <-> 10  0.293 0.024 0.772
#> 5 <-> 6   0.612 0.249 0.902
#> 5 <-> 7   0.380 0.099 0.750
#> 5 <-> 8   0.385 0.100 0.750
#> 5 <-> 9   0.287 0.025 0.787
#> 5 <-> 10  0.895 0.618 0.991
#> 6 <-> 7   0.262 0.062 0.574
#> 6 <-> 8   0.497 0.135 0.877
#> 6 <-> 9   0.312 0.082 0.661
#> 6 <-> 10  0.310 0.075 0.652
#> 7 <-> 8   0.106 0.008 0.394
#> 7 <-> 9   0.704 0.205 0.975
#> 7 <-> 10  0.294 0.028 0.777
#> 8 <-> 9   0.384 0.102 0.750
#> 8 <-> 10  0.870 0.548 0.991
#> 9 <-> 10  0.501 0.140 0.857

# Full model fitted with Stan, which allows extra effects
fit_edge <- bison_model(
  (event | effort) ~ dyad(node_1_id, node_2_id),
  data = df,
  model_type = "binary",
  priors = get_default_priors("binary")
)
#> Running MCMC with 4 parallel chains...
#> 
#> Chain 1 finished in 0.2 seconds.
#> Chain 2 finished in 0.3 seconds.
#> Chain 3 finished in 0.2 seconds.
#> Chain 4 finished in 0.2 seconds.
#> 
#> All 4 chains finished successfully.
#> Mean chain execution time: 0.2 seconds.
#> Total execution time: 0.4 seconds.
#> 
```
