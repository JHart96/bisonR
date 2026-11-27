# Regression with brms

[`bison_brm()`](https://bisonr.org/reference/bison_brm.md) fits
regression models with the [brms](https://paulbuerkner.com/brms/)
package in which a network variable, such as an edge weight, a node
metric or a global metric, is the response or a predictor. The
uncertainty in the network is included in the results.

[`bison_brm()`](https://bisonr.org/reference/bison_brm.md) takes a
number of posterior draws from the edge model (each draw is one possible
network that is consistent with the data). It calculates the network
variable for each draw, fits the brms model to each resulting dataset
and combines the results. This is similar to multiple imputation, a
method for missing data in which a model is fitted to several plausible
versions of a dataset and the results are pooled. It takes the same
arguments as
[`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html), as
well as the edge model(s) and data, and returns a `brmsfit` object.

## Formula syntax

Network variables are specified with `bison()` in the formula. Inside
`bison()`, give the metric name (see
[`vignette("network_metrics")`](https://bisonr.org/articles/network_metrics.md))
followed by the column(s) in the data that identify the dyad, node or
network for each row:

| Standard formula | [`bison_brm()`](https://bisonr.org/reference/bison_brm.md) formula | Row of the data frame is a… |
|----|----|----|
| `edge_weight ~ x` | `bison(edge_weight(node_1_id, node_2_id)) ~ x` | dyad |
| `y ~ edge_weight` | `y ~ bison(edge_weight(node_1_id, node_2_id))` | dyad |
| `strength ~ x` | `bison(node_strength(node_id)) ~ x` | node |
| `y ~ strength` | `y ~ bison(node_strength(node_id))` | node |
| `degree ~ x` | `bison(node_degree[0.2](node_id)) ~ x` | node |
| `cv ~ x` | `bison(global_cv(bison_network)) ~ x` | network |

When a list of edge models is given,
[`bison_brm()`](https://bisonr.org/reference/bison_brm.md) adds a factor
`bison_network` to the data, giving the position of each edge model in
the list. It can be used like any other variable, for example to compare
networks or as a random effect.

Some notes:

- **Scale of edge weights.** In
  [`bison_brm()`](https://bisonr.org/reference/bison_brm.md), edge
  weights are on the link scale: log-odds (logit) for binary models and
  log rates for count models. Node and global metrics are calculated
  from edge weights on the natural scale.
- **Scaling.** `scale = "z"` converts the network variable to z-scores
  within each posterior draw (subtracting the mean and dividing by the
  standard deviation), and `scale = "centre"` only subtracts the mean.
  Standardising can make priors easier to set. For edge and node
  metrics, scaling is done separately within each network. This removes
  differences in the average level between networks, so use
  `scale = "none"` to compare levels between networks. For global
  metrics, scaling is done across networks.
- **Selecting rows.** `subset` selects rows of the data with a logical
  expression, e.g. `subset = sex == "F"`. Network metrics are still
  calculated from the whole network.
- **Number of draws.** `num_draws` sets the number of posterior draws
  used. A model is fitted to each, so more draws take longer. We suggest
  at least 100 for final analyses.
- **R-hat warnings.** R-hat is a convergence diagnostic that compares
  the MCMC chains (independent runs of the sampler). Because the results
  of separate models are combined, brms may warn about high R-hat
  values. This is expected, and convergence of the individual models can
  be checked with `fit$rhats`.

## Dyadic regression

``` r

library(bisonR)

sim_data <- simulate_bison_model("binary", aggregated = TRUE)
df_edge <- sim_data$df_sim
fit_edge <- bison_model(
  (event | effort) ~ dyad(node_1_id, node_2_id),
  data = df_edge,
  model_type = "binary",
  priors = get_default_priors("binary")
)
```

### Edge weight as the response

Edges that share a node are not independent. A multi-membership random
effect, `(1 | mm(node_1_id, node_2_id))`, accounts for this by
estimating an effect for each node that applies to every dyad the node
is part of. It should usually be included in dyadic regressions.
[`bison_brm()`](https://bisonr.org/reference/bison_brm.md) gives a
warning if it is missing.

``` r

fit_dyadic <- bison_brm(
  bison(edge_weight(node_1_id, node_2_id)) ~ age_diff + (1 | mm(node_1_id, node_2_id)),
  fit_edge,
  df_edge
)
summary(fit_dyadic)
```

### Edge weight as a predictor

``` r

bison_brm(
  outcome ~ bison(edge_weight(node_1_id, node_2_id)) + (1 | mm(node_1_id, node_2_id)),
  fit_edge,
  df_edge
)
```

### Comparing edge weights across networks

Give a list of edge models and a matching list of data frames.
`bison_network` identifies the network of each row:

``` r

bison_brm(
  bison(edge_weight(node_1_id, node_2_id)) ~ bison_network + (1 | mm(node_1_id, node_2_id)),
  list(fit_edge_control, fit_edge_treatment),
  list(df_edge_control, df_edge_treatment)
)
```

## Nodal regression

The data frame has one row per node, with a column containing the node
names used in the edge model:

``` r

df_nodal <- data.frame(
  node = levels(df_edge$node_1_id),
  trait = rnorm(10)
)
```

### Node metric as a predictor

``` r

bison_brm(
  trait ~ bison(node_eigen(node)),
  fit_edge,
  df_nodal
)
```

### Node metric as the response

``` r

bison_brm(
  bison(node_strength(node)) ~ trait,
  fit_edge,
  df_nodal
)
```

### Controlling for the network across several networks

``` r

bison_brm(
  bison(node_eigen(node)) ~ trait + (1 | bison_network),
  list(fit_edge_1, fit_edge_2),
  list(df_nodal_1, df_nodal_2)
)
```

## Global regression

Global metrics need a list of edge models, one per network. The data can
be a list of one-row data frames, or a single data frame with one row
per edge model in the same order:

``` r

df_global <- data.frame(
  bison_network = factor(1:2),
  condition = c("before", "after")
)
bison_brm(
  bison(global_cv(bison_network)) ~ condition,
  list(fit_edge_before, fit_edge_after),
  df_global
)
```

Each network contributes one observation to a global regression, so with
only a few networks the results should be interpreted with care.

## Priors

Priors describe which parameter values are plausible before the data are
seen.
[`bison_brm_get_prior()`](https://bisonr.org/reference/bison_brm_get_prior.md)
returns the default brms priors for a model as a `brmsprior` object,
which can be modified in the usual way:

``` r

prior <- bison_brm_get_prior(
  bison(edge_weight(node_1_id, node_2_id)) ~ age_diff + (1 | mm(node_1_id, node_2_id)),
  fit_edge,
  df_edge
)
prior$prior[1] <- "normal(0, 1)" # Prior on all population-level effects

bison_brm(
  bison(edge_weight(node_1_id, node_2_id)) ~ age_diff + (1 | mm(node_1_id, node_2_id)),
  fit_edge,
  df_edge,
  prior = prior
)
```

## More

Other brms features, such as different response families or smooth
terms, can also be used with
[`bison_brm()`](https://bisonr.org/reference/bison_brm.md). See the
[brms documentation](https://paulbuerkner.com/brms/) for details.
