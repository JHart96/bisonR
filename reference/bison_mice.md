# Prepare network data for multiple imputation with the mice package

Combines a regression data frame with posterior draws of a network
metric, producing one completed copy of the data for each draw. This is
the step [`bison_brm()`](https://bisonr.org/reference/bison_brm.md) uses
to carry the uncertainty in the network into a regression. Use
`bison_mice()` directly to fit the regression with software other than
brms, or to inspect the data that
[`bison_brm()`](https://bisonr.org/reference/bison_brm.md) would use.

## Usage

``` r
bison_mice(
  edgemodel_list,
  data_list,
  param_names,
  target_name,
  metric_name,
  num_draws = 100,
  z_score = FALSE,
  scale = c("none", "centre", "z")
)
```

## Arguments

- edgemodel_list:

  A fitted `bison_model` object, or a list of them for an analysis of
  several networks.

- data_list:

  A data frame of regression variables, or a list of data frames with
  one per edge model, in the same order. For edge and node metrics, each
  row refers to a dyad or node in the corresponding edge model. For
  global metrics, give one data frame with one row per edge model.

- param_names:

  Names of the columns in the data that identify the targets of the
  metric: two columns (e.g. `c("node_1_id", "node_2_id")`) for edge
  weights and one (e.g. `"node_id"`) for node metrics. The values must
  match the node names in the edge model. Not used for global metrics.

- target_name:

  The level of the metric: `"edge"`, `"node"` or `"global"`.

- metric_name:

  Name of the metric, without the target prefix, e.g. `"weight"` (the
  only edge metric), `"strength"`, `"degree[0.2]"` or `"cv"`. See
  [`extract_metric()`](https://bisonr.org/reference/extract_metric.md)
  for the metrics available.

- num_draws:

  Number of posterior draws to use, which is the number of imputed
  datasets. 100 is a reasonable value for a final analysis; small values
  (such as 5) are useful for testing code.

- z_score:

  If `TRUE`, standardise the network variable within each posterior
  draw. Equivalent to `scale = "z"`.

- scale:

  How to transform the network variable within each posterior draw:
  `"none"` (default), `"centre"` (subtract the mean) or `"z"` (subtract
  the mean and divide by the standard deviation). For edge and node
  metrics this is done across the rows of each data frame separately;
  global metrics are transformed across networks.

## Value

A `mids` object from the mice package, holding the original data (with
the network variable missing) and one completed dataset per posterior
draw.

## Details

**Multiple imputation** is a method for missing data: the analysis is
repeated on several plausible versions of the data and the results are
combined, so the final estimates include the uncertainty about the
missing values. bisonR treats a network metric as a variable whose true
values are unknown but described by the posterior distribution of the
edge model (the distribution of values consistent with the data and
priors). Each posterior draw, a single sample from that distribution,
gives one version of the data.

The [mice](https://cran.r-project.org/package=mice) package is an R
package for multiple imputation. It stores the versions of a dataset in
a `mids` object ("multiply imputed data set"). Regression packages that
support mice, such as brms
([`brms::brm_multiple()`](https://paulbuerkner.com/brms/reference/brm_multiple.html)),
and functions in mice itself
([with()](https://amices.org/mice/reference/with.mids.html) and
[`mice::pool()`](https://amices.org/mice/reference/pool.html)) can fit a
model to each version and combine the results. Use
[`mice::complete()`](https://tidyr.tidyverse.org/reference/complete.html)
to get the data frames back.

The network metric is added as a column named after the metric, such as
`bison_node_strength` or `bison_edge_weight`. A factor `bison_network`
gives the position of the edge model in `edgemodel_list` for each row.
Edge weights are on the link scale (log-odds for binary models, log
rates for count models), which is unbounded and usually better suited to
regression. Node and global metrics are calculated from edge weights on
the natural scale, as in
[`extract_metric()`](https://bisonr.org/reference/extract_metric.md).

## See also

[`bison_brm()`](https://bisonr.org/reference/bison_brm.md), which calls
this function and fits the regression;
[`extract_metric()`](https://bisonr.org/reference/extract_metric.md) for
posterior draws of a metric on its own.

## Examples

``` r
sim_data <- simulate_bison_model("binary", aggregated = TRUE)
df <- sim_data$df_sim
fit_edge <- bison_model(
  (event | effort) ~ dyad(node_1_id, node_2_id),
  data = df,
  model_type = "binary_conjugate",
  priors = get_default_priors("binary_conjugate")
)
# One imputed dataset per posterior draw of node strength
df_nodal <- data.frame(node = levels(df$node_1_id))
imputed <- bison_mice(fit_edge, df_nodal, "node", "node", "strength", num_draws = 5)
head(mice::complete(imputed, 1))
#>   node bison_network bison_node_strength
#> 1    1             1            3.211276
#> 2    2             1            3.559097
#> 3    3             1            4.185086
#> 4    4             1            5.222036
#> 5    5             1            4.413708
#> 6    6             1            3.411393
```
