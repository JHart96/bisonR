# Posterior predictive checks

Compares data simulated from a fitted model with the observed data, to
check whether the model describes the data well. If the model fits, the
observed data should look like a typical simulated dataset. Use it after
[`bison_model()`](https://bisonr.org/reference/bison_model.md), before
interpreting the results.

## Usage

``` r
plot_predictions(obj, ...)

# S3 method for class 'bison_model'
plot_predictions(
  obj,
  num_draws = 20,
  type = c("density", "point"),
  draw_data = TRUE,
  ...
)
```

## Arguments

- obj:

  A fitted `bison_model` object.

- ...:

  Arguments passed to methods.

- num_draws:

  Number of simulated datasets to compare against.

- type:

  One or more of `"density"`, `"stats"` and `"point"`.

- draw_data:

  Whether to show the observed data (set to `FALSE` for prior predictive
  checks).

## Value

`NULL`, invisibly. The function is called for its plots.

## Details

Three types of check are available, chosen with `type`:

- `"density"`: the distribution of observed event *rates* per dyad
  (events divided by observation effort), compared with the same
  distribution in `num_draws` simulated datasets. Rates are used rather
  than event counts because counts largely reflect observation effort.

- `"stats"`: summary statistics of the per-dyad rates in each simulated
  dataset (proportion of dyads with no events, mean, standard deviation
  and maximum), with the observed value shown as a vertical line. These
  can show excess zeros or over- or under-dispersion between dyads.

- `"point"`: posterior medians and intervals of the edge weights plotted
  against simple point estimates (events / effort). This is not a
  predictive check, but shows how estimates are pulled towards the prior
  when there are few data.

The default shows the `"density"` and `"point"` plots.

## See also

[`prior_predictive_check()`](https://bisonr.org/reference/prior_predictive_check.md)
for the same plots using only the priors, and
[`plot_trace()`](https://bisonr.org/reference/plot_trace.md).

## Examples

``` r
sim_data <- simulate_bison_model("binary", aggregated = TRUE)
fit_edge <- bison_model(
  (event | effort) ~ dyad(node_1_id, node_2_id),
  data = sim_data$df_sim,
  model_type = "binary_conjugate",
  priors = get_default_priors("binary_conjugate")
)
plot_predictions(fit_edge, type = "density")

plot_predictions(fit_edge, type = "stats")
```
