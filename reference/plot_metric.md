# Plot the posterior distribution of a network metric

Plots the posterior distribution of a network metric from
[`extract_metric()`](https://bisonr.org/reference/extract_metric.md),
which shows the most plausible values of the metric and how uncertain it
is. Continuous metrics are shown as smoothed densities; metrics that
only take whole-number values (such as degree) are shown as bar charts.
A matrix with several columns (for example, one per node) gives one plot
per column.

## Usage

``` r
plot_metric(x, metric_name = "Network Metric", ...)
```

## Arguments

- x:

  A matrix of metric samples from
  [`extract_metric()`](https://bisonr.org/reference/extract_metric.md).

- metric_name:

  Label for the x axis. With several columns, the column name or number
  is added to it.

- ...:

  Additional arguments passed to the plotting function
  ([`graphics::barplot()`](https://rdrr.io/r/graphics/barplot.html) or
  [`graphics::curve()`](https://rdrr.io/r/graphics/curve.html)), such as
  `xlim`.

## Value

`NULL`, invisibly. The function is called for its plots.

## See also

[`extract_metric()`](https://bisonr.org/reference/extract_metric.md)

## Examples

``` r
sim_data <- simulate_bison_model("binary", aggregated = TRUE)
fit_edge <- bison_model(
  (event | effort) ~ dyad(node_1_id, node_2_id),
  data = sim_data$df_sim,
  model_type = "binary_conjugate",
  priors = get_default_priors("binary_conjugate")
)
plot_metric(extract_metric(fit_edge, "global_cv", num_draws = 200), "CV of edge weights")
```
