# Get the mean edge weight of each mixture component

For a mixture with a given number of components, returns the posterior
median and credible interval of each component's mean edge weight, on
the natural scale (probabilities for binary models, rates for count
models). This describes how strong each type of relationship is. A
credible interval contains the true value with the given probability,
according to the model.

## Usage

``` r
get_component_means(object, num_components, ci = 0.9)
```

## Arguments

- object:

  A `bison_mixture` object from
  [`bison_mixture()`](https://bisonr.org/reference/bison_mixture.md).

- num_components:

  The number of components in the mixture to use, between 1 and the
  `num_components` used in
  [`bison_mixture()`](https://bisonr.org/reference/bison_mixture.md).

- ci:

  Width of the credible interval, e.g. `0.9` for a 90% interval.

## Value

A data frame with one row per component (`Component 1`, ...), in
increasing order of their mean, and columns for the posterior median
(`50%`) and the lower and upper bounds of the credible interval (named
after their quantiles, e.g. `5%` and `95%`).

## See also

[`bison_mixture()`](https://bisonr.org/reference/bison_mixture.md),
[`plot_mixture()`](https://bisonr.org/reference/plot_mixture.md)

## Examples

``` r
sim_data <- simulate_bison_model_mixture(
  "count", num_components = 2, component_weights = c(0.5, 0.5)
)
fit_edge <- bison_model(
  (event | effort) ~ dyad(node_1_id, node_2_id),
  data = sim_data$df_sim,
  model_type = "count_conjugate",
  priors = get_default_priors("count_conjugate")
)
fit_mixture <- bison_mixture(fit_edge, num_components = 3, verbose = FALSE)
get_component_means(fit_mixture, 2)
#>                     50%          5%         95%
#> Component 1  0.07992664  0.07158064  0.08945183
#> Component 2 11.80658784 11.66211708 11.93652237
```
