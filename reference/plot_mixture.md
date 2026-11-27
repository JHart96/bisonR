# Plot the components of an edge mixture model

Shows the distribution of posterior median edge weights (on the natural
scale) as a histogram, with each edge coloured by the component it most
probably belongs to. The posterior median and credible interval of each
component's mean are shown as a vertical line and shaded band. Use it to
see how well the components separate the edge weights.

## Usage

``` r
plot_mixture(object, num_components = NULL, ci = 0.9, breaks = 30)
```

## Arguments

- object:

  A `bison_mixture` object from
  [`bison_mixture()`](https://bisonr.org/reference/bison_mixture.md).

- num_components:

  Number of components to show. Defaults to the number with the highest
  posterior probability.

- ci:

  Width of the credible interval for the component means.

- breaks:

  Number of histogram bins (passed to
  [`graphics::hist()`](https://rdrr.io/r/graphics/hist.html)).

## Value

The data used for the plot, invisibly: a data frame with each dyad, its
posterior median edge weight and its most probable component.

## See also

[`bison_mixture()`](https://bisonr.org/reference/bison_mixture.md),
[`get_component_means()`](https://bisonr.org/reference/get_component_means.md)

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
plot_mixture(fit_mixture)
```
