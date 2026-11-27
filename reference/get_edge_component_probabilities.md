# Get the probability of each edge belonging to each component

For a mixture with a given number of components, returns the posterior
probability that each edge belongs to each component (type of
relationship). The probabilities include the uncertainty in both the
edge weights and the mixture. Use it to classify edges, for example as
weak or strong bonds.

## Usage

``` r
get_edge_component_probabilities(object, num_components)
```

## Arguments

- object:

  A `bison_mixture` object from
  [`bison_mixture()`](https://bisonr.org/reference/bison_mixture.md).

- num_components:

  The number of components in the mixture to use, between 1 and the
  `num_components` used in
  [`bison_mixture()`](https://bisonr.org/reference/bison_mixture.md).
  This is often the most probable number, from
  [`get_network_component_probabilities()`](https://bisonr.org/reference/get_network_component_probabilities.md).

## Value

A data frame with one row per dyad, columns `node_1` and `node_2`, and
one column per component (`P(component 1)`, ...). Each row sums to 1.
Components are in increasing order of their mean.

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
head(get_edge_component_probabilities(fit_mixture, 2))
#>   node_1 node_2 P(component 1) P(component 2)
#> 1      1      2              0              1
#> 2      1      3              1              0
#> 3      1      4              0              1
#> 4      1      5              1              0
#> 5      1      6              1              0
#> 6      1      7              0              1
```
