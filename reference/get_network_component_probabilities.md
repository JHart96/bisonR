# Get the probability of each number of components in the network

Returns the posterior probability of each number of components (types of
relationship) considered by
[`bison_mixture()`](https://bisonr.org/reference/bison_mixture.md), from
1 up to the maximum. Use it to judge how many types of edge the network
contains.

## Usage

``` r
get_network_component_probabilities(object)
```

## Arguments

- object:

  A `bison_mixture` object from
  [`bison_mixture()`](https://bisonr.org/reference/bison_mixture.md).

## Value

A data frame with columns `num_components` and `probability`, with one
row for each number of components. The probabilities sum to 1.

## See also

[`bison_mixture()`](https://bisonr.org/reference/bison_mixture.md),
[`get_edge_component_probabilities()`](https://bisonr.org/reference/get_edge_component_probabilities.md)

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
get_network_component_probabilities(fit_mixture)
#>   num_components  probability
#> 1              1 1.413129e-29
#> 2              2 7.843249e-01
#> 3              3 2.156751e-01
```
