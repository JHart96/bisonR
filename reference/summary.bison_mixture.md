# Summarise an edge mixture model

Summarises the results of
[`bison_mixture()`](https://bisonr.org/reference/bison_mixture.md).
Printing the summary shows the most probable number of components and
its probability, the probability of each number of components, and, for
the most probable number, the posterior median of each component mean
(on the natural scale) and the probability of the first few edges
belonging to each component.

## Usage

``` r
# S3 method for class 'bison_mixture'
summary(object, ...)
```

## Arguments

- object:

  A `bison_mixture` object from
  [`bison_mixture()`](https://bisonr.org/reference/bison_mixture.md).

- ...:

  Not used.

## Value

A `summary.bison_mixture` object, which is printed as described above.

## See also

[`get_network_component_probabilities()`](https://bisonr.org/reference/get_network_component_probabilities.md),
[`get_edge_component_probabilities()`](https://bisonr.org/reference/get_edge_component_probabilities.md),
[`get_component_means()`](https://bisonr.org/reference/get_component_means.md)
for the same results as data frames.

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
summary(fit_mixture)
#> === Fitted dyadic mixture model ===
#> Maximum number of components: 3
#> Most probable number of components: 2
#> Posterior probability: 100%
#> === Component probabilities ===
#>        1 2 3
#> P(K=k) 0 1 0
#> === Component means (K = 2) ===
#>            1      2
#> median 0.093 11.778
#> === Edge component probabilities (K = 2) ===
#>         1 2
#> 1 <-> 2 1 0
#> 1 <-> 3 0 1
#> 1 <-> 4 1 0
#> 1 <-> 5 1 0
#> 1 <-> 6 0 1
#> 1 <-> 7 0 1
#> ... (99 more rows; see get_edge_component_probabilities())
```
