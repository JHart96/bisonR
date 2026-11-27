# Simulate observation data with several types of edge

Simulates observations of social events in an undirected network of 15
nodes, where each edge belongs to one of several components (types of
relationship) with different mean edge weights. Use it to try out
[`bison_mixture()`](https://bisonr.org/reference/bison_mixture.md), or
to check that it recovers the true number of components.

## Usage

``` r
simulate_bison_model_mixture(model_type, num_components, component_weights)
```

## Arguments

- model_type:

  Type of data to simulate, `"binary"` or `"count"`.

- num_components:

  Number of components in the simulation (at least 2).

- component_weights:

  Vector giving the proportion of edges in each component, e.g.
  `c(0.5, 0.5)`. Must sum to 1.

## Value

A list with `df_sim`, the simulated observations, and `df_true`, the
true edge weights on the link scale. `df_sim` has columns `node_1_id`,
`node_2_id`, `event` and `effort` (observation effort). It also has a
`duration` column, identical to `effort`, for code written for earlier
versions; new code should use `effort`.

## Details

Component means are evenly spaced between -2.5 and 2.5 on the link scale
(log-odds for binary data, log rates for count data), and edge weights
within a component have a standard deviation of 0.5. Each dyad is
observed between 1 and 20 times, and the data are summed to one row per
dyad.

## See also

[`bison_mixture()`](https://bisonr.org/reference/bison_mixture.md),
[`simulate_bison_model()`](https://bisonr.org/reference/simulate_bison_model.md)

## Examples

``` r
sim_data <- simulate_bison_model_mixture(
  "binary", num_components = 2, component_weights = c(0.5, 0.5)
)
head(sim_data$df_sim)
#>   node_1_id node_2_id event effort duration
#> 1         1         2    15     17       17
#> 2         1         3    17     18       18
#> 3         1         4    19     19       19
#> 4         1         5     0     20       20
#> 5         1         6     0      2        2
#> 6         1         7     1     16       16
```
