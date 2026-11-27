# Simulate observation data for an edge model

Simulates observations of social events in an undirected network with
known edge weights. Use it to try out bisonR functions, to learn the
format of the data, or to check that a model recovers the true values
before using it on real data.

## Usage

``` r
simulate_bison_model(
  model_type,
  aggregated,
  location_effect = TRUE,
  age_diff_effect = TRUE,
  num_nodes = 10,
  max_obs = 5
)
```

## Arguments

- model_type:

  Type of data to simulate, `"binary"` (whether an event occurred in
  each sampling period) or `"count"` (the number of events).

- aggregated:

  If `TRUE`, observations are summed to one row per dyad, with the
  location column removed. If `FALSE`, there is one row per observation.

- location_effect:

  Whether location affects the probability or rate of events.

- age_diff_effect:

  Whether the age difference between the two nodes affects the
  probability or rate of events (by 0.2 per unit of age difference on
  the link scale).

- num_nodes:

  Number of nodes in the simulated network.

- max_obs:

  Maximum number of observations per dyad. Each dyad gets between 1 and
  `max_obs`.

## Value

A list with `df_sim`, the simulated observations, and `df_true`, the
true edge weights on the link scale. `df_sim` has columns `node_1_id`,
`node_2_id`, `event` and `effort` (observation effort), and predictors.
It also has a `duration` column, identical to `effort`, for code written
for earlier versions; new code should use `effort`. The true edge
weights in `df_true` don't include the effects of age difference or
location.

## Details

Each dyad's true edge weight is drawn from a standard normal
distribution on the link scale (log-odds for binary data, log rates for
count data). Each node gets an age, and each observation takes place at
one of five locations, each with its own effect. Binary data have one
sampling period (effort 1) per observation; count data have an
observation effort drawn uniformly between 1 and `max_obs`.

## See also

[`simulate_bison_model_mixture()`](https://bisonr.org/reference/simulate_bison_model_mixture.md),
[`bison_model()`](https://bisonr.org/reference/bison_model.md)

## Examples

``` r
sim_data <- simulate_bison_model("binary", aggregated = FALSE)
head(sim_data$df_sim)  # Observations
#>   node_1_id node_2_id event effort age_diff    age_1    age_2 location duration
#> 1         1         2     1      1 2.946436 23.37606 20.42962        4        1
#> 2         1         2     1      1 2.946436 23.37606 20.42962        4        1
#> 3         1         2     1      1 2.946436 23.37606 20.42962        3        1
#> 4         1         2     1      1 2.946436 23.37606 20.42962        1        1
#> 5         1         2     0      1 2.946436 23.37606 20.42962        3        1
#> 6         1         3     0      1 2.688601 23.37606 20.68746        5        1
head(sim_data$df_true) # True edge weights (on the logit scale)
#>   node_1_id node_2_id edge_weight age_diff    age_1    age_2
#> 1         1         2   1.3930343 2.946436 23.37606 20.42962
#> 2         1         3  -0.5495812 2.688601 23.37606 20.68746
#> 3         1         4   1.1436925 2.818835 23.37606 20.55722
#> 4         1         5  -0.7115225 9.347903 23.37606 14.02815
#> 5         1         6  -0.6320064 2.483173 23.37606 20.89288
#> 6         1         7  -0.0657349 7.461969 23.37606 15.91409
```
