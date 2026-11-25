# Simulates data

Simulates data

## Usage

``` r
simulate_bison_model(
  model_type,
  aggregated,
  location_effect = TRUE,
  age_diff_effect = TRUE
)
```

## Arguments

- model_type:

  Type of data to simulate, one of "binary", "count", and "duration".

- aggregated:

  TRUE/FALSE indicating whether to aggregate over dyads or not.

- location_effect:

  TRUE/FALSE indicating whether there is an effect of location on
  observations.

- age_diff_effect:

  TRUE/FALSE indicating whether there is an effect of age difference on
  observations.

## Value

A dataframe of the format accepted by edge weight models.
