# Simulate data for an edge model with mixture components

Simulate data for an edge model with mixture components

## Usage

``` r
simulate_bison_model_mixture(model_type, num_components, component_weights)
```

## Arguments

- model_type:

  Type of data to simulate, either "binary" or "count".

- num_components:

  Number of mixture components in the simulation (\>=2).

- component_weights:

  Vector of component weightings for each component. Must sum to 1.

## Value

A dataframe of the format accepted by edge weight models.
