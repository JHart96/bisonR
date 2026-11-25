# Generate a mice mids object

Generate a mice mids object

## Usage

``` r
bison_mice(
  edgemodel_list,
  data_list,
  param_names,
  target_name,
  metric_name,
  num_draws = 100,
  z_score = FALSE
)
```

## Arguments

- edgemodel_list:

  A bisonR edge model (or a list of them)

- data_list:

  A dataframe of regression variables compatible with a brms model (or a
  list of them).

- param_names:

  List of column names for network features (e.g. node_1_id, node_2_id
  or just node_id)

- target_name:

  Name of network feature (edge, node, global)

- metric_name:

  Name of network metric to calculate

- num_draws:

  Number of draws from the network posterior

- z_score:

  Whether to Z-score bison variable or not.

## Value

A mice mids object
