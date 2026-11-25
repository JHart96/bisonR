# Fit a brms model on BISoN network data

Fit a brms model on BISoN network data

## Usage

``` r
bison_brm(
  formula,
  edgemodel_list,
  data_list,
  num_draws = 100,
  z_score = FALSE,
  cores = 4,
  chains = 2,
  ...
)
```

## Arguments

- formula:

  A formula compatible with a brms model

- edgemodel_list:

  A bisonR edge model (or a list of them)

- data_list:

  A dataframe of regression variables compatible with a brms model (or a
  list of them)

- num_draws:

  Number of draws from the network posterior to use when fitting model

- z_score:

  Whether to Z-score bison variable or not.

- cores:

  Number of computational cores to use.

- chains:

  Number of chains to run per model.

- ...:

  Additional arguments to be passed to brm(), such as family or priors

## Value

A fitted brms model
