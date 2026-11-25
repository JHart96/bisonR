# Get priors for BISoN brm model

Wrapper for brms::get_prior() for use with bison_brm.

## Usage

``` r
bison_brm_get_prior(formula, edgemodel_list, data_list, ...)
```

## Arguments

- formula:

  A formula compatible with a brms model

- edgemodel_list:

  A bisonR edge model (or a list of them)

- data_list:

  A dataframe of regression variables compatible with a brms model (or a
  list of them)

- ...:

  Additional parameters to be passed to brms::get_prior()

## Value

The output from brms::get_prior()
