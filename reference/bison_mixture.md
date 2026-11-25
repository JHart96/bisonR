# Fit a dyadic mixture model

Fit a dyadic mixture model

## Usage

``` r
bison_mixture(edgemodel, num_components = 5, criterion = "ICL", verbose = TRUE)
```

## Arguments

- edgemodel:

  A fitted edge weight model.

- num_components:

  Maximum number of components to fit.

- criterion:

  "ICL" or "BIC" describing the criterion to use for model comparison.

- verbose:

  TRUE/FALSE indicating whether to output messages and progress when
  fitting.

## Value

A bison_mixture object
