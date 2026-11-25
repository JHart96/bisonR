# Return default priors for a given model

Return default priors for a given model

## Usage

``` r
get_default_priors(model_type)
```

## Arguments

- model_type:

  A character specifying the type of model to retrieve priors for
  ("binary", "count", "duration", etc)

## Value

A list of default priors that can be used in the
[`bison_model()`](https://bisonr.org/reference/bison_model.md) function.
