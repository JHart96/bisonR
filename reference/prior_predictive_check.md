# Run prior predictive check for a model specification

Run prior predictive check for a model specification

## Usage

``` r
prior_predictive_check(
  formula,
  data,
  model_type,
  plot_type = "density",
  options = NULL,
  priors = NULL
)
```

## Arguments

- formula:

  Formula for the model

- data:

  Data for the model

- model_type:

  Type of model to fit

- plot_type:

  Type of predictive plot to show ("density" or "marginal")

- options:

  List of additional information needed for the predictive plot. Check
  the examples for full details.

- priors:

  List of priors for the model
