# Prior checks

Prior checks

## Usage

``` r
prior_check(priors, model_type, type = "density")
```

## Arguments

- priors:

  List of priors for a model, can be retrieved using
  `get_default_prior()`.

- model_type:

  Type of model the priors will be used for (same as the argument for
  `get_default_prior()`).

- type:

  Type of prior check to run, `"value"` or `"prediction"`. Details
  below.

## Details

The parameter `type` determines what type of prior check to run.
`type="value"` will plot the prior probability over the parameter value
on the original scale. `type="predictive"` will run a prior predictive
plot, where predictions from the model are generated using only prior
probabilities (the model when not updated from the data).
