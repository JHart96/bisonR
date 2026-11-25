# Run leave-one-out model comparison

Run leave-one-out model comparison

## Usage

``` r
model_comparison(models, ...)
```

## Arguments

- models:

  List of fitted models to compare

- ...:

  Additional arguments to be passed to
  [`loo::loo_model_weights`](https://mc-stan.org/loo/reference/loo_model_weights.html)

## Value

The output of the `loo_compare` function for the fitted models.
