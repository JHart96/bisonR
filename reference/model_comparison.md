# Compare models with leave-one-out cross-validation

Compares fitted edge models by their expected predictive accuracy,
estimated with Pareto-smoothed importance sampling leave-one-out
cross-validation (PSIS-LOO) from the loo package, and returns model
weights. Leave-one-out cross-validation measures how well a model
predicts each observation when that observation is left out of the fit;
PSIS-LOO approximates this without refitting the model. Use it to
compare models with different terms, for example to test whether edge
weights differ between dyads by comparing a model with a `dyad()` term
against a null model with only an intercept.

## Usage

``` r
model_comparison(models, ...)
```

## Arguments

- models:

  A list of two or more models fitted with
  [`bison_model()`](https://bisonr.org/reference/bison_model.md),
  ideally named. All models must be fitted to the same data (the same
  rows, in the same order).

- ...:

  Additional arguments passed to
  [`loo::loo_model_weights()`](https://mc-stan.org/loo/reference/loo_model_weights.html),
  such as `method = "pseudobma"`.

## Value

The model weights returned by
[`loo::loo_model_weights()`](https://mc-stan.org/loo/reference/loo_model_weights.html),
one per model, named after the elements of `models`. By default these
are stacking weights.

## Details

The model weights describe how much each model contributes to the best
predictions: a model with a weight close to 1 predicts the data better
than the others. loo may warn about high Pareto k values, which mean the
approximation is unreliable for some observations; see
[`loo::loo()`](https://mc-stan.org/loo/reference/loo.html).

If the models have different numbers of posterior draws (for example
because a chain failed), the same number of draws is taken at random
from each model.

## See also

[`bison_model()`](https://bisonr.org/reference/bison_model.md)

## Examples

``` r
# Test for non-random edge weights by comparing a model with edge weights
# against a null model where every dyad has the same rate of events.
sim_data <- simulate_bison_model("binary", aggregated = FALSE)
priors <- get_default_priors("binary")
fit_edge <- bison_model((event | effort) ~ dyad(node_1_id, node_2_id),
                        data = sim_data$df_sim, model_type = "binary", priors = priors)
#> Running MCMC with 4 parallel chains...
#> 
#> Chain 1 finished in 0.5 seconds.
#> Chain 2 finished in 0.5 seconds.
#> Chain 3 finished in 0.5 seconds.
#> Chain 4 finished in 0.5 seconds.
#> 
#> All 4 chains finished successfully.
#> Mean chain execution time: 0.5 seconds.
#> Total execution time: 0.7 seconds.
#> 
fit_null <- bison_model((event | effort) ~ 1,
                        data = sim_data$df_sim, model_type = "binary", priors = priors)
#> Running MCMC with 4 parallel chains...
#> 
#> Chain 1 finished in 0.3 seconds.
#> Chain 2 finished in 0.4 seconds.
#> Chain 3 finished in 0.3 seconds.
#> Chain 4 finished in 0.3 seconds.
#> 
#> All 4 chains finished successfully.
#> Mean chain execution time: 0.3 seconds.
#> Total execution time: 0.6 seconds.
#> 
model_comparison(list(non_random = fit_edge, random = fit_null))
#> Warning: Some Pareto k diagnostic values are too high. See help('pareto-k-diagnostic') for details.
#> Method: stacking
#> ------
#>            weight
#> non_random 0.543 
#> random     0.457 
```
