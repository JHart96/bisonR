# Get the estimates of other effects in an edge model

Edge weight models can include effects other than the edge weights:
fixed effects, random effects, partial pooling and zero inflation (see
[`bison_model()`](https://bisonr.org/reference/bison_model.md)). This
function returns their posterior medians and credible intervals in a
data frame. A fixed effect is a single estimated effect of a variable,
such as age difference. A random effect estimates a separate value for
each level of a grouping variable, such as location, with the values
assumed to come from a shared normal distribution.

## Usage

``` r
get_effects(obj, ci = 0.9)
```

## Arguments

- obj:

  A fitted `bison_model` object, fitted with Stan (conjugate models only
  estimate edge weights).

- ci:

  Width of the credible interval, e.g. `0.9` for a 90% interval.

## Value

A data frame with columns `effect`, `type`, `median`, and the lower and
upper bounds of the credible interval (named after their quantiles, e.g.
`5%` and `95%`). It has no rows if the model has no effects other than
edge weights.

## Details

The `type` column says what each row is:

- `"fixed"`: a fixed effect, such as `age` for a numeric column, or
  `sex = F` for one level of a factor. These are on the same scale as
  the edge weights before they are transformed: log-odds for binary
  models and log rates for count models.

- `"random"`: the random effect for one level of a group, such as
  `location = 2`, on the same scale as fixed effects.

- `"random sd"`: the standard deviation of the random effects in a
  group.

- `"pooling sd"`: with `partial_pooling = TRUE`, the standard deviation
  of the edge weights.

- `"zero inflation"`: with `zero_inflated = TRUE`, the probability that
  an observation is a structural zero.

## See also

[`get_edgelist()`](https://bisonr.org/reference/get_edgelist.md) for the
edge weights, and
[`summary.bison_model()`](https://bisonr.org/reference/summary.bison_model.md),
which shows both.

## Examples

``` r
sim_data <- simulate_bison_model("binary", aggregated = FALSE)
df <- sim_data$df_sim
df$location <- factor(df$location)
fit_edge <- bison_model(
  (event | effort) ~ dyad(node_1_id, node_2_id) + age_diff + (1 | location),
  data = df,
  model_type = "binary",
  priors = get_default_priors("binary")
)
#> Running MCMC with 4 parallel chains...
#> 
#> Chain 1 finished in 1.0 seconds.
#> Chain 2 finished in 1.0 seconds.
#> Chain 3 finished in 1.0 seconds.
#> Chain 4 finished in 1.0 seconds.
#> 
#> All 4 chains finished successfully.
#> Mean chain execution time: 1.0 seconds.
#> Total execution time: 1.3 seconds.
#> 
get_effects(fit_edge)
#>         effect      type median     5%    95%
#> 1     age_diff     fixed  0.267  0.087  0.444
#> 2 location = 1    random -0.355 -1.414  0.583
#> 3 location = 2    random  0.043 -0.964  1.043
#> 4 location = 3    random  1.291  0.118  2.722
#> 5 location = 4    random -1.527 -2.649 -0.423
#> 6 location = 5    random -0.052 -1.107  0.997
#> 7 sd(location) random sd  1.099  0.507  1.918
```
