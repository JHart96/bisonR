# Get the default priors for an edge model

Returns the default priors for a type of edge model, as a starting point
for choosing priors. A prior describes which values of a parameter are
plausible before seeing the data. The defaults may not suit every study,
so check them with
[`prior_check()`](https://bisonr.org/reference/prior_check.md) and
[`prior_predictive_check()`](https://bisonr.org/reference/prior_predictive_check.md)
and change them where needed before passing them to
[`bison_model()`](https://bisonr.org/reference/bison_model.md).

## Usage

``` r
get_default_priors(model_type)
```

## Arguments

- model_type:

  The type of model: `"binary"`, `"count"`, `"binary_conjugate"` or
  `"count_conjugate"`.

## Value

A named list of priors, in the format used by the `priors` argument of
[`bison_model()`](https://bisonr.org/reference/bison_model.md).

## Details

Each prior is a string naming a distribution and its parameters, such as
`"normal(0, 2.5)"`. The entries are:

- `edge`: the edge weights. For `"binary"` and `"count"` models this is
  a normal distribution on the link scale (log-odds for binary models,
  log rates for count models). With `partial_pooling = TRUE`, the second
  value is instead the scale of a half-normal prior on the spread of the
  edge weights. For `"binary_conjugate"` it is a beta distribution on
  the probability, and for `"count_conjugate"` a gamma distribution
  (shape, rate) on the rate.

- `fixed`: fixed effects, as a normal distribution on the link scale.

- `random_std`: the standard deviation of random effects, as a
  half-normal distribution with the given scale.

- `zero_prob`: in zero-inflated models, the probability that a row of
  data is a structural zero (a zero regardless of the edge weight), as a
  beta distribution.

Conjugate models only have an `edge` prior. Priors for terms that a
model doesn't include, such as `fixed` in a model without fixed effects,
have no effect.

## See also

[`prior_check()`](https://bisonr.org/reference/prior_check.md),
[`prior_predictive_check()`](https://bisonr.org/reference/prior_predictive_check.md),
[`bison_model()`](https://bisonr.org/reference/bison_model.md)

## Examples

``` r
priors <- get_default_priors("binary")
priors
#> $edge
#> [1] "normal(0, 2.5)"
#> 
#> $fixed
#> [1] "normal(0, 2.5)"
#> 
#> $random_std
#> [1] "half-normal(1)"
#> 
#> $zero_prob
#> [1] "beta(1, 1)"
#> 
# Priors are strings describing distributions, and can be changed:
priors$edge <- "normal(-1, 1.5)"
```
