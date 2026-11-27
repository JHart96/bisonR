# Summarise an edge weight model

Summarises a fitted edge model: the posterior median and credible
interval of each edge weight, and of any other effects in the model. A
credible interval contains the true value with the given probability,
according to the model. Use
[`summary()`](https://rdrr.io/r/base/summary.html) to check the
estimates after fitting with
[`bison_model()`](https://bisonr.org/reference/bison_model.md).

## Usage

``` r
# S3 method for class 'bison_model'
summary(object, ci = 0.9, transform = TRUE, ...)
```

## Arguments

- object:

  A fitted `bison_model` object.

- ci:

  Width of the credible interval, e.g. `0.9` for a 90% interval.

- transform:

  If `TRUE` (the default), edge weights are reported on the natural
  scale (probabilities for binary models, rates for count models). If
  `FALSE`, they're on the link scale used inside the model (log-odds for
  binary models, log rates for count models). In models with an
  intercept or other effects, the natural-scale values are for an
  observation where those effects are 0.

- ...:

  Not used.

## Value

A `summary.bison_model` object. Printing it shows, for each dyad, the
posterior median edge weight and the lower and upper bounds of the
credible interval, followed by any other effects in the model (see
[`get_effects()`](https://bisonr.org/reference/get_effects.md)).

## See also

[`get_edgelist()`](https://bisonr.org/reference/get_edgelist.md) to get
the edge weights as a data frame.

## Examples

``` r
sim_data <- simulate_bison_model("binary", aggregated = TRUE)
fit_edge <- bison_model(
  (event | effort) ~ dyad(node_1_id, node_2_id),
  data = sim_data$df_sim,
  model_type = "binary_conjugate",
  priors = get_default_priors("binary_conjugate")
)
summary(fit_edge, ci = 0.95)
#> === Fitted BISoN edge model ===
#> Data type: binary
#> Formula: (event | effort) ~ dyad(node_1_id, node_2_id)
#> Number of nodes: 10
#> Number of dyads: 45
#> Directed: FALSE
#> === Edge list summary ===
#>          median  2.5% 97.5%
#> 1 <-> 2   0.708 0.159 0.989
#> 1 <-> 3   0.703 0.150 0.988
#> 1 <-> 4   0.866 0.476 0.995
#> 1 <-> 5   0.495 0.103 0.907
#> 1 <-> 6   0.868 0.472 0.995
#> 1 <-> 7   0.869 0.472 0.994
#> 1 <-> 8   0.711 0.141 0.988
#> 1 <-> 9   0.736 0.360 0.960
#> 1 <-> 10  0.287 0.012 0.854
#> 2 <-> 3   0.423 0.122 0.783
#> 2 <-> 4   0.687 0.277 0.951
#> 2 <-> 5   0.313 0.048 0.715
#> 2 <-> 6   0.307 0.047 0.720
#> 2 <-> 7   0.157 0.007 0.611
#> 2 <-> 8   0.582 0.225 0.887
#> 2 <-> 9   0.720 0.148 0.987
#> 2 <-> 10  0.497 0.096 0.903
#> 3 <-> 4   0.499 0.099 0.916
#> 3 <-> 5   0.294 0.013 0.846
#> 3 <-> 6   0.727 0.345 0.954
#> 3 <-> 7   0.612 0.195 0.936
#> 3 <-> 8   0.292 0.013 0.840
#> 3 <-> 9   0.613 0.199 0.933
#> 3 <-> 10  0.707 0.164 0.988
#> 4 <-> 5   0.839 0.410 0.993
#> 4 <-> 6   0.490 0.146 0.855
#> 4 <-> 7   0.736 0.355 0.956
#> 4 <-> 8   0.843 0.401 0.993
#> 4 <-> 9   0.703 0.154 0.989
#> 4 <-> 10  0.575 0.219 0.875
#> 5 <-> 6   0.703 0.161 0.987
#> 5 <-> 7   0.800 0.272 0.991
#> 5 <-> 8   0.714 0.154 0.988
#> 5 <-> 9   0.497 0.090 0.911
#> 5 <-> 10  0.489 0.091 0.896
#> 6 <-> 7   0.420 0.114 0.776
#> 6 <-> 8   0.869 0.485 0.995
#> 6 <-> 9   0.682 0.290 0.945
#> 6 <-> 10  0.496 0.094 0.905
#> 7 <-> 8   0.708 0.184 0.986
#> 7 <-> 9   0.203 0.009 0.701
#> 7 <-> 10  0.112 0.005 0.472
#> 8 <-> 9   0.707 0.141 0.988
#> 8 <-> 10  0.206 0.008 0.724
#> 9 <-> 10  0.264 0.043 0.650
```
