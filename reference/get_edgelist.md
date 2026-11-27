# Get the edge list with uncertainty

Returns the estimated edge weight of every dyad as a data frame: the
posterior median and a credible interval, which contains the true value
with the given probability according to the model. Use it to report or
export edge weights, or to work with them in other software.

## Usage

``` r
get_edgelist(obj, ci = 0.9, transform = TRUE)
```

## Arguments

- obj:

  A fitted `bison_model` object.

- ci:

  Width of the credible interval, e.g. `0.9` for a 90% interval.

- transform:

  If `TRUE` (the default), edge weights are on the natural scale
  (probabilities for binary models, rates for count models). If `FALSE`,
  they're on the link scale used inside the model (log-odds for binary
  models, log rates for count models). In models with an intercept or
  other effects, the natural-scale values are for an observation where
  those effects are 0.

## Value

A data frame with one row per dyad and columns `node_1`, `node_2`,
`median`, and the lower and upper bounds of the credible interval (named
after their quantiles, e.g. `5%` and `95%`). Values are rounded to three
decimal places.

## See also

[`draw_edgelist_samples()`](https://bisonr.org/reference/draw_edgelist_samples.md)
for individual posterior draws of the edge list,
[`summary.bison_model()`](https://bisonr.org/reference/summary.bison_model.md),
[`plot_network()`](https://bisonr.org/reference/plot_network.md).

## Examples

``` r
sim_data <- simulate_bison_model("binary", aggregated = TRUE)
fit_edge <- bison_model(
  (event | effort) ~ dyad(node_1_id, node_2_id),
  data = sim_data$df_sim,
  model_type = "binary_conjugate",
  priors = get_default_priors("binary_conjugate")
)
head(get_edgelist(fit_edge))
#>   node_1 node_2 median    5%   95%
#> 1      1      2  0.711 0.216 0.975
#> 2      1      3  0.498 0.194 0.802
#> 3      1      4  0.383 0.100 0.754
#> 4      1      5  0.733 0.422 0.935
#> 5      1      6  0.684 0.350 0.929
#> 6      1      7  0.585 0.274 0.853
```
