# Draw samples of the edge list from the posterior

Returns randomly chosen posterior draws of all the edge weights. Each
draw is one plausible version of the network, consistent with the data.
Use it to repeat an analysis in other software on many versions of the
network, so the results include the uncertainty in the edge weights.

## Usage

``` r
draw_edgelist_samples(obj, num_draws, transform = FALSE)
```

## Arguments

- obj:

  A fitted `bison_model` object.

- num_draws:

  Number of posterior draws to return. Draws are chosen without
  replacement, unless more are requested than the model has.

- transform:

  If `TRUE`, edge weights are on the natural scale (probabilities or
  rates); if `FALSE` (the default, for backwards compatibility) they're
  on the link scale (log-odds or log rates).

## Value

A data frame with one row per dyad, columns `node_1` and `node_2`, and
one column per posterior draw (`draw.1`, `draw.2`, ...).

## See also

[`get_edgelist()`](https://bisonr.org/reference/get_edgelist.md) for a
summary of each edge weight,
[`bison_to_igraph()`](https://bisonr.org/reference/bison_to_igraph.md)
for draws as igraph networks, and
[`extract_metric()`](https://bisonr.org/reference/extract_metric.md) for
draws of network metrics.

## Examples

``` r
sim_data <- simulate_bison_model("binary", aggregated = TRUE)
fit_edge <- bison_model(
  (event | effort) ~ dyad(node_1_id, node_2_id),
  data = sim_data$df_sim,
  model_type = "binary_conjugate",
  priors = get_default_priors("binary_conjugate")
)
head(draw_edgelist_samples(fit_edge, num_draws = 3, transform = TRUE))
#>   node_1 node_2    draw.1     draw.2    draw.3
#> 1      1      2 0.2976116 0.03244854 0.6726602
#> 2      1      3 0.8843165 0.18874521 0.8969270
#> 3      1      4 0.1601845 0.43872967 0.4599058
#> 4      1      5 0.2104690 0.73684290 0.6061320
#> 5      1      6 0.8727320 0.46594499 0.3865661
#> 6      1      7 0.7500169 0.90475837 0.4474396
```
