# Convert a fitted edge model into a list of igraph networks

Each network is one draw from the posterior distribution of networks,
with edge weights on the natural scale (probabilities or rates). The
vertices are always all the nodes in the model, in the same order as the
model's node levels, so results can be matched back to node names.

## Usage

``` r
bison_to_igraph(obj, num_draws = 1000)
```

## Arguments

- obj:

  A fitted `bison_model` object.

- num_draws:

  Number of networks to draw from the posterior. Draws are chosen at
  random, without replacement unless more are requested than the model
  has.

## Value

A list of `igraph` objects, with edge weights in the `weight` edge
attribute. Edges are only included for dyads in the model.

## Details

Use it to calculate network measures that
[`extract_metric()`](https://bisonr.org/reference/extract_metric.md)
doesn't provide, with functions from the igraph package or other
software that reads igraph objects. Calculating a measure for every
network gives its posterior distribution, which includes the uncertainty
in the edge weights.

## See also

[`extract_metric()`](https://bisonr.org/reference/extract_metric.md),
[`draw_edgelist_samples()`](https://bisonr.org/reference/draw_edgelist_samples.md)

## Examples

``` r
sim_data <- simulate_bison_model("binary", aggregated = TRUE)
fit_edge <- bison_model(
  (event | effort) ~ dyad(node_1_id, node_2_id),
  data = sim_data$df_sim,
  model_type = "binary_conjugate",
  priors = get_default_priors("binary_conjugate")
)
nets <- bison_to_igraph(fit_edge, num_draws = 10)
# Posterior distribution of the strength of each node
strengths <- sapply(nets, igraph::strength)
```
