# Fit an edge mixture model

Looks for clusters (components) in the edge weights, which can be
interpreted as different types of social relationship, such as weak and
strong bonds. Use it after
[`bison_model()`](https://bisonr.org/reference/bison_model.md) to ask
whether the network contains distinct types of edge, and which edges
belong to each type.

## Usage

``` r
bison_mixture(
  edgemodel,
  num_components = 5,
  criterion = c("ICL", "BIC"),
  verbose = TRUE,
  num_draws = 200
)
```

## Arguments

- edgemodel:

  A fitted `bison_model` object.

- num_components:

  Maximum number of components to consider.

- criterion:

  `"ICL"` (default) or `"BIC"`: the information criterion used to weigh
  the different numbers of components. BIC (Bayesian information
  criterion) measures how well the mixture fits the edge weights. ICL
  (integrated completed likelihood) also penalises components that
  overlap, so it favours clearly separated clusters.

- verbose:

  Whether to show a progress bar.

- num_draws:

  Number of posterior draws to fit mixtures to (at most the number
  available). More draws give more precise probabilities but take
  longer.

## Value

A `bison_mixture` object. Use
[`summary()`](https://rdrr.io/r/base/summary.html) on it for an
overview,
[`get_network_component_probabilities()`](https://bisonr.org/reference/get_network_component_probabilities.md)
for the probability of each number of components,
[`get_edge_component_probabilities()`](https://bisonr.org/reference/get_edge_component_probabilities.md)
for the component of each edge,
[`get_component_means()`](https://bisonr.org/reference/get_component_means.md)
for the component means, and
[`plot_mixture()`](https://bisonr.org/reference/plot_mixture.md) to plot
them. Its `draw_ids` element gives the rows of `edgemodel$edge_samples`
used for each mixture draw, so component assignments can be matched to
the edge weights of the same draw.

## Details

A Gaussian mixture model describes the edge weights as coming from
several normal distributions (the components), each with its own mean
and spread. It is fitted with the mclust package to the edge weights (on
the link scale) in each of `num_draws` posterior draws, for every number
of components from 1 to `num_components`. In each draw, the numbers of
components are weighted by an information criterion, a measure of fit
that penalises complexity. Averaging over draws gives the probability of
each number of components, and of each edge belonging to each component,
taking into account the uncertainty in the edge weights.

Components are labelled in increasing order of their mean, so component
1 is always the weakest type of edge.

## See also

[`get_network_component_probabilities()`](https://bisonr.org/reference/get_network_component_probabilities.md),
[`get_edge_component_probabilities()`](https://bisonr.org/reference/get_edge_component_probabilities.md),
[`get_component_means()`](https://bisonr.org/reference/get_component_means.md),
[`plot_mixture()`](https://bisonr.org/reference/plot_mixture.md)

## Examples

``` r
sim_data <- simulate_bison_model_mixture(
  "binary", num_components = 2, component_weights = c(0.5, 0.5)
)
fit_edge <- bison_model(
  (event | effort) ~ dyad(node_1_id, node_2_id),
  data = sim_data$df_sim,
  model_type = "binary_conjugate",
  priors = get_default_priors("binary_conjugate")
)
fit_mixture <- bison_mixture(fit_edge, num_components = 3, verbose = FALSE)
summary(fit_mixture)
#> === Fitted dyadic mixture model ===
#> Maximum number of components: 3
#> Most probable number of components: 1
#> Posterior probability: 93.7%
#> === Component probabilities ===
#>            1     2     3
#> P(K=k) 0.937 0.056 0.007
#> === Component means (K = 1) ===
#>            1
#> median 0.439
#> === Edge component probabilities (K = 1) ===
#>         1
#> 1 <-> 2 1
#> 1 <-> 3 1
#> 1 <-> 4 1
#> 1 <-> 5 1
#> 1 <-> 6 1
#> 1 <-> 7 1
#> ... (99 more rows; see get_edge_component_probabilities())
```
