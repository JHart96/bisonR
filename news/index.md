# Changelog

## bisonR 1.0.0

bisonR 1.0.0 is the first stable release. It fixes a number of bugs,
mostly affecting edge cases and less common uses of the package, adds
checks on input data and priors, and updates the documentation, tests
and CI. From this release, bisonR follows semantic versioning: breaking
changes will only be made in major releases.

### Bug fixes

- Fixed effects with names ending in `1` (e.g. `age1`) were dropped from
  the model ([\#80](https://github.com/JHart96/bisonR/issues/80)).
  Effort columns with digits in their names, and column names containing
  `.`, are now also handled correctly.
- In undirected networks, data with rows for both A-B and B-A are now
  treated as the same dyad.
- Conjugate models now use all rows for each dyad when data are
  disaggregated. Previously only the last row was used.
- Node metrics from
  [`extract_metric()`](https://bisonr.org/reference/extract_metric.md)
  and
  [`bison_to_igraph()`](https://bisonr.org/reference/bison_to_igraph.md)
  are now always returned in the order of the node levels, with node
  names as column names. Previously the order could differ from the one
  [`bison_brm()`](https://bisonr.org/reference/bison_brm.md) assumed
  when the input data weren’t sorted
  ([\#63](https://github.com/JHart96/bisonR/issues/63)).
- The zero-inflation probability is now bounded between 0 and 1
  ([\#77](https://github.com/JHart96/bisonR/issues/77)), and posterior
  predictions and log-likelihoods now include the zero inflation.
- `extract_metric(standardise = TRUE)` now centres within each posterior
  draw, as documented.
- Mixture model components are now ordered by their means in each
  posterior draw, to avoid label switching.
- `global_diameter` now uses 1/weight as the edge length, consistent
  with betweenness and closeness.
- `global_clustering[t]` is 0 when there are no connected triples,
  consistent with `node_clustering[t]`.
- [`simulate_bison_model()`](https://bisonr.org/reference/simulate_bison_model.md)
  put the true `age_diff`, `age_1` and `age_2` values in the wrong
  columns of `df_true`. The simulated data in `df_sim` and the true edge
  weights were not affected.
- [`simulate_bison_model_mixture()`](https://bisonr.org/reference/simulate_bison_model_mixture.md)
  no longer returns `age_diff`, `age_1`, `age_2` and `location` columns,
  which were never simulated and contained other values.
- [`prior_check()`](https://bisonr.org/reference/prior_check.md) now
  shows the correct density for binary edge priors on the probability
  scale.
- `node_degree[t]` now counts edges with weight at or above `t`, as
  documented.
- [`model_comparison()`](https://bisonr.org/reference/model_comparison.md)
  now works when models have different numbers of posterior draws, and
  accounts for autocorrelation in MCMC draws.
- [`get_edge_component_probabilities()`](https://bisonr.org/reference/get_edge_component_probabilities.md)
  and
  [`get_component_means()`](https://bisonr.org/reference/get_component_means.md)
  check `num_components`,
  [`simulate_bison_model()`](https://bisonr.org/reference/simulate_bison_model.md)
  checks `model_type`, and
  [`bison_brm()`](https://bisonr.org/reference/bison_brm.md) accepts a
  `backend` argument.

### Breaking changes

- bisonR now requires R 4.1.0 or later.
- cmdstanr is no longer attached by
  [`library(bisonR)`](https://bisonr.org).
- brms and mice have moved to Suggests and are only needed for
  [`bison_brm()`](https://bisonr.org/reference/bison_brm.md),
  [`bison_brm_get_prior()`](https://bisonr.org/reference/bison_brm_get_prior.md)
  and [`bison_mice()`](https://bisonr.org/reference/bison_mice.md).
- For undirected networks, `node_1` in
  [`get_edgelist()`](https://bisonr.org/reference/get_edgelist.md) is
  now the node that comes first in the node levels.
- Numeric node IDs are written in full (never in scientific notation),
  so integer and decimal columns with the same IDs give the same nodes.
  When node columns aren’t factors with the same levels, nodes are
  ordered numerically if all names are numbers, and otherwise
  alphabetically in a locale-independent order (upper case before lower
  case), so results are the same on every computer. If only one node
  column is a factor, its unused levels are now kept.
- `extract_metric(obj, "edge_weight")` returns edge weights on the
  natural scale by default. Use `transform = FALSE` for the link scale.
  [`bison_brm()`](https://bisonr.org/reference/bison_brm.md) still uses
  the link scale.
- Conjugate models now give an error if the formula includes fixed
  effects, random effects or an intercept.
- Random effects are now centred on 0, and the `random_mean` prior has
  been removed (prior lists that include it give a warning). Previously
  each random effect had an estimated mean, which overlapped with the
  edge weights so that their overall level depended on the priors.
  Models with random effects but no `dyad()` term or factor fixed effect
  need an explicit intercept, e.g. `(event | effort) ~ 1 + (1 | group)`,
  and give a warning without one.
- [`bison_save()`](https://bisonr.org/reference/bison_save.md) writes a
  single `.rds` file, and its `filename` argument is now `file`.
  [`bison_load()`](https://bisonr.org/reference/bison_load.md) can still
  load models saved by earlier versions.
- [`convert_gbi_to_bison()`](https://bisonr.org/reference/convert_gbi_to_bison.md)
  uses the column names of the GBI matrix as node names, if it has them,
  rather than column numbers.
- [`get_edge_component_probabilities()`](https://bisonr.org/reference/get_edge_component_probabilities.md)
  columns are named `P(component 1)`, `P(component 2)`, …, and
  [`get_component_means()`](https://bisonr.org/reference/get_component_means.md)
  rows `Component 1`, `Component 2`, ….
- [`simulate_bison_model()`](https://bisonr.org/reference/simulate_bison_model.md)
  count data have more variation between edges, and a smaller
  age-difference effect.
- Simulated data now have an `effort` column for observation effort, and
  the documentation uses `(event | effort)`. The `duration` column is
  kept as a copy of `effort`, so existing code still works
  ([\#64](https://github.com/JHart96/bisonR/issues/64)).

### Model fitting

- Edge weights and random effects use a non-centred parameterisation,
  which reduces divergent transitions with `partial_pooling = TRUE`.
- [`bison_model()`](https://bisonr.org/reference/bison_model.md) gains a
  `chains` argument, and passes `...` to the Stan sampler (e.g. `seed`,
  `adapt_delta`). The fixed `step_size = 0.1` has been removed.
- Input data are checked before fitting: missing columns, missing
  values, non-integer or negative events, more events than sampling
  periods, and self-loops
  ([\#69](https://github.com/JHart96/bisonR/issues/69)).
- Priors are checked against the model type, and missing priors are
  filled in with the defaults
  ([\#70](https://github.com/JHart96/bisonR/issues/70)). Standard
  deviations, scales and shapes must be positive.
- Formulas can use `- 1` to exclude an intercept.
- [`bison_model()`](https://bisonr.org/reference/bison_model.md) warns
  when a zero-inflated binary model has one sampling period per row,
  where zero inflation can’t be separated from the edge weights.
- Compiled Stan models are cached between sessions.
- Conjugate models now return predictions and log-likelihoods, so they
  work with
  [`plot_predictions()`](https://bisonr.org/reference/plot_predictions.md),
  [`prior_predictive_check()`](https://bisonr.org/reference/prior_predictive_check.md)
  and
  [`model_comparison()`](https://bisonr.org/reference/model_comparison.md).
- New [`print()`](https://rdrr.io/r/base/print.html) method for fitted
  models.
- New [`get_effects()`](https://bisonr.org/reference/get_effects.md)
  returns the estimates of fixed effects, random effects, partial
  pooling and zero inflation, with readable labels such as `sex = F`.
  [`summary()`](https://rdrr.io/r/base/summary.html) now shows these
  after the edge weights, and works for models without a `dyad()` term.
- [`bison_model()`](https://bisonr.org/reference/bison_model.md) warns
  if any MCMC chains fail, and if `partial_pooling = TRUE` is used
  without a `dyad()` term.
- [`extract_metric()`](https://bisonr.org/reference/extract_metric.md)
  warns when nodes have no dyads in the data (their values are fixed
  rather than estimated), and when eigenvector centrality or closeness
  is calculated on a network with more than one component. Malformed
  thresholds give an error.
- For undirected models,
  [`bison_model()`](https://bisonr.org/reference/bison_model.md) gives a
  message when the data have rows for the same pair in both orders,
  which are combined into one edge
  ([\#71](https://github.com/JHart96/bisonR/issues/71)).
- [`prior_predictive_check()`](https://bisonr.org/reference/prior_predictive_check.md)
  gains a `num_draws` argument.

### Plots

- [`plot_predictions()`](https://bisonr.org/reference/plot_predictions.md)
  is now an S3 generic. The `"density"` check compares event rates per
  dyad rather than event counts, uses random posterior draws, and draws
  the observed data last in a separate colour
  ([\#92](https://github.com/JHart96/bisonR/issues/92)). New
  `type = "stats"` option.
- [`plot_network()`](https://bisonr.org/reference/plot_network.md) draws
  light nodes with dark labels sized to fit the node names, and node
  settings passed through `...` (such as `vertex.color` or
  `vertex.label.cex`) now replace the defaults. It also shows all nodes,
  gains a `threshold_on` argument, rescales edge widths for count
  models, and adds a caption explaining the line widths
  ([\#93](https://github.com/JHart96/bisonR/issues/93)). Thanks to
  [@cebus-imitator](https://github.com/cebus-imitator) for the threshold
  fixes ([\#94](https://github.com/JHart96/bisonR/issues/94),
  [\#95](https://github.com/JHart96/bisonR/issues/95)).
- [`plot_trace()`](https://bisonr.org/reference/plot_trace.md) shows
  each edge in a separate panel.
- [`plot_metric()`](https://bisonr.org/reference/plot_metric.md) uses
  bar charts for discrete metrics
  ([\#57](https://github.com/JHart96/bisonR/issues/57)).
- New [`plot_mixture()`](https://bisonr.org/reference/plot_mixture.md)
  to plot the components of an edge mixture model
  ([\#55](https://github.com/JHart96/bisonR/issues/55)). The mixture
  models article now compares the ICL and BIC criteria
  ([\#54](https://github.com/JHart96/bisonR/issues/54)).
- [`bison_mixture()`](https://bisonr.org/reference/bison_mixture.md)
  gains `num_draws` (default 200) and records the posterior draws it
  used in `draw_ids`. `criterion` is no longer case-sensitive.

### `bison_brm()`

- Warns when a dyadic regression doesn’t include a multi-membership
  term.
- Gives clearer errors when the data contain nodes or dyads that aren’t
  in the edge model
  ([\#68](https://github.com/JHart96/bisonR/issues/68)).
- Supports thresholded metrics such as `node_degree[0.5]`.
- New `subset` argument to select rows of the regression data,
  e.g. `subset = sex == "F"`
  ([\#58](https://github.com/JHart96/bisonR/issues/58)).
- New `scale` argument (`"none"`, `"centre"` or `"z"`) to transform the
  network variable within each posterior draw
  ([\#72](https://github.com/JHart96/bisonR/issues/72)).
  `z_score = TRUE` is equivalent to `scale = "z"`. Global metrics are
  scaled across networks.
- Global regressions check that there is one row of data per edge model,
  and that any `bison_network` column matches the order of the edge
  models.
- Gives an error if a network variable isn’t finite (for example,
  closeness of a node with no dyads).
- [`bison_brm_get_prior()`](https://bisonr.org/reference/bison_brm_get_prior.md)
  gains `scale`, `subset` and `z_score` arguments.

### Other changes

- [`bison_save()`](https://bisonr.org/reference/bison_save.md) and
  [`bison_load()`](https://bisonr.org/reference/bison_load.md) are now
  exported ([\#82](https://github.com/JHart96/bisonR/issues/82)).
- [`bison_to_igraph()`](https://bisonr.org/reference/bison_to_igraph.md)
  works with `num_draws = 1`
  ([\#73](https://github.com/JHart96/bisonR/issues/73)).
- [`convert_gbi_to_bison()`](https://bisonr.org/reference/convert_gbi_to_bison.md)
  is faster and uses the column names of the GBI matrix as node names.
- [`draw_edgelist_samples()`](https://bisonr.org/reference/draw_edgelist_samples.md)
  gains a `transform` argument
  ([\#60](https://github.com/JHart96/bisonR/issues/60)).
- [`simulate_bison_model()`](https://bisonr.org/reference/simulate_bison_model.md)
  gains `num_nodes` and `max_obs` arguments.
- All exported functions have examples.
- Updated articles, and new articles on preparing data, priors, model
  diagnostics, additional effects in edge models, directed networks,
  group-based data, count data, nodal regression, comparing networks and
  working with posterior networks
  ([\#42](https://github.com/JHart96/bisonR/issues/42),
  [\#74](https://github.com/JHart96/bisonR/issues/74),
  [\#81](https://github.com/JHart96/bisonR/issues/81),
  [\#83](https://github.com/JHart96/bisonR/issues/83),
  [\#84](https://github.com/JHart96/bisonR/issues/84)-#93).
- Expanded unit tests.
- Updated documentation website.
- Rewritten help pages for all functions, with more context on what each
  function is for and when to use it.
- Michael Weiss, Daniel Franks, Lauren Brent and Samuel Ellis are listed
  as authors, and Joe Wilde as a contributor.
- The GitHub repository is JHart96/bisonR and the website is
  <https://bisonr.org>.
- Removed unused dependencies.

## bisonR 0.4.6

- Added node_clustering network metric.
- Added global_clustering network metric.
- Added global_diameter network metric.

## bisonR 0.4.5

- Improved GBI conversion algorithm.

## bisonR 0.4.4

- Added argument pass-through in model_comparison.
- Fixed bug with degree computation.

## bisonR 0.4.3

- Various arguments to improve user control.
- Changed network plots.
- Added unit tests for metrics.
- Fixed documentation bugs.

## bisonR 0.4.2

- Calculates component means in mixture models.

## bisonR 0.4.1

- Fixed Github actions issues.
- Fixed bug in mixture model.

## bisonR 0.4.0

- Support for Z-score bison data in bison_brms.
- Example of using bison for multiple imputation.
- Duration to binary conversion utility function.
- bisonR to igraph posterior extraction utility function.
- Bug fixes.

## bisonR 0.3.0

- Integrated regressions with brms.
- Extended unit tests.
- Improved documentation.
- Partial pooling.
- Zero inflated models.

## bisonR 0.2.0

- Conjugate models
- Prior checks
- Model comparison
- Directed networks
- Support for missing dyads
- Mixture models
- Refactored data processing
- GBI data conversion
- Global network metrics
- Improved documentation
- Non-random edge weight analysis
- Threshold degree metrics

## bisonR 0.1.0

- Basic feature set implemented (edge models, dyadic and nodal
  regression, basic checks, etc).
- All functions are documented.
- Unit tests configured.
- R CMD check passing.
- Getting started page complete.

## bisonR 0.0.0.9000

- Initialised package.
