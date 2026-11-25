# bisonR 1.0.0

bisonR 1.0.0 is the first stable release. It fixes a number of bugs, mostly affecting edge cases and less common uses of the package, adds checks on input data and priors, and updates the documentation, tests and CI. From this release, bisonR follows semantic versioning: breaking changes will only be made in major releases.

## Bug fixes

* Fixed effects with names ending in `1` (e.g. `age1`) were dropped from the model (#80). Effort columns with digits in their names, and column names containing `.`, are now also handled correctly.
* In undirected networks, data with rows for both A-B and B-A are now treated as the same dyad.
* Conjugate models now use all rows for each dyad when data are disaggregated. Previously only the last row was used.
* Node metrics from `extract_metric()` and `bison_to_igraph()` are now always returned in the order of the node levels, with node names as column names. Previously the order could differ from the one `bison_brm()` assumed when the input data weren't sorted (#63).
* The zero-inflation probability is now bounded between 0 and 1 (#77), and posterior predictions and log-likelihoods now include the zero inflation.
* `extract_metric(standardise = TRUE)` now centres within each posterior draw, as documented.
* Mixture model components are now ordered by their means in each posterior draw, to avoid label switching.
* `global_diameter` now uses 1/weight as the edge length, consistent with betweenness and closeness.
* `global_clustering[t]` is 0 when there are no connected triples, consistent with `node_clustering[t]`.
* `prior_check()` now shows the correct density for binary edge priors on the probability scale.
* `node_degree[t]` now counts edges with weight at or above `t`, as documented.
* `model_comparison()` now works when models have different numbers of posterior draws, and accounts for autocorrelation in MCMC draws.

## Breaking changes

* bisonR now requires R 4.1.0 or later.
* cmdstanr is no longer attached by `library(bisonR)`.
* brms and mice have moved to Suggests and are only needed for `bison_brm()`, `bison_brm_get_prior()` and `bison_mice()`.
* For undirected networks, `node_1` in `get_edgelist()` is now the node that comes first in the node levels.
* Numeric node IDs are written in full (never in scientific notation), so integer and decimal columns with the same IDs give the same nodes. When node columns aren't factors with the same levels, nodes are ordered numerically if all names are numbers, and otherwise alphabetically in a locale-independent order (upper case before lower case), so results are the same on every computer. If only one node column is a factor, its unused levels are now kept.
* `extract_metric(obj, "edge_weight")` returns edge weights on the natural scale by default. Use `transform = FALSE` for the link scale. `bison_brm()` still uses the link scale.
* Conjugate models now give an error if the formula includes fixed effects, random effects or an intercept.
* Random effects are now centred on 0, and the `random_mean` prior has been removed (prior lists that include it give a warning). Previously each random effect had an estimated mean, which overlapped with the edge weights so that their overall level depended on the priors. Models with random effects but no `dyad()` term or factor fixed effect need an explicit intercept, e.g. `(event | effort) ~ 1 + (1 | group)`, and give a warning without one.
* `bison_save()` writes a single `.rds` file, and its `filename` argument is now `file`. `bison_load()` can still load models saved by earlier versions.
* `convert_gbi_to_bison()` uses the column names of the GBI matrix as node names, if it has them, rather than column numbers.
* `get_edge_component_probabilities()` columns are named `P(component 1)`, `P(component 2)`, ..., and `get_component_means()` rows `Component 1`, `Component 2`, ....
* `simulate_bison_model()` count data have more variation between edges, and a smaller age-difference effect.
* Simulated data now have an `effort` column for observation effort, and the documentation uses `(event | effort)`. The `duration` column is kept as a copy of `effort`, so existing code still works (#64).

## Model fitting

* Edge weights and random effects use a non-centred parameterisation, which reduces divergent transitions with `partial_pooling = TRUE`.
* `bison_model()` gains a `chains` argument, and passes `...` to the Stan sampler (e.g. `seed`, `adapt_delta`). The fixed `step_size = 0.1` has been removed.
* Input data are checked before fitting: missing columns, missing values, non-integer or negative events, more events than sampling periods, and self-loops (#69).
* Priors are checked against the model type, and missing priors are filled in with the defaults (#70). Standard deviations, scales and shapes must be positive.
* Formulas can use `- 1` to exclude an intercept.
* `bison_model()` warns when a zero-inflated binary model has one sampling period per row, where zero inflation can't be separated from the edge weights.
* Compiled Stan models are cached between sessions.
* Conjugate models now return predictions and log-likelihoods, so they work with `plot_predictions()`, `prior_predictive_check()` and `model_comparison()`.
* New `print()` method for fitted models.
* `bison_model()` warns if any MCMC chains fail, and if `partial_pooling = TRUE` is used without a `dyad()` term.
* `extract_metric()` warns when nodes have no dyads in the data (their values are fixed rather than estimated), and when eigenvector centrality or closeness is calculated on a network with more than one component. Malformed thresholds give an error.
* For undirected models, `bison_model()` gives a message when the data have rows for the same pair in both orders, which are combined into one edge (#71).
* `prior_predictive_check()` gains a `num_draws` argument.

## Plots

* `plot_predictions()` is now an S3 generic. The `"density"` check compares event rates per dyad rather than event counts, uses random posterior draws, and draws the observed data last in a separate colour (#92). New `type = "stats"` option.
* `plot_network()` shows all nodes, gains a `threshold_on` argument, rescales edge widths for count models, and adds a caption explaining the line widths (#93). Thanks to @cebus-imitator for the threshold fixes (#94, #95).
* `plot_trace()` shows each edge in a separate panel.
* `plot_metric()` uses bar charts for discrete metrics (#57).
* New `plot_mixture()` to plot the components of an edge mixture model (#55). The mixture models article now compares the ICL and BIC criteria (#54).
* `bison_mixture()` gains `num_draws` (default 200) and records the posterior draws it used in `draw_ids`. `criterion` is no longer case-sensitive.

## `bison_brm()`

* Warns when a dyadic regression doesn't include a multi-membership term.
* Gives clearer errors when the data contain nodes or dyads that aren't in the edge model (#68).
* Supports thresholded metrics such as `node_degree[0.5]`.
* New `subset` argument to select rows of the regression data, e.g. `subset = sex == "F"` (#58).
* New `scale` argument (`"none"`, `"centre"` or `"z"`) to transform the network variable within each posterior draw (#72). `z_score = TRUE` is equivalent to `scale = "z"`. Global metrics are scaled across networks.
* Global regressions check that there is one row of data per edge model, and that any `bison_network` column matches the order of the edge models.
* Gives an error if a network variable isn't finite (for example, closeness of a node with no dyads).
* `bison_brm_get_prior()` gains `scale`, `subset` and `z_score` arguments.

## Other changes

* `bison_save()` and `bison_load()` are now exported (#82).
* `bison_to_igraph()` works with `num_draws = 1` (#73).
* `convert_gbi_to_bison()` is faster and uses the column names of the GBI matrix as node names.
* `draw_edgelist_samples()` gains a `transform` argument (#60).
* `simulate_bison_model()` gains `num_nodes` and `max_obs` arguments.
* All exported functions have examples.
* Updated articles, and new articles on preparing data, priors, model diagnostics, additional effects in edge models, directed networks, group-based data, count data, nodal regression, comparing networks and working with posterior networks (#42, #74, #81, #83, #84-#93).
* Expanded unit tests.
* Updated documentation website.
* The GitHub repository is JHart96/bisonR and the website is <https://bisonr.org>.
* Removed unused dependencies.

# bisonR 0.4.6
* Added node_clustering network metric.
* Added global_clustering network metric.
* Added global_diameter network metric.

# bisonR 0.4.5
* Improved GBI conversion algorithm.

# bisonR 0.4.4
* Added argument pass-through in model_comparison.
* Fixed bug with degree computation.

# bisonR 0.4.3
* Various arguments to improve user control.
* Changed network plots.
* Added unit tests for metrics.
* Fixed documentation bugs.

# bisonR 0.4.2
* Calculates component means in mixture models.

# bisonR 0.4.1
* Fixed Github actions issues.
* Fixed bug in mixture model.

# bisonR 0.4.0
* Support for Z-score bison data in bison_brms.
* Example of using bison for multiple imputation.
* Duration to binary conversion utility function.
* bisonR to igraph posterior extraction utility function.
* Bug fixes.

# bisonR 0.3.0
* Integrated regressions with brms.
* Extended unit tests.
* Improved documentation.
* Partial pooling.
* Zero inflated models.

# bisonR 0.2.0

* Conjugate models
* Prior checks
* Model comparison
* Directed networks
* Support for missing dyads
* Mixture models
* Refactored data processing
* GBI data conversion
* Global network metrics
* Improved documentation
* Non-random edge weight analysis
* Threshold degree metrics

# bisonR 0.1.0

* Basic feature set implemented (edge models, dyadic and nodal regression, basic checks, etc).
* All functions are documented.
* Unit tests configured.
* R CMD check passing.
* Getting started page complete.

# bisonR 0.0.0.9000

* Initialised package.
