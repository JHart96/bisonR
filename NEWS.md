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

## Breaking changes

* cmdstanr is no longer attached by `library(bisonR)`.
* brms and mice have moved to Suggests and are only needed for `bison_brm()`.
* For undirected networks, `node_1` in `get_edgelist()` is now the node that comes first in the node levels.
* `extract_metric(obj, "edge_weight")` returns edge weights on the natural scale by default. Use `transform = FALSE` for the link scale. `bison_brm()` still uses the link scale.
* Conjugate models now give an error if the formula includes fixed effects, random effects or an intercept.
* `bison_save()` writes a single `.rds` file.
* `simulate_bison_model()` count data have more variation between edges, and a smaller age-difference effect.

## Model fitting

* Edge weights and random effects use a non-centred parameterisation, which reduces divergent transitions with `partial_pooling = TRUE`.
* `bison_model()` gains a `chains` argument, and passes `...` to the Stan sampler (e.g. `seed`, `adapt_delta`). The fixed `step_size = 0.1` has been removed.
* Input data are checked before fitting: missing columns, missing values, non-integer or negative events, more events than sampling periods, and self-loops (#69).
* Priors are checked against the model type, and missing priors are filled in with the defaults (#70).
* Compiled Stan models are cached between sessions.
* Conjugate models now return predictions and log-likelihoods, so they work with `plot_predictions()`, `prior_predictive_check()` and `model_comparison()`.
* New `print()` method for fitted models.

## Plots

* `plot_predictions()` is now an S3 generic. The `"density"` check compares event rates per dyad rather than event counts, uses random posterior draws, and draws the observed data last in a separate colour (#92). New `type = "stats"` option.
* `plot_network()` shows all nodes, gains a `threshold_on` argument, rescales edge widths for count models, and adds a caption explaining the line widths (#93). Thanks to @cebus-imitator for the threshold fixes (#94, #95).
* `plot_trace()` shows each edge in a separate panel.
* `plot_metric()` uses bar charts for discrete metrics (#57).

## `bison_brm()`

* Warns when a dyadic regression doesn't include a multi-membership term.
* Gives clearer errors when the data contain nodes or dyads that aren't in the edge model (#68).
* Supports thresholded metrics such as `node_degree[0.5]`.

## Other changes

* `bison_save()` and `bison_load()` are now exported (#82).
* `bison_to_igraph()` works with `num_draws = 1` (#73).
* `convert_gbi_to_bison()` is faster and uses the column names of the GBI matrix as node names.
* `draw_edgelist_samples()` gains a `transform` argument (#60).
* `simulate_bison_model()` gains `num_nodes` and `max_obs` arguments.
* All exported functions have examples.
* Updated articles, including a new article on preparing data (#42, #74, #81, #83, #84-#93).
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
