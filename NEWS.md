# bisonR 0.5.0

bisonR 0.5.0 is a correctness and reliability release. Several bugs could give wrong results **without any error or warning**, so we recommend re-running analyses fitted with earlier versions if any of the situations below apply to you.

## Bugs that could silently change results

* Fixed effects whose names end in `1` (e.g. `age1`) were dropped from the model and switched on an intercept (#80). Effort columns with digits in their names (e.g. `obs2`) were read as `NA`, and column names containing `.` were rejected. The formula is now parsed as an R expression, so any column name works.
* In undirected networks, data containing both `A, B` and `B, A` rows created a second, prior-only edge for the same pair. Both orders now refer to the same dyad.
* Conjugate models (`binary_conjugate`, `count_conjugate`) only used the *last* observation of each dyad when data were disaggregated (one row per observation). They now use all observations, and give the same answer for aggregated and disaggregated data.
* Node-level metrics from `extract_metric()` and `bison_to_igraph()` were ordered by the order nodes first appeared in the data, while `bison_brm()` assumed the order of the node factor levels. Unless the data were sorted, regressions on node metrics could attach values to the wrong individuals (#63). Metrics are now always in node-level order, with node names as column names, and nodes without any edges are kept.
* Zero-inflated models: the zero-inflation probability wasn't bounded above by 1 (#77), and posterior predictions and log-likelihoods (used by `plot_predictions()` and `model_comparison()`) ignored the zero inflation.
* `extract_metric(standardise = TRUE)` subtracted the overall mean instead of centring within each posterior draw as documented.
* Mixture models could suffer from label switching between posterior draws, making component means and membership probabilities unreliable. Components are now ordered by their means in every draw.
* `global_diameter` treated edge weights as distances (so stronger edges counted as *longer*). It now uses 1/weight, like betweenness and closeness.

## Breaking changes

* `library(bisonR)` no longer attaches cmdstanr. Call cmdstanr functions with `cmdstanr::`.
* brms and mice are now optional (Suggests). Install them to use `bison_brm()`.
* Undirected dyads are stored with the node that comes first in the node levels as `node_1`, which can change the order of `node_1` and `node_2` in `get_edgelist()` for some datasets.
* `extract_metric(obj, "edge_weight")` now returns natural-scale edge weights by default (probabilities or rates), consistent with `get_edgelist()`. Use `transform = FALSE` for the link scale. `bison_brm()` still uses link-scale edge weights.
* Conjugate models reject formulas with fixed or random effects or intercepts, which were previously ignored silently.
* `bison_save()` now writes a single `.rds` file (`bison_save(fit, "model.rds")`).
* Count data from `simulate_bison_model()` now have realistic variation between edges, and the age-difference effect is smaller (0.2 per year on the link scale).

## Model fitting

* Edge weights and random effects use a non-centred parameterisation, which removes the divergent transitions that were common with `partial_pooling = TRUE`.
* `bison_model()` gains `chains` and `...` (passed to Stan, e.g. `seed`, `adapt_delta`). The fixed `step_size = 0.1` was removed.
* Data are checked before fitting, with clear errors for missing columns, missing values, non-integer or negative events, more events than sampling periods, and self-loops (#69).
* Priors are checked against the model type, and missing priors are filled in from the defaults (#70).
* Compiled Stan models are cached between R sessions, so models are only compiled once per bisonR and CmdStan version. bisonR gives a clear message if CmdStan isn't installed.
* Conjugate models now provide predictions and log-likelihoods, so `plot_predictions()`, `prior_predictive_check()` and `model_comparison()` work with them.
* New `print()` method for fitted models.

## Diagnostics and plots

* `plot_predictions()` was redesigned. The `"density"` check previously compared raw event counts per dyad, which mostly reflect sampling effort and matched almost any model. It now compares event *rates*, with a shared smoothing bandwidth, randomly chosen posterior draws, and the observed data drawn on top in a contrasting colour (#92). A new `type = "stats"` check compares the proportion of zero dyads, and the mean, SD and maximum rate, with the observed values.
* `plot_network()` draws every node, can threshold on the median or lower bound (`threshold_on`), scales count-model edges sensibly, curves directed edges, explains the two line widths in a caption (#93), and returns its layout so several networks can share node positions. Thanks to @cebus-imitator for the threshold fixes (#94, #95).
* `plot_trace()` shows each edge in its own panel, titled with the dyad name.
* `plot_metric()` shows discrete metrics (e.g. degree) as bar charts (#57).

## Regression with `bison_brm()`

* Warns if a dyadic regression has no multi-membership term (`(1 | mm(node_1, node_2))`) to account for edges sharing nodes.
* Clear errors when the regression data contain nodes or dyads that aren't in the edge model (#68).
* Thresholded metrics such as `bison(node_degree[0.5](node))` can now be used.

## Other improvements

* `bison_save()` and `bison_load()` are exported and work for all model types (#82).
* `bison_to_igraph()` works with `num_draws = 1` (#73).
* `convert_gbi_to_bison()` is faster and uses column names of the GBI matrix as node names.
* `draw_edgelist_samples()` gains a `transform` argument (#60).
* `simulate_bison_model()` gains `num_nodes` and `max_obs` arguments.
* Every exported function has runnable examples.
* Documentation: rewritten *Getting started* article (addressing the feedback in #84 to #93), a new *Preparing your data* article covering effort, aggregation, directed networks and group-by-individual data (#42, #74, #81), definitions of every network metric (#83), and corrected brms examples that include multi-membership terms.
* Tests now check results numerically (e.g. against analytic posteriors and simulated true values) instead of only checking that functions run.

## Project

* Refreshed documentation site at <https://bisonr.org>: new home page, typography and colour scheme, dark mode, grouped function reference and article menu.
* Development has moved back to <https://github.com/JHart96/bisonR>, and the documentation is at <https://bisonr.org>.
* `main` is protected: changes arrive through pull requests that must pass R CMD check on all platforms and bump the version. See the contributing guide.
* Removed unused dependencies (tidyverse, formula.tools, BH, RcppEigen and others).

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
