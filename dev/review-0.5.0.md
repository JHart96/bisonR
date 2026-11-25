# bisonR code review and 0.5.0 proposals

Reviewed November 2026 against `bison-team/bisonR@16a166a` (the current `main`) and `JHart96/bisonR` (`main@ccb4556`, `dev@56a598d`). Sources: all of `R/`, `inst/stan/`, `tests/`, `vignettes/`, the CI workflows, git history, and every issue, PR and discussion on both repos. GitHub Projects couldn't be read because the `gh` token lacks the `read:project` scope.

Each finding is tagged with how it was checked:
**[ran]** means I reproduced it by sourcing the package R files in R 4.5.3 (igraph, mclust and a few other packages aren't installed locally, so I couldn't load the whole package).
**[read]** means I found it by reading the code.

---

## 0. Summary

bisonR is small (about 1,700 lines of R and 3 Stan models) and its core idea is sound. However, several silent-wrong-answer bugs affect common use cases, the tests can't catch them, and CI hasn't run successfully on either repo for a long time. I'd suggest making 0.5.0 a **correctness and trust release** rather than a feature release:

1. Fix the silent-wrong-answer bugs in §1. Each one produces plausible-looking but wrong output.
2. Replace the predictive and network plots with honest versions (§2).
3. Get CI green and rewrite the tests to check numbers, not just "no warning" (§4).
4. Rebuild the docs around tasks people actually do (§5).
5. Reconcile the two repos so there's one canonical home with redirects (§7).
6. Defer most new features to 0.6, apart from a couple of small, high-value ones (§8).

---

## 1. Correctness bugs (highest priority)

| # | Bug | Where | Impact | Check |
|---|-----|-------|--------|-------|
| 1.1 | **The formula parser swallows any term that starts with `0` or ends with `1`.** The regex `^0\|1$` means "starts with 0 **or** ends with 1". A fixed effect called `age1` disappears *and* switches the intercept on. | `R/bison_model.R:518` | Fixed effects silently dropped (this is JHart96#80 / bison-team#19). | **[ran]** `~ dyad(a,b) + age1 + sex` gives `fixed = "sex"`, `intercept = TRUE` |
| 1.2 | **A divisor variable whose name contains a digit is treated as a numeric constant.** `str_match(name, "\\d+")` matches `obs2`, then `as.numeric("obs2")` returns `NA`. | `R/bison_model.R:451` | All divisors become `NA`, so Stan errors or the conjugate model returns NaN. | **[ran]** |
| 1.3 | **Variable names containing `.` are rejected** (`body.mass`), even though dotted names are idiomatic in R. | `R/bison_model.R:532` | Warning, then the term is dropped. | **[ran]** |
| 1.4 | **Undirected data containing both A–B and B–A rows creates a phantom dyad.** `dyad_to_idx` uses `unique(cbind(n1, n2))` without sorting the pair, so B–A gets its own dyad ID. `get_dyad_ids()` then maps every row to the first match, which leaves the second dyad with no data. | `R/bison_model.R:350-359` | Extra prior-only "edges" appear in the edgelist, plots, metrics and regressions, and the dyad count is wrong. | **[ran]** 3 rows over 2 real dyads gives `num_dyads: 3` |
| 1.5 | **Conjugate models use only the last observation per dyad when data are disaggregated.** The loop over rows *overwrites* `chain[, dyad_id]` instead of accumulating events and divisors. The chain also gets `num_rows` columns instead of `num_dyads`, so the extra columns become `qlogis(0) = -Inf`. | `R/conjugate_models.R:3-20, 36-50` | Wrong posteriors. With extra columns, `get_edgelist()`/`summary()` fail and `extract_metric("edge_weight")` reshapes the data incorrectly. | **[ran]** 3/3 events gives posterior mean 0.66 (that's Beta(2,1), from the last row only) instead of 0.80 |
| 1.6 | Conjugate predictions pass a `divisor` of length `num_rows` alongside `chain[i,]` of length `num_rows` to `rbinom(num_edges, …)`. This only works when the data are aggregated. | `R/conjugate_models.R:32, 54` | Wrong posterior or prior predictive output for disaggregated data. | [read] |
| 1.7 | **Zero-inflation probability has no upper bound.** `vector<lower=0>` should be `<lower=0, upper=1>`. | `inst/stan/binary.stan:37`, `count.stan:37` | Causes the "Probability parameter is 1.015" rejections (JHart96#77 / bison-team#18). Also, with `priors_only=TRUE` the `zero_prob` prior is skipped, leaving it improper. | [read] |
| 1.8 | **`event_pred` and `log_lik` ignore zero inflation.** | `binary.stan:95-101`, `count.stan:95-101` | Posterior predictive plots and LOO comparisons for ZI models describe a *different model* from the one that was fitted. | [read] |
| 1.9 | **igraph vertex order isn't tied to `node_to_idx`.** `graph_from_edgelist()` orders vertices by first appearance in the edgelist. `extract_metric()` returns columns in that order, but `bison_mice()` indexes them by `node_to_idx` position. Nodes that appear in no dyad are dropped altogether, which makes the `num_targets` dimensions mismatch. | `R/metrics.R:23`, `R/utils.R:68`, `R/bison_mice.R` | **Node metrics get assigned to the wrong individuals in regressions** whenever the data aren't sorted by node. This is JHart96#63 / bison-team#7, still open. | [read]. It works by coincidence on the simulated data because they're sorted. |
| 1.10 | `standardise=TRUE` subtracts the *grand* mean, but the docs say "mean-center within sample". | `R/metrics.R:16, 50` | Wrong centring, so between-draw variation is left in. | [read] |
| 1.11 | **`edge_weight` is on the link (logit or log) scale in `extract_metric()` and `bison_brm()`, while every other metric is on the natural scale.** There's no option to choose and it isn't documented. | `R/metrics.R:11` | Users regress on logit edge weights without knowing it. Coefficients get misread (JHart96#60 asks for a transform argument). | [read] |
| 1.12 | `plot_trace()` guards against conjugate models with `model_type %in% c("binary_conjugate", …)`, but `model_type` has already had `_conjugate` stripped, so the guard never fires. It then crashes on `NULL$draws`. | `R/prediction.R:14` | Unhelpful error. | **[ran]** |
| 1.13 | `bison_mixture()` hard-codes `num_samples <- 200`, which overrides the actual draw count. | `R/bison_mixture.R:16` | Errors if fewer than 200 draws, otherwise ignores the rest. | [read] |
| 1.14 | Mixture component labels aren't aligned across draws. `Mclust()` doesn't guarantee component order, and labels aren't sorted by mean before averaging. | `R/bison_mixture.R:36-50` | There's a risk of label switching, which would make the component means and per-edge membership probabilities meaningless. Sort the components by mean on every draw. | [read] |
| 1.15 | `print.summary.bison_mixture` prints `x$num_components`, a field that doesn't exist, so the line comes out blank. It also rounds component means *before* `plogis`/`exp`. | `R/bison_mixture.R:110, 128` | Cosmetic, but wrong numbers. | [read] |
| 1.16 | `bison_save()`/`bison_load()` are missing `@export`, so they aren't in NAMESPACE. Their roxygen examples have a syntax error (`"testsave)`, `#' #' @examples`), and they fail on conjugate models (`fit = NULL`). | `R/saveload_bisonobj.R` | The 0.4.6 feature can't be used. | [read] |
| 1.17 | `summary()` labels directed dyads with `<->`. | `R/bison_model.R:149` | Misleading. | [read] |
| 1.18 | `extract_metric("edge_weight")` samples draws without replacement and errors if `num_draws` exceeds the draws available. `draw_edgelist_samples()` samples *with* replacement. | `R/metrics.R:13` | Inconsistent. | [read] |

The `count` simulation also uses `exp(0.01 * predictor)` (`R/simulations.R:47, 101`), so simulated count edge weights barely vary and `df_true$edge_weight` isn't on the scale of the fitted model. This matters for docs and tests (see §4).

---

## 2. Misleading outputs

### 2.1 `plot_predictions(type = "density")`

The posterior predictive plot is the one you flagged. I think it's misleading for structural reasons, not just cosmetic ones.

- **It mostly shows sampling effort, not model fit.** It plots the density of *summed events per dyad*. For binary data that sum is dominated by how many times each dyad was observed, which the model takes as given. Almost any model, including a badly wrong one, reproduces this plot, so it nearly always "passes".
- **It applies a KDE to integer counts,** with a separate bandwidth for each curve. The smoothing therefore differs between the observed and predicted curves, and small-count data (most binary datasets) look smooth and well matched when they aren't.
- **It uses the first `num_draws` rows of `event_preds`** (consecutive draws from chain 1) instead of a random sample of draws.
- For zero-inflated models the predictions come from the non-ZI likelihood (1.8).
- Issue JHart96#92 asks for the observed data to be drawn in a contrasting colour on top of the predictions.

**Proposal:**
- Plot *rates* (`event / divisor`) per dyad, using bars or a rootogram for discrete data.
- Add targeted `ppc_stat`-style checks: proportion of zero dyads, max, SD and CV of rates. Those are what actually expose over- or under-dispersion and zero inflation.
- Sample draws at random.
- Optionally hand off to `bayesplot` (`ppc_bars`, `ppc_rootogram`, `ppc_stat`), which already solves these problems.
- Make the docs say plainly what each check can and can't detect.

### 2.2 `plot_predictions(type = "point")`

This plot is fine as a sanity check, but it's labelled as a *predictive* check, which it isn't. Rename it (for example `plot_shrinkage()`) or document it as "BISoN vs naive estimates".

### 2.3 `plot_network()`

- The two repos have **three different versions** of this function. The fixes from cebus-imitator (JHart96#94 merged, JHart96#95 and bison-team#34 open) conflict with your own fix in bison-team#25.
- Thresholding on the lower bound versus the median is a real design choice. Expose it as `threshold_on = c("median", "lower")` instead of picking one silently.
- Nodes with no remaining edges vanish from the plot. Build the graph from `node_to_idx` so all nodes always appear.
- Edge width is the raw weight × `lwd`. For count models (rates above 1) this blows up, so scale widths by the maximum upper bound.
- The meaning of the overlapping lines isn't explained anywhere (JHart96#93). Add a legend or caption saying that the dark line is the lower bound and the translucent line is the upper bound.

### 2.4 `plot_metric()`

`plot_metric()` applies a logspline density to discrete metrics such as degree and diameter (JHart96#57). Use a histogram or bar chart for integer-valued metrics.

---

## 3. Statistical and methodological concerns

- **Dyadic regressions in the docs leave out node-level non-independence.** The old `dyadic_regression()` (removed in `32802db`) had multi-membership effects **on by default** (`mm=TRUE`). When it was swapped for brms, that safeguard was lost. `getting_started.Rmd` and `bison_brms.Rmd` now show `bison(edge_weight(n1, n2)) ~ age_diff` with no `(1 | mm(node_1_id, node_2_id))`, which is the mistake the BISoN papers warn about. Restore this as the documented default and consider a `bison_brm()` warning when an edge-level model has no `mm()` term.
- **`bison_brm()` uses `brm_multiple()`.** That fits one Stan model per draw, so `num_draws=100` means 100 fits, and the combined R-hat is inflated by design. JHart96#59 is "expected behaviour", but it needs documenting so users don't chase it. A faster and more principled option is the multivariate-normal (MVN) approximation of the edge posterior used as measurement error (see §8). Issues JHart96#21 and #67 point in that direction.
- **Fixed-effect factors are dummy-coded with one column per level**, plus an optional intercept. That model is only identified through the priors. Use treatment contrasts by default, or document the behaviour. A cleaner long-term fix is to build design matrices with `model.matrix()` and `lme4`-style random-effect parsing, rather than custom regex (this would also fix 1.1–1.3).
- **Partial pooling uses a centred parameterisation** (`edge_weight ~ normal(mu, edge_sigma)`), which is prone to funnel geometry. Switch to non-centred. The hard-coded `step_size=0.1` and `chains=4` should be arguments, and `mc_cores` currently also scales the number of conjugate samples (`iter_sampling * mc_cores`), which is surprising.
- **Conjugate and Stan models use different prior families** (Beta/Gamma versus Normal on the link scale), so the same `"binary"` analysis gives different answers depending on the backend. Document this and offer matched defaults.
- **The defaults are weakly informative but questionable for association data.** `normal(0, 2.5)` on the logit scale puts a lot of mass near 0 and 1, which `prior_check()` shows clearly. Explain the reasoning behind the defaults, or move to `normal(0, 1.5)`-style defaults with a vignette on choosing priors.
- **`global_density` is the mean edge weight across *observed* dyads only.** Unobserved dyads are absent rather than estimated from the prior. That's defensible, but users should know.

---

## 4. Tests and CI

### CI

- **CI doesn't run at all on either repo.** The latest runs fail with "workflow file issue", because the actions are deprecated (`actions/checkout@v2` and the matching r-lib versions). Before that, bison-team#28 notes `formatR` missing for `tidy=TRUE` in vignettes.
- **Fix:**
  - Regenerate the workflows with `usethis::use_github_action("check-standard")` and `"test-coverage"` and `"pkgdown"`.
  - Pin and cache CmdStan with `actions/cache` keyed on the version.
  - Add `formatR` to Suggests, or drop `tidy=TRUE`.

### Tests

The test suite is almost entirely smoke tests (`expect_warning(..., regexp = NA)`). The two numerical recovery tests are effectively vacuous:

- **The interval check is broken.** `mean(true < est_ub & est > est_lb)` (`test-edge_model.R:74, 109`) checks that the *median* is above the lower bound, which is always true. It should be `true > est_lb & true < est_ub`.
- **The count test compares the wrong scales.** It compares log-scale "true" weights with natural-scale estimates, and the simulated weights hardly vary anyway (§1, last paragraph).

**Proposed test plan for 0.5.0:**

1. **Unit tests for the data layer** that don't need Stan, so they're fast:
   - formula parsing, covering every term type, names with digits and dots, `0`/`1` intercepts and the `age1` regression
   - divisor as a constant versus a column
   - undirected A–B/B–A canonicalisation
   - directed versus undirected `dyad_to_idx`
   - unobserved factor levels
   - `get_dyad_ids()` against a vectorised reference
   - GBI conversion edge cases (empty groups, a single individual, constraints)
2. **Conjugate correctness:** the analytic posterior mean and variance for aggregated data should equal those for the same data disaggregated. This catches 1.5.
3. **Stan versus conjugate agreement:** on the same data with matched priors, the posterior medians should agree within tolerance.
4. **Parameter recovery with correct coverage:** binary, count, directed, with fixed effects, with random effects and with ZI, each at a fixed seed with `iter_sampling` kept small. Use `skip_on_cran()` and tag these tests as slow.
5. **Node alignment:** shuffle the input rows, then check that `extract_metric("node_strength")` columns match the strength computed by hand for each *named* node. This catches 1.9.
6. **Metric correctness:** use a small hand-built network with known strength, betweenness and clustering values.
7. **Save and load round trip.**
8. **Snapshot tests for `print`/`summary`** (`expect_snapshot()`), and `vdiffr` for plots.
9. **Compile Stan models once per test session** (helper file), with `skip_if_not(cmdstan available)`.

---

## 5. Documentation

What's there today:

- 4 vignettes, written in tutorial-paragraph style.
- **None of the exported functions has an `@examples` block.** The only two examples are the broken, unexported save and load functions.
- The vignettes contain errors:
  - They reference `draw_network_metric_samples()` and `draw_node_metric_samples()`, which don't exist.
  - They link to a `BISoN introduction` vignette that doesn't exist.
  - The brms table uses `bison(strength(node_id))`, which won't parse; it should be `node_strength`.
  - They link `summary()` to the cmdstanr docs.
- The committed `docs/` folder on `main` is a stale 0.1.0 site (the live site deploys from `gh-pages`), and it still advertises `dyadic_regression`/`nodal_regression`.
- 14 hackathon doc issues (JHart96#81 and #83–93, plus bison-team#8) are open and unaddressed.

**Proposed structure for bisonr.org (pkgdown):**

- **Get started:**
  - what BISoN is and why: "the edge model *estimates* edges; downstream models *use* them" (#88, #90)
  - installation, with CmdStan troubleshooting merged with the Windows notes
  - a 10-minute end-to-end walkthrough
- **Preparing data:**
  - required columns and what `event` and `divisor` (observation effort, not event duration) mean (#85, #86, #89)
  - aggregated versus disaggregated data, and "aggregate to the coarsest variable" (#42)
  - directed data (#81)
  - node ID factors (#42)
  - GBI and group data (#74)
  - focal follows, scans, count data and duration data
- **Edge models:**
  - binary versus count versus conjugate, and when to use each
  - choosing and checking priors
  - diagnostics: trace, R-hat and ESS, and posterior predictive checks with an honest explanation of what each detects
  - fixed and random effects in the edge model, and how they change what the edge weight means
  - partial pooling and zero inflation
- **Analyses (one article per question):**
  - Is my network non-random?
  - Does a dyadic trait predict edge weight? (with `mm()`)
  - Does a node trait predict centrality?
  - Comparing two or more networks (global metrics)
  - Edge mixture models and social "types"
  - Using posterior networks with other tools (igraph, `asnipe`, `sna`)
- **Reference:**
  - network metric definitions, stating that they're weighted versions, with thresholded versions for unweighted (#83, bison-team#8)
  - a table of the scale each output is on
  - FAQ
- **Every exported function** should get `@examples`, wrapped in `\donttest{}` where Stan is needed, plus `@seealso` links and a grouped `reference:` index in `_pkgdown.yml`.

I'd suggest using one small real dataset shipped in `data/` (with permission, for example a published association dataset), rather than simulations alone. Real data make the docs far more convincing.

---

## 6. Package engineering and best practice

- **`require()` at file top level** (`bison_model.R:1-3`, `bison_brm.R:1`) attaches packages as a side effect at build time. Remove these lines and use `pkg::fn` with `Imports`.
- **DESCRIPTION:**
  - These are in `Imports` but never used in `R/`: `tidyverse`, `knitr`, `rmarkdown`, `testthat`, `covr`, `BH`, `RcppEigen`, `formula.tools`. Move the dev tools to `Suggests` and drop the rest. `tidyverse` alone pulls in about 80 packages.
  - `brms` and `mice` are heavy and only needed for `bison_brm()`. Move them to `Suggests` with `rlang::check_installed()`.
  - `cmdstanr` in `Depends` also attaches it for users. Move it to `Imports`, or to `Suggests` if a CRAN release is ever wanted.
  - `RoxygenNote` is 7.1.2. Regenerate the docs; NAMESPACE is currently out of sync (the save and load functions are missing).
- **Model compilation:** `build_stan_model()` compiles into `tempdir()`, so every R session recompiles (about 1 minute). Cache executables in `tools::R_user_dir("bisonR", "cache")`, keyed on the Stan file hash and CmdStan version, or use the `instantiate` package to compile at install time.
- **Input validation** (JHart96#68–71): check that `model_type` is valid (use `match.arg`), that binary events are integers with `event <= divisor`, that count divisors are positive, that no row is a self-loop, that priors match the model type, and that regression data rows match the edge model's nodes. Use `cli::cli_abort()` with actionable messages.
- **Performance:** `get_dyad_ids()` is an O(rows × dyads) R loop, and `convert_gbi_to_bison()` grows a data frame with `bind_rows` inside a double loop, which is O(n²) copies. Both can be vectorised (`match()` on a canonical `"i_j"` key, and `combn` or a matrix product for the GBI). Large field datasets will notice the difference.
- **S3 hygiene:**
  - Make `plot_predictions()` a real generic (`UseMethod`).
  - Add `print.bison_model` and `as.data.frame` / `as_draws` methods.
  - Store a `bisonR_version` in fitted objects so saved models can be migrated.
  - Stop storing the `CmdStanModel` object, which points at a temp path, or make `save` and `load` handle it.
- **Duplication:** the binary and count Stan files are about 95% identical. Merge them into one file with a `family` switch, or generate both from a template.
- **Repo hygiene:**
  - Remove the committed `.Rhistory` and `.DS_Store` files (7 of them) and the stale `docs/` folder.
  - Remove `inst/stan/duration_old.stan`, and either remove the broken `duration.stan` or move it to a branch.
  - Add `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`, a PR template and `CITATION.cff` / `inst/CITATION`, which links the papers and the Zenodo DOI (bison-team#29).
- **Style:** add `lintr` and `styler` configuration and a pre-commit hook. The codebase currently mixes `=` and `<-`.

---

## 7. Reconciling `JHart96/bisonR` and `bison-team/bisonR`

### Current state

- **`bison-team/bisonR` is not a GitHub fork.** It's an independent copy with 0 stars and no releases. It's 8 commits ahead of the shared base (`f1d1b07`): issue templates, save and load, and your `plot_network` fix.
- **`JHart96/bisonR`** has the community:
  - 9 stars and 3 forks
  - 95 issue and PR numbers
  - Discussions (#23, with 8 comments)
  - the GitHub Releases and the Zenodo DOI (`v0.4.3-zenodo`)
  - all existing install instructions: `remotes::install_github("JHart96/bisonR")` is in the README, the vignettes and presumably in published methods sections
- **`JHart96/bisonR` is 2 commits ahead.** cebus-imitator's PR #94 conflicts with the bison-team fix and is superseded by it. cebus-imitator is also the only active outside contributor; they have open PRs on *both* repos.
- **Both sites are live:** `jhart96.github.io/bisonR` (rebuilt Nov 2025) and `bisonr.org` (from bison-team `gh-pages`, Apr 2025).
- **Most bison-team issues duplicate JHart96 issues.** Only bison-team #26, #28 and #29 and PR #34 are unique.

### Option A (recommended): transfer the original repo into the org

GitHub permanently redirects transferred repos for git, the web *and* the API. That means `remotes::install_github("JHart96/bisonR")` keeps working, and the stars, issues, discussions, releases and forks all come along.

1. Merge the 8 bison-team commits into `JHart96/bisonR` main. The histories share a base, so this is a normal merge. Resolve `plot_network` with the §2.3 rewrite.
2. Manually re-create bison-team #26, #28 and #29 in JHart96, and ask cebus-imitator to re-target PR #34.
3. Rename `bison-team/bisonR` to `bisonR-staging` or archive it. Transfer requires the target name to be free.
4. Transfer `JHart96/bisonR` to `bison-team`.
5. Set the Pages custom domain `bisonr.org` on the transferred repo.
6. GitHub Pages URLs don't redirect after a transfer, so `jhart96.github.io/bisonR` will 404. To keep old links alive, publish a tiny `JHart96/JHart96.github.io` user site whose `bisonR/**` pages meta-refresh to `bisonr.org`.
7. Re-enable the Zenodo GitHub integration for the new owner so future releases mint new versions under the same concept DOI.
8. **Never** create a new repo called `JHart96/bisonR` afterwards, because doing so silently breaks the redirect.

### Option B: keep bison-team as canonical and archive the original

1. Merge the JHart96-only commits into bison-team.
2. Push a final `0.4.7` to JHart96 whose only change is an `.onAttach` message: "bisonR has moved to github.com/bison-team/bisonR, reinstall with `remotes::install_github('bison-team/bisonR')`". Add a README banner and pin an announcement issue and discussion.
3. Replace JHart96 `gh-pages` with redirect stubs to bisonr.org.
4. Archive the repo (read-only).

The cost is that stars, issue history and Discussions stay behind, and anyone with the old URL in a script never gets 0.5.0 unless they see the message.

### For either option

- Update every `JHart96` reference: README badges, install lines, `_pkgdown.yml` `url:` (still `jhart96.github.io`), vignette install chunks and `DESCRIPTION` `BugReports:` (missing; add it).
- Announce the move in the 0.5.0 NEWS.

---

## 8. Lost, stalled and roadmap features

Collected from git history, closed-but-unfinished issues and code that's still in the tree.

| Feature | Evidence | Status | Suggestion |
|---|---|---|---|
| **Duration model** (exponential durations + Poisson event counts) | Issue #41 closed; "Removed duration model for now" (`334ed56`, Sep 2022); `inst/stan/duration.stan` still shipped but unreachable and broken (it references `sender_receiver`, indexes `dyad_ids[i]` over `num_edges`, and uses an unbounded `rate ~ normal`) | Lost | High value for focal-follow and grooming data. Target 0.6 with a proper design. Remove the dead file in 0.5. |
| **Sender/receiver (correlated directed) effects** | Issue #44 closed; "Started sender/receiver effects" (`728dd95`); half-written in `duration.stan` | Lost | Natural for directed networks: an MVN prior on (A→B, B→A) plus node sender and receiver effects. Target 0.6. |
| **Multi-membership by default in dyadic regression** | Old `dyadic_regression(mm=TRUE)`, removed in `32802db` | Lost in the brms migration | **Restore for 0.5** through the docs and a `bison_brm()` warning (§3). |
| **MVN posterior approximation** for downstream regression | Issue #21 closed (accuracy check); issue #67 "Replace imputation with se()" closed unimplemented; used in the BISoN example repos | Stalled | A single brms fit with `me()`/`mi()` or `fcor`-style MVN error. That's 100× faster than `brm_multiple` with valid diagnostics. Strong 0.5 or 0.6 candidate. |
| **INLA backend** | "Added INLA to edge_model function" (`7b6e185`, May 2022) | Dropped | Probably not worth reviving. The conjugate and Stan paths cover it. |
| **Cox / survival model** | Issue #34 closed with no code | Never started | Defer, or fold into the duration model. |
| `transform` argument for draws | JHart96#60 / bison-team#5 | Open | **0.5**. Fixes the 1.11 inconsistency. |
| `scale` option for `bison_brm` | #72 / #15 | Open | 0.5, and generalise `z_score`. |
| Row filtering in `bison_brm` | #58 / #3 | Open | 0.5, small. |
| Mixture component plots and ICL comparison | #55, #54 / #1 | Open | 0.5 together with the label-switching fix. |
| `plot_network` filter option | #66 / #10 | Open | Covered by the `plot_network` rewrite. |
| Auto-detect directedness | #71 / #14 | Open | Warn only; never switch silently. |
| Simulation column names | #64 / #9 | Open | 0.5. Use `effort` instead of `duration` to remove the event-duration confusion (#89). |
| Save and load | #82 / #21 | Half-done | Finish and export (1.16). |

Ideas that came out of the review rather than from existing issues:

- `bison_model()` support for **multiple networks in one fit**, for example per-season networks with shared priors.
- A `tidy()`/`broom` interface and a `posterior::as_draws()` method, so the ecosystem tools (bayesplot, tidybayes) just work.
- A **`bison_data()` helper** that validates and standardises input. Most user pain is at the data-preparation step.
- Precompiled models or `instantiate`, so first use doesn't involve a long compile.

---

## 9. Proposed 0.5.0 milestone

**Must (correctness and infrastructure)**
- CI green on all three operating systems; CmdStan cached; coverage reported
- Fix 1.1–1.12 and 1.16, each with a regression test from §4
- Zero-inflation bounds and ZI-aware `event_pred`/`log_lik`
- Node-ordering fix (1.9), with a shuffled-row test
- `transform=` argument and a documented scale for all extracted quantities
- Clean up DESCRIPTION: dependencies, `require()`, roxygen regenerated
- Resolve the three `plot_network` variants; merge or close the cebus-imitator PRs with thanks
- Repo reconciliation (§7), including install-URL updates and a deprecation message

**Should (trust and usability)**
- Rewrite the predictive checks (§2.1) and the network plot (§2.3)
- Restore multi-membership guidance and the `bison_brm()` warning
- Mixture label alignment, component plot, and using the real draw count
- Input validation with clear errors (#68–71)
- `@examples` on every exported function; restructure the docs (§5), closing the 14 hackathon issues
- Cache compiled models across sessions
- Vectorise `get_dyad_ids()` and the GBI conversion

**Could (pull in if time allows, otherwise 0.6)**
- MVN measurement-error regression path as a faster alternative to `brm_multiple`
- Non-centred partial pooling; expose `chains`/`step_size`/`adapt_delta`
- Merge the binary and count Stan models
- `as_draws()`/`tidy()` methods

**0.6+**
- Duration model, sender/receiver effects, multi-network fits

### Breaking changes to flag in NEWS

- The formula parser is stricter: names containing dots are now accepted, but names that used to be silently dropped will now count.
- Undirected dyads are canonicalised, so dyad counts and IDs can change for some datasets.
- The default scale of `edge_weight` in `extract_metric()` changes if you choose natural-scale defaults.
- Conjugate results change for disaggregated data. They're now correct.
- `cmdstanr` is no longer attached on `library(bisonR)`.

### A note on versioning

bison-team#29 proposes CalVer (`YYYY.M.COMMIT`). Since you're planning `0.5.0`, I'd stay on SemVer until 1.0. CalVer would jump from `0.4.6` to `2026.x`, which works for ordering but signals nothing about breaking changes, and the release above has several.
