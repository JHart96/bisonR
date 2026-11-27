
utils::globalVariables(c(".imp"))

#' Fit a regression with a network metric
#'
#' Fits a Bayesian regression with the brms package in which a network metric,
#' such as node strength or edge weight, is a response or predictor. The
#' uncertainty in the metric, estimated by the edge model, is carried into the
#' regression by multiple imputation. Use it after [bison_model()] to test
#' hypotheses about edges, nodes or whole networks, for example whether age
#' predicts node strength.
#'
#' @details
#' The network metric is written in the formula as a `bison()` term containing
#' the metric and the data columns that identify its targets:
#' * `bison(edge_weight(node_1_id, node_2_id))` for dyadic regression;
#' * `bison(node_strength(node_id))`, `bison(node_degree[0.2](node_id))` and so
#'   on for nodal regression;
#' * `bison(global_cv(bison_network))` and so on for global metrics, with one
#'   row of data per network.
#'
#' The `bison()` term can be on either side of the formula, and the rest of the
#' formula uses ordinary brms syntax, including random effects such as
#' `(1 | mm(node_1_id, node_2_id))`. In the fitted model, the network variable is
#' named after the metric, such as `bison_node_strength`. See [extract_metric()]
#' for the metrics available.
#'
#' Multiple imputation repeats an analysis on several plausible versions of the
#' data and combines the results. Here, each version uses one posterior draw of
#' the network metric (see [bison_mice()]). [brms::brm_multiple()] fits the
#' regression to each version and combines the posterior draws of all the fits,
#' so `num_draws` regressions are fitted and fitting time grows with
#' `num_draws`. Models are fitted with the cmdstanr backend unless `backend` is
#' given in `...`.
#'
#' Edge weights are used on the link scale (log-odds for binary models, log
#' rates for count models). Node and global metrics are calculated from edge
#' weights on the natural scale.
#'
#' With a list of edge models, a factor `bison_network` is added to the data,
#' giving the position of each edge model in the list. It can be used as a
#' predictor or random effect to compare networks.
#'
#' @param formula A brms formula containing one `bison()` term (see *Details*).
#' @param edgemodel_list A fitted `bison_model` object, or a list of them for
#'   an analysis of several networks.
#' @param data_list A data frame of regression variables, or a list of data
#'   frames with one per edge model, in the same order. For edge and node
#'   metrics, each row refers to a dyad or node, and the columns named in the
#'   `bison()` term must hold node names that match the edge model. For global
#'   metrics, give one data frame with one row per edge model.
#' @param num_draws Number of posterior draws of the network metric to use,
#'   which is the number of regressions fitted. 100 is a reasonable value for
#'   a final analysis; smaller values (such as 5 or 10) are useful while
#'   developing a model.
#' @param z_score If `TRUE`, standardise the network variable within each
#'   posterior draw. Equivalent to `scale = "z"`.
#' @param cores Number of processor cores to use.
#' @param chains Number of MCMC chains to run for each regression. MCMC (Markov
#'   chain Monte Carlo) draws samples from the posterior; a chain is one
#'   independent run of the sampler.
#' @param scale How to transform the network variable within each posterior
#'   draw: `"none"` (default), `"centre"` (subtract the mean) or `"z"` (subtract
#'   the mean and divide by the standard deviation). For edge and node metrics,
#'   transformations are applied across the rows used in the model, separately
#'   for each network, which removes differences in the average level between
#'   networks; use `"none"` to compare levels between networks. Global metrics
#'   have one value per network, so they are transformed across the networks.
#' @param subset An optional logical expression, evaluated in each data frame,
#'   selecting the rows to use, e.g. `subset = sex == "F"`. Network metrics are
#'   still calculated from the whole network; only the rows of the regression
#'   are selected. Not available for global metrics.
#' @param ... Additional arguments passed to [brms::brm_multiple()], such as
#'   `family` or `prior`.
#'
#' @return A fitted brms model (a `brmsfit` object) combining the fits to all
#'   imputed datasets. Use `summary()`, [brms::conditional_effects()] or
#'   [brms::hypothesis()] on it as for any brms model.
#'
#' @seealso [bison_brm_get_prior()] to see and set the priors of the
#'   regression, [bison_mice()] for the imputed data, and [extract_metric()].
#' @examplesIf requireNamespace("brms", quietly = TRUE) && requireNamespace("mice", quietly = TRUE) && !identical(Sys.getenv("BISONR_TEST_BRMS"), "false")
#' sim_data <- simulate_bison_model("binary", aggregated = TRUE)
#' df <- sim_data$df_sim
#' fit_edge <- bison_model(
#'   (event | effort) ~ dyad(node_1_id, node_2_id),
#'   data = df,
#'   model_type = "binary_conjugate",
#'   priors = get_default_priors("binary_conjugate")
#' )
#' # Does age difference predict edge weight? The multi-membership term accounts
#' # for edges that share a node not being independent.
#' fit_dyadic <- bison_brm(
#'   bison(edge_weight(node_1_id, node_2_id)) ~ age_diff + (1 | mm(node_1_id, node_2_id)),
#'   fit_edge,
#'   df,
#'   num_draws = 5 # Use more (e.g. 100) in a real analysis
#' )
#' summary(fit_dyadic)
#' @export
bison_brm <- function(formula, edgemodel_list, data_list, num_draws=100, z_score=FALSE, cores=4, chains=2,
                      scale = c("none", "centre", "z"), subset = NULL, ...) {
  check_suggested("brms", "bison_brm")
  check_suggested("mice", "bison_brm")
  scale <- match.arg(scale)
  # Parse formula
  parsed_formula <- parse_bison_brms_formula(formula)
  brms_formula <- parsed_formula$brms_formula
  data_list <- subset_bison_data(data_list, substitute(subset), parent.frame(), parsed_formula$target_name)

  num_models <- if (inherits(edgemodel_list, "bison_model")) 1 else length(edgemodel_list)
  check_network_scaling(if (isTRUE(z_score)) "z" else scale, num_models, parsed_formula$target_name, brms_formula)

  if (parsed_formula$target_name == "edge" && !grepl("mm\\(", brms_formula)) {
    warning(
      "Dyadic regressions should normally account for the non-independence of edges that share ",
      "a node, for example by adding a multi-membership term: + (1 | mm(",
      parsed_formula$param_names[1], ", ", parsed_formula$param_names[2], ")).",
      call. = FALSE
    )
  }

  # Generate mice object
  mice_obj <- bison_mice(
    edgemodel_list,
    data_list,
    parsed_formula$param_names,
    parsed_formula$target_name,
    parsed_formula$metric_name,
    num_draws,
    z_score,
    scale = scale
  )

  # Run brms imputation. brm_multiple() warns that E-BFMI can't be computed for
  # the combined fit; that warning isn't informative here, so it is muffled.
  brm_args <- list(...)
  if (is.null(brm_args$backend)) brm_args$backend <- "cmdstanr"
  withCallingHandlers(
    do.call(brms::brm_multiple, c(list(brms_formula, mice_obj, cores = cores, chains = chains), brm_args)),
    warning = function(w) {
      if (grepl("E-BFMI not computed", conditionMessage(w), fixed = TRUE)) invokeRestart("muffleWarning")
    }
  )
}

#' List the priors of a bison_brm() regression
#'
#' Shows the parameters of a [bison_brm()] regression that can be given priors,
#' and their default priors, using [brms::get_prior()]. Use it before fitting
#' to choose priors, which can then be set with the `prior` argument of
#' [bison_brm()]. A prior is the distribution that describes what is plausible
#' for a parameter before seeing the data.
#'
#' @param formula A brms formula containing one `bison()` term, as for
#'   [bison_brm()].
#' @param edgemodel_list A fitted `bison_model` object, or a list of them.
#' @param data_list A data frame of regression variables, or a list of them,
#'   as for [bison_brm()].
#' @param scale,subset,z_score As for [bison_brm()]. The scale of the network
#'   variable affects the default priors, so use the same values as in
#'   [bison_brm()].
#' @param ... Additional arguments passed to [brms::get_prior()], such as
#'   `family`.
#'
#' @return A data frame of priors (a `brmsprior` object) from
#'   [brms::get_prior()], with one row per parameter or class of parameters.
#'
#' @seealso [bison_brm()], [brms::set_prior()]
#'
#' @examplesIf requireNamespace("brms", quietly = TRUE) && requireNamespace("mice", quietly = TRUE) && !identical(Sys.getenv("BISONR_TEST_BRMS"), "false")
#' sim_data <- simulate_bison_model("binary", aggregated = TRUE)
#' df <- sim_data$df_sim
#' fit_edge <- bison_model(
#'   (event | effort) ~ dyad(node_1_id, node_2_id),
#'   data = df,
#'   model_type = "binary_conjugate",
#'   priors = get_default_priors("binary_conjugate")
#' )
#' bison_brm_get_prior(
#'   bison(edge_weight(node_1_id, node_2_id)) ~ age_diff + (1 | mm(node_1_id, node_2_id)),
#'   fit_edge,
#'   df
#' )
#' @export
bison_brm_get_prior <- function(formula, edgemodel_list, data_list, scale = c("none", "centre", "z"),
                                subset = NULL, z_score = FALSE, ...) {
  scale <- match.arg(scale)
  if (isTRUE(z_score)) {
    if (!scale %in% c("none", "z")) {
      warning("`z_score = TRUE` overrides `scale = \"", scale, "\"`; using scale = \"z\".", call. = FALSE)
    }
    scale <- "z"
  }
  check_suggested("brms", "bison_brm_get_prior")
  check_suggested("mice", "bison_brm_get_prior")
  # Parse formula
  parsed_formula <- parse_bison_brms_formula(formula)
  brms_formula <- parsed_formula$brms_formula
  data_list <- subset_bison_data(data_list, substitute(subset), parent.frame(), parsed_formula$target_name)

  # Generate mice object
  mice_obj <- bison_mice(
    edgemodel_list,
    data_list,
    parsed_formula$param_names,
    parsed_formula$target_name,
    parsed_formula$metric_name,
    num_draws=2,
    scale = scale
  )

  data <- mice::complete(mice_obj)

  brms::get_prior(brms_formula, data, ...)
}

parse_bison_brms_formula <- function(formula) {
  # Parse formula string
  formula_string <- deparse1(formula)[1]

  # Locate bison(...) term
  locs <- stringr::str_locate(formula_string, "bison\\([a-zA-Z0-9_.,\\[\\]\\(\\)\\s]*\\)")
  if (is.na(locs[1])) {
    stop("The formula must contain a bison(...) term, e.g. bison(node_strength(node_id)).", call. = FALSE)
  }
  formula_string_pre <- substr(formula_string, 1, locs[1] - 1)
  formula_string_post <- substr(formula_string, locs[2] + 1, nchar(formula_string))

  # Extract bison(...) term
  bison_term <- stringr::str_extract(formula_string, "bison\\([a-zA-Z0-9_.,\\[\\]\\(\\)\\s]*\\)")
  bison_term <- stringr::str_replace_all(bison_term, "\\s", "")

  # Extract function name and arguments
  bison_term_split <- stringr::str_split(bison_term, "bison\\(|\\(|\\)")[[1]]
  bison_term_split <- bison_term_split[nchar(bison_term_split) > 0]
  metric_name <- bison_term_split[1]
  params <- bison_term_split[2]

  # Extract parameter names
  param_names <- stringr::str_split(params, ",")[[1]]
  metric_parts <- stringr::str_split(metric_name, "_", n = 2)[[1]]
  if (length(metric_parts) != 2 || !metric_parts[1] %in% c("edge", "node", "global")) {
    stop(
      "Network metrics in bison() must start with edge_, node_ or global_, e.g. ",
      "bison(node_strength(node_id)), but got \"", metric_name, "\".",
      call. = FALSE
    )
  }
  target_name <- metric_parts[1]
  metric_name <- metric_parts[2]

  # Generate bison term for brms model.
  new_bison_term <- bison_term_name(target_name, metric_name)
  brms_formula <- paste0(formula_string_pre, new_bison_term, formula_string_post)
  if (inherits(tryCatch(stats::as.formula(brms_formula), error = function(e) e), "error")) {
    stop(
      "Couldn't read the bison() term in the formula. Write it without surrounding brackets, ",
      "e.g. y ~ bison(node_strength(node)) + x.",
      call. = FALSE
    )
  }

  return(list(
    param_names=param_names,
    target_name=target_name,
    metric_name=metric_name,
    new_bison_term=new_bison_term,
    brms_formula=brms_formula
  ))
}

# Name of the data column holding a network metric in bison_brm() models, e.g.
# "bison_node_strength" or "bison_node_degree_0_5" for node_degree[0.5].
bison_term_name <- function(target_name, metric_name) {
  name <- paste0("bison_", target_name, "_", metric_name)
  name <- gsub("[^A-Za-z0-9_]+", "_", name)
  sub("_+$", "", name)
}

# Select rows of the regression data with a logical expression evaluated in
# each data frame.
subset_bison_data <- function(data_list, subset_expr, env, target_name) {
  if (is.null(subset_expr)) {
    return(data_list)
  }
  if (target_name == "global") {
    stop("`subset` can't be used with global metrics; select the edge models and rows directly.", call. = FALSE)
  }
  select_rows <- function(df) {
    keep <- eval(subset_expr, df, env)
    if (!is.logical(keep) || length(keep) != nrow(df)) {
      stop("`subset` must be a logical expression with one value per row of the data.", call. = FALSE)
    }
    keep[is.na(keep)] <- FALSE
    if (!any(keep)) {
      stop("`subset` doesn't select any rows.", call. = FALSE)
    }
    df[keep, , drop = FALSE]
  }
  if (is.data.frame(data_list)) {
    return(select_rows(data_list))
  }
  lapply(data_list, select_rows)
}

# Scaling edge or node metrics within each network removes differences in
# their average level between networks, so warn if networks are compared.
check_network_scaling <- function(scale, num_models, target_name, brms_formula) {
  if (scale != "none" && num_models > 1 && target_name != "global" &&
      grepl("bison_network", brms_formula, fixed = TRUE)) {
    warning(
      "With scale = \"", scale, "\", the network variable is scaled separately within each network, ",
      "which removes differences in its average between networks, so the effect of bison_network ",
      "on its level will be close to 0. Use scale = \"none\" to compare levels between networks.",
      call. = FALSE
    )
  }
  invisible(NULL)
}
