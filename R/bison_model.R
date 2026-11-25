#' Fit an edge weight model
#'
#' Fits a BISoN edge weight model, which estimates the strength of the social
#' connection between each pair of nodes (dyad) *with uncertainty*, from
#' observations of social events and the effort spent observing each dyad.
#'
#' @param formula A formula of the form `(event | effort) ~ dyad(node_1, node_2)`.
#'   The left-hand side gives the column of social events (0/1 or counts) and the
#'   observation effort (a column, or a single number), and the right-hand side
#'   gives the terms that explain the events. See *Formula syntax* below.
#' @param data A data frame with one row per observation (disaggregated) or one
#'   row per dyad (aggregated).
#' @param model_type The type of edge model: `"binary"` (events are present/absent
#'   in each sampling period, edge weights are probabilities), `"count"` (events
#'   are counts, edge weights are rates per unit effort), or the fast conjugate
#'   versions `"binary_conjugate"` and `"count_conjugate"`, which only support a
#'   `dyad()` term.
#' @param directed `TRUE` or `FALSE`: whether the network is directed. For
#'   undirected networks, rows for A-B and B-A refer to the same dyad.
#' @param partial_pooling Whether to estimate the spread of edge weights from the
#'   data, so that information is shared between edges.
#' @param zero_inflated Whether to use a zero-inflated model for excess zeroes.
#' @param priors A list of priors, in the format returned by [get_default_priors()].
#' @param refresh How often the sampler prints progress (0 for silent).
#' @param mc_cores Number of chains to run in parallel.
#' @param iter_sampling Number of post-warmup draws per chain.
#' @param iter_warmup Number of warmup iterations per chain.
#' @param priors_only If `TRUE`, ignore the data and sample from the priors (used
#'   for prior predictive checks).
#' @param chains Number of MCMC chains.
#' @param ... Further arguments passed to the `$sample()` method of
#'   [cmdstanr::CmdStanModel], such as `seed` or `adapt_delta`.
#'
#' @section Formula syntax:
#' * `dyad(node_1, node_2)`: estimate an edge weight for each dyad. The two
#'   columns hold node names or IDs. If they are both factors with the same
#'   levels, nodes that never appear in the data are kept.
#' * `x`: a fixed effect of column `x`. Factors get one parameter per level.
#' * `(1 | g)`: a random intercept for each level of column `g`.
#' * `1` or `0`: include or exclude a global intercept (excluded by default).
#'
#' Including extra terms changes what the edge weights mean: they become the
#' dyad's tendency to interact *after* accounting for those effects.
#'
#' @return A `bison_model` object. Use [summary.bison_model()] or [get_edgelist()] to see the
#'   edge weights, [plot_network()] to visualise them, [plot_trace()] and
#'   [plot_predictions()] for diagnostics, and [extract_metric()] or [bison_brm()]
#'   for downstream analyses.
#'
#' @seealso [get_default_priors()], [prior_check()], [summary.bison_model()]
#'
#' @examples
#' sim_data <- simulate_bison_model("binary", aggregated = TRUE)
#' df <- sim_data$df_sim
#'
#' # Fast conjugate model, no MCMC needed
#' fit_edge <- bison_model(
#'   (event | duration) ~ dyad(node_1_id, node_2_id),
#'   data = df,
#'   model_type = "binary_conjugate",
#'   priors = get_default_priors("binary_conjugate")
#' )
#' summary(fit_edge)
#'
#' # Full model fitted with Stan, which allows extra effects
#' fit_edge <- bison_model(
#'   (event | duration) ~ dyad(node_1_id, node_2_id),
#'   data = df,
#'   model_type = "binary",
#'   priors = get_default_priors("binary")
#' )
#' @export
bison_model <- function(formula, data,
                        model_type = c("binary", "count", "binary_conjugate", "count_conjugate"),
                        directed = FALSE, partial_pooling = FALSE, zero_inflated = FALSE,
                        priors = NULL, refresh = 0, mc_cores = 4, iter_sampling = 1000,
                        iter_warmup = 1000, priors_only = FALSE, chains = 4, ...) {
  model_type <- match.arg(model_type)
  full_model_type <- model_type

  # If user-specified priors haven't been set, use the defaults
  if (is.null(priors)) {
    message(
      "No priors set by user, using default priors instead. ",
      "We recommend setting and checking priors explicitly for reliable inference."
    )
    priors <- get_default_priors(model_type)
  }
  priors <- validate_priors(priors, model_type)

  use_conjugate_model <- grepl("_conjugate$", model_type)
  model_type <- sub("_conjugate$", "", model_type)

  model_info <- get_bison_model_data(formula, data, directed, model_type)
  model_spec <- model_info$model_spec

  if (use_conjugate_model) {
    if (is.null(model_spec$node_1_name)) {
      stop("Conjugate models need a dyad(node_1, node_2) term.", call. = FALSE)
    }
    if (model_spec$intercept || length(model_spec$fixed) > 0 || length(model_spec$random) > 0) {
      stop(
        "Conjugate models only support dyad() terms. Use model_type = \"", model_type,
        "\" to include intercepts, fixed effects or random effects.",
        call. = FALSE
      )
    }
    if (partial_pooling || zero_inflated) {
      stop("Partial pooling and zero inflation aren't available for conjugate models.", call. = FALSE)
    }
  }

  # Set the priors in model data
  prior_parameters <- extract_prior_parameters(priors)
  model_data <- c(model_info$model_data, prior_parameters)
  model_data$priors_only <- as.integer(priors_only)
  model_data$partial_pooling <- as.integer(partial_pooling)
  model_data$zero_inflated <- as.integer(zero_inflated)

  if (use_conjugate_model) {
    model <- fit_conjugate_model(
      model_type, model_data,
      num_samples = iter_sampling * chains,
      priors_only = priors_only
    )
    chain <- model$chain
    event_preds <- model$event_preds
    log_lik <- model$log_lik
    fit <- NULL
    stan_model <- NULL
  } else {
    stan_model <- build_stan_model(model_type)

    fit <- stan_model$sample(
      data = model_data,
      refresh = refresh,
      chains = chains,
      parallel_chains = mc_cores,
      iter_sampling = iter_sampling,
      iter_warmup = iter_warmup,
      ...
    )

    if (model_info$num_dyads > 0) {
      chain <- as_plain_matrix(fit$draws("edge_weight", format = "matrix"))
    } else {
      chain <- NULL
    }
    event_preds <- as_plain_matrix(fit$draws("event_pred", format = "matrix"))
    log_lik <- as_plain_matrix(fit$draws("log_lik", format = "matrix"))
  }

  obj <- list(
    chain = chain,
    event_preds = event_preds,
    edge_samples = chain,
    num_nodes = model_info$num_nodes,
    num_dyads = model_info$num_dyads,
    node_to_idx = model_info$node_to_idx,
    dyad_to_idx = model_info$dyad_to_idx,
    dyad_names = model_info$dyad_names,
    model_info = model_info,
    directed = directed,
    fit = fit,
    formula = formula,
    model_type = model_type,
    full_model_type = full_model_type,
    priors = priors,
    model_data = model_data,
    input_data = data,
    stan_model = stan_model,
    conjugate = use_conjugate_model,
    log_lik = log_lik,
    bisonR_version = as.character(utils::packageVersion("bisonR"))
  )
  class(obj) <- "bison_model"
  obj
}

# Convert a posterior draws_matrix into a plain numeric matrix, so that
# subsetting behaves like ordinary R matrices.
as_plain_matrix <- function(x) {
  out <- matrix(as.numeric(x), nrow = nrow(x), ncol = ncol(x))
  colnames(out) <- colnames(x)
  out
}

#' @export
print.bison_model <- function(x, ...) {
  cat(
    "BISoN edge model (", x$full_model_type, ")\n",
    "Formula: ", deparse1(x$formula), "\n",
    "Nodes: ", x$num_nodes, ", dyads: ", x$num_dyads,
    ", directed: ", x$directed, "\n",
    "Use summary() to see the estimated edge weights.\n",
    sep = ""
  )
  invisible(x)
}

#' Summarise an edge weight model
#'
#' @param object A fitted `bison_model` object.
#' @param ci Width of the credible interval, e.g. `0.9` for a 90% interval.
#' @param transform If `TRUE` (the default), edge weights are reported on the
#'   natural scale (probabilities for binary models, rates for count models).
#'   If `FALSE`, they're on the internal link scale (logit or log).
#' @param ... Not used.
#'
#' @return A `summary.bison_model` object. Printing it shows, for each dyad, the
#'   posterior median edge weight and the lower and upper bounds of the
#'   credible interval.
#'
#' @examples
#' sim_data <- simulate_bison_model("binary", aggregated = TRUE)
#' fit_edge <- bison_model(
#'   (event | duration) ~ dyad(node_1_id, node_2_id),
#'   data = sim_data$df_sim,
#'   model_type = "binary_conjugate",
#'   priors = get_default_priors("binary_conjugate")
#' )
#' summary(fit_edge, ci = 0.95)
#' @export
summary.bison_model <- function(object, ci = 0.90, transform = TRUE, ...) {
  summary_obj <- list()

  summary_obj$description <- paste0(
    "=== Fitted BISoN edge model ===",
    "\nData type: ", object$model_type,
    "\nFormula: ", deparse1(object$formula),
    "\nNumber of nodes: ", object$num_nodes,
    "\nNumber of dyads: ", object$num_dyads,
    "\nDirected: ", object$directed,
    "\n=== Edge list summary ===\n"
  )

  summary_obj$edgelist <- get_edgelist(object, ci = ci, transform = transform)
  summary_obj$dyad_names <- object$dyad_names

  class(summary_obj) <- "summary.bison_model"
  summary_obj
}

#' Print an edge model summary
#'
#' @param x A `summary.bison_model` object.
#' @param ... Additional arguments passed to [print()].
#'
#' @seealso [summary.bison_model()]
#' @keywords internal
#' @export
print.summary.bison_model <- function(x, ...) {
  cat(x$description)
  summary_matrix <- as.matrix(x$edgelist[, 3:5])
  rownames(summary_matrix) <- x$dyad_names
  print(summary_matrix, ...)
  invisible(x)
}

#' Get the edge list with uncertainty
#'
#' @param obj A fitted `bison_model` object.
#' @param ci Width of the credible interval, e.g. `0.9` for a 90% interval.
#' @param transform If `TRUE` (the default), edge weights are on the natural
#'   scale (probabilities for binary models, rates for count models). If `FALSE`,
#'   they're on the internal link scale (logit or log).
#'
#' @return A data frame with one row per dyad and columns `node_1`, `node_2`,
#'   `median`, and the lower and upper bounds of the credible interval (named
#'   after their quantiles, e.g. `5%` and `95%`).
#'
#' @examples
#' sim_data <- simulate_bison_model("binary", aggregated = TRUE)
#' fit_edge <- bison_model(
#'   (event | duration) ~ dyad(node_1_id, node_2_id),
#'   data = sim_data$df_sim,
#'   model_type = "binary_conjugate",
#'   priors = get_default_priors("binary_conjugate")
#' )
#' head(get_edgelist(fit_edge))
#' @export
get_edgelist <- function(obj, ci = 0.9, transform = TRUE) {
  if (is.null(obj$dyad_to_idx)) {
    stop("This model has no edge weights (it has no dyad() term).", call. = FALSE)
  }
  check_ci(ci)
  node_names <- names(obj$node_to_idx)
  lb <- 0.5 * (1 - ci)
  ub <- 1 - lb

  edge_samples <- obj$edge_samples
  if (transform) {
    edge_samples <- edge_transform(obj)
  }

  quantiles <- apply(edge_samples, 2, quantile, probs = c(0.5, lb, ub), names = FALSE)

  edgelist <- data.frame(
    node_1 = node_names[obj$dyad_to_idx[, 1]],
    node_2 = node_names[obj$dyad_to_idx[, 2]],
    median = round(quantiles[1, ], 3),
    lb = round(quantiles[2, ], 3),
    ub = round(quantiles[3, ], 3),
    row.names = NULL
  )
  colnames(edgelist)[4] <- paste0(format_percent(lb), "%")
  colnames(edgelist)[5] <- paste0(format_percent(ub), "%")
  edgelist
}

check_ci <- function(ci) {
  if (!is.numeric(ci) || length(ci) != 1 || ci <= 0 || ci >= 1) {
    stop("`ci` must be a single number between 0 and 1, e.g. 0.9.", call. = FALSE)
  }
}

format_percent <- function(p) {
  as.character(round(p * 100, 2))
}

#' Draw samples of the edge list from the posterior
#'
#' @param obj A fitted `bison_model` object.
#' @param num_draws Number of posterior draws to return.
#' @param transform If `TRUE`, edge weights are on the natural scale; if `FALSE`
#'   (the default, for backwards compatibility) they're on the link scale.
#'
#' @return A data frame with columns `node_1` and `node_2` followed by one
#'   column per posterior draw (`draw.1`, `draw.2`, ...).
#'
#' @examples
#' sim_data <- simulate_bison_model("binary", aggregated = TRUE)
#' fit_edge <- bison_model(
#'   (event | duration) ~ dyad(node_1_id, node_2_id),
#'   data = sim_data$df_sim,
#'   model_type = "binary_conjugate",
#'   priors = get_default_priors("binary_conjugate")
#' )
#' head(draw_edgelist_samples(fit_edge, num_draws = 3, transform = TRUE))
#' @export
draw_edgelist_samples <- function(obj, num_draws, transform = FALSE) {
  if (is.null(obj$dyad_to_idx)) {
    stop("This model has no edge weights (it has no dyad() term).", call. = FALSE)
  }
  node_names <- names(obj$node_to_idx)
  edge_samples <- if (transform) edge_transform(obj) else obj$chain
  draw_ids <- sample_draw_ids(nrow(edge_samples), num_draws)
  draws <- t(edge_samples[draw_ids, , drop = FALSE])
  colnames(draws) <- paste0("draw.", seq_len(num_draws))
  data.frame(
    node_1 = factor(node_names[obj$dyad_to_idx[, 1]], levels = node_names),
    node_2 = factor(node_names[obj$dyad_to_idx[, 2]], levels = node_names),
    draws,
    row.names = NULL,
    check.names = FALSE
  )
}

# Choose posterior draws at random: without replacement where possible,
# with replacement if more draws are requested than are available.
sample_draw_ids <- function(num_available, num_draws) {
  if (num_draws < 1) {
    stop("`num_draws` must be at least 1.", call. = FALSE)
  }
  sample.int(num_available, size = num_draws, replace = num_draws > num_available)
}
