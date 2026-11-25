#' Compare models with leave-one-out cross-validation
#'
#' Compares fitted models by their expected predictive accuracy, estimated with
#' Pareto-smoothed importance sampling leave-one-out cross-validation (PSIS-LOO)
#' from the loo package, and returns model weights.
#'
#' @param models A list of fitted models, ideally named. All models must be fitted
#'   to the same data (the same rows, in the same order).
#' @param ... Additional arguments passed to [loo::loo_model_weights()], such as
#'   `method = "pseudobma"`.
#'
#' @details If the models have different numbers of posterior draws (for example
#'   because a chain failed), the same number of draws is taken at random from
#'   each model.
#'
#' @return The model weights returned by [loo::loo_model_weights()].
#' @examplesIf !is.null(tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL))
#' # Test for non-random edge weights by comparing a model with edge weights
#' # against a null model where every dyad has the same rate of events.
#' sim_data <- simulate_bison_model("binary", aggregated = FALSE)
#' priors <- get_default_priors("binary")
#' fit_edge <- bison_model((event | effort) ~ dyad(node_1_id, node_2_id),
#'                         data = sim_data$df_sim, model_type = "binary", priors = priors)
#' fit_null <- bison_model((event | effort) ~ 1,
#'                         data = sim_data$df_sim, model_type = "binary", priors = priors)
#' model_comparison(list(non_random = fit_edge, random = fit_null))
#' @export
model_comparison <- function(models, ...) {
  if (!is.list(models) || length(models) < 2) {
    stop("`models` must be a list of at least two fitted models.", call. = FALSE)
  }
  log_liks <- lapply(models, function(model) model$log_lik)
  if (any(vapply(log_liks, is.null, logical(1)))) {
    stop("All models must contain log-likelihoods (fitted with bison_model()).", call. = FALSE)
  }

  num_obs <- vapply(log_liks, ncol, integer(1))
  if (length(unique(num_obs)) > 1) {
    stop(
      "Models must be fitted to the same data, but they have different numbers of rows (",
      paste(num_obs, collapse = ", "), ").",
      call. = FALSE
    )
  }

  num_draws <- vapply(log_liks, nrow, integer(1))
  if (length(unique(num_draws)) > 1) {
    common <- min(num_draws)
    message(
      "Models have different numbers of posterior draws (", paste(num_draws, collapse = ", "),
      "), so ", common, " draws from each are used."
    )
    log_liks <- lapply(log_liks, function(x) x[sort(sample.int(nrow(x), common)), , drop = FALSE])
  }

  loos <- lapply(log_liks, function(x) loo::loo(x, r_eff = NA))
  results <- loo::loo_model_weights(loos, ...)
  if (!is.null(names(models))) {
    names(results) <- names(models)
  }
  results
}
