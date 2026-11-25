#' Run leave-one-out model comparison
#'
#' @param models List of fitted models to compare
#' @param ... Additional arguments to be passed to `loo::loo_model_weights`
#'
#' @return The output of the `loo_compare` function for the fitted models.
#' @examples
#' # Test for non-random edge weights by comparing a model with edge weights
#' # against a null model where every dyad has the same rate of events.
#' sim_data <- simulate_bison_model("binary", aggregated = FALSE)
#' priors <- get_default_priors("binary")
#' fit_edge <- bison_model((event | duration) ~ dyad(node_1_id, node_2_id),
#'                         data = sim_data$df_sim, model_type = "binary", priors = priors)
#' fit_null <- bison_model((event | duration) ~ 1,
#'                         data = sim_data$df_sim, model_type = "binary", priors = priors)
#' model_comparison(list(non_random = fit_edge, random = fit_null))
#' @export
model_comparison <- function(models, ...) {
  loos <- lapply(models, function(model) loo::loo(model$log_lik, r_eff=NA))

  results_matrix <- loo::loo_model_weights(loos, ...)
  if (!is.null(names(models))) {
    names(results_matrix) <- names(models)
  }
  results_matrix
}
