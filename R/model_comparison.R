#' Run leave-one-out model comparison
#'
#' @param models List of fitted models to compare
#'
#' @return The output of the `loo_compare` function for the fitted models.
#' @export
model_comparison <- function(models) {
  loos <- lapply(models, function(model) loo::loo(model$log_lik, r_eff=NA))

  results_matrix <- loo::loo_model_weights(loos)
  if (!is.null(names(models))) {
    names(results_matrix) <- names(models)
  }
  results_matrix
}
