#' Run leave-one-out model comparison
#'
#' @param models List of fitted models to compare
#' @param ... Additional arguments to be passed to the `loo_compare` function.
#'
#' @return The output of the `loo_compare` function for the fitted models.
#' @export
model_comparison <- function(models, ...) {
  loos <- lapply(models, function(model) loo::loo(model$log_lik))
  loo::loo_compare(loos, ...)
}
