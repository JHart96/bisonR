#' Run prior predictive check for a model specification
#'
#' @param formula Formula for the model
#' @param data Data for the model
#' @param model_type Type of model to fit
#' @param priors List of priors for the model
#' @param plot_type Type of predictive plot to show ("density" or "marginal")
#' @param options List of additional information needed for the predictive plot. Check the examples for full details.
#'
#' @export
prior_predictive_check <- function(formula, data, model_type, plot_type="density", options=NULL, priors=NULL) {
  if (model_type %in% c("binary", "count", "binary_conjugate", "count_conjugate")) {
    directed <- options$directed
    if (is.null(directed)) {
      directed=FALSE
    }
    bison_model.prior_predictive_check(formula, data, model_type, directed, priors, plot_type)
  }
}

bison_model.prior_predictive_check <- function(formula, data, model_type, directed, priors, plot_type) {
  fit <- bison_model(
    formula=formula,
    data=data,
    data_type=model_type,
    directed=directed,
    priors=priors,
    priors_only=TRUE
  )
  plot_predictions(fit, num_draws=20, draw_data=FALSE, type=plot_type)
}
