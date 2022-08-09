#' Run prior predictive check for a model specification
#'
#' @param formula Formula for the model
#' @param data Data for the model
#' @param model_type Type of model to fit
#' @param priors List of priors for the model
#'
#' @export
prior_predictive_check <- function(formula, data, model_type, options=NULL, priors=NULL) {
  if (model_type %in% c("binary", "count", "binary_conjugate", "count_conjugate")) {
    directed <- options$directed
    if (is.null(directed)) {
      directed=FALSE
    }
    edge_model.prior_predictive_check(formula, data, model_type, directed, priors)
  }
  if (model_type %in% c("dyadic_regression")) {
    edgemodel <- options$edgemodel
    mm <- options$mm
    if (is.null(mm)) {
      mm <- TRUE
    }
    if (is.null(edgemodel)) {
      stop("Predictive checks for dyadic regressions require an edge model and to be provided.")
    }
    dyadic_regression.prior_predictive_check(formula, data, edgemodel, mm, priors)
  }
  if(model_type %in% c("nodal_regression")) {
    edgemodel <- options$edgemodel
    if (is.null(edgemodel)) {
      stop("Predictive checks for nodal regressions require an edge model to be provided.")
    }
    nodal_regression.prior_predictive_check(formula, data, edgemodel, priors)
  }
}

edge_model.prior_predictive_check <- function(formula, data, model_type, directed, priors) {
  fit <- edge_model(
    formula=formula,
    data=data,
    data_type=model_type,
    directed=directed,
    priors=priors,
    priors_only=TRUE
  )
  plot_predictions(fit, num_draws=20, draw_data=FALSE, type="density")
}

dyadic_regression.prior_predictive_check <- function(formula, data, edgemodel, mm, priors) {
  fit <- dyadic_regression(
    formula=formula,
    edgemodel=edgemodel,
    df=data,
    mm=mm,
    priors=priors,
    priors_only=TRUE
  )
  plot_predictions(fit, num_draws=20, draw_data=FALSE, type="density")
}

nodal_regression.prior_predictive_check <- function(formula, data, edgemodel, priors) {
  fit <- nodal_regression(
    formula=formula,
    edgemodel=edgemodel,
    df=data,
    priors=priors,
    priors_only=TRUE
  )
  plot_predictions(fit, num_draws=20, draw_data=FALSE, type="density")
}


