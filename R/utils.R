#' Traceplot of MCMC chains
#'
#' @param obj
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot_trace <- function(obj, ...) {
  if (class(obj) == "edge_model") {
    plot_trace.edge_model(obj, ...)
  } else if (class(obj) == "dyadic_model") {
    plot_trace.dyadic_model(obj, ...)
  } else if (class(obj) == "nodal_model") {
    plot_trace.nodal_model(obj, ...)
  }
}

#' Posterior predictive checks
#'
#' @param obj
#'
#' @return
#' @export
#'
#' @examples
plot_predictions <- function(obj, ...) {
  if (class(obj) == "edge_model") {
    plot_predictions.edge_model(obj, ...)
  } else if (class(obj) == "dyadic_model") {
    plot_predictions.dyadic_model(obj, ...)
  } else if (class(obj) == "nodal_model") {
    plot_predictions.nodal_model(obj, ...)
  }
}

#' Calculate contrasts between parameters
#'
#' @param obj
#'
#' @return
#' @export
#'
#' @examples
get_contrasts <- function(obj, parameter_1, parameter_2, ci=0.90) {
  beta_samples <- obj$chain
  colnames(beta_samples) <- colnames(obj$model_data$design_fixed)
  beta_diff <- beta_samples[, parameter_1] - beta_samples[, parameter_2]
  beta_quant <- quantile(beta_diff, probs=c(0.5, 0.5 * (1 - ci), ci + 0.5 * (1 - ci)))
  return(round(beta_quant, 3))
}

build_stan_model <- function(model_name) {
  model_filepath <- system.file("stan", paste0(model_name, ".stan"), package="bisonR")
  model <- cmdstanr::cmdstan_model(model_filepath)
}
