#' Traceplot of MCMC chains
#'
#' @param obj A fitted S3 object for edge weight, dyadic regression, or nodal regression models.
#' @param ... Additional arguments to be passed to trace plots.
#'
#' @details Plots MCMC chains for a fitted model. If the model has too many parameters to plot, only the first 12
#' will be plotted. To change the parameters being plotted, pass the `par_ids` argument
#'
#' @export
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
#' @param obj A fitted S3 object for edge weight, dyadic regression, or nodal regression models.
#' @param ... Additional arguments to be passed to predictive plots.
#'
#' @details
#' Plots the densities of summary statistics of observed data against predictions made by a model. These plots
#' can be used to check the predictive performance of a model.
#'
#' @export
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
#' @param obj Fitted S3 dyadic or nodal model.
#' @param parameter_1 Name of parameter 1, as shown in the model summary.
#' @param parameter_2 Name of parameter 2, as shown in the model summary.
#' @param ci Credible interval width, calculated by quantiles.
#'
#' @return A vector of the median, lower bound, and upper bound of the posterior contrast between parameters 1 and 2.
#' @export
get_contrasts <- function(obj, parameter_1, parameter_2, ci=0.90) {
  beta_samples <- obj$chain
  colnames(beta_samples) <- colnames(obj$model_data$design_fixed)
  beta_diff <- beta_samples[, parameter_1] - beta_samples[, parameter_2]
  beta_quant <- quantile(beta_diff, probs=c(0.5, 0.5 * (1 - ci), ci + 0.5 * (1 - ci)))
  return(round(beta_quant, 3))
}

edge_transform <- function(obj) {
  edge_samples <- obj$edge_samples
  if (obj$data_type == "binary") {
    return(plogis(edge_samples))
  }
  if (obj$data_type == "count") {
    return(exp(edge_samples))
  }
  if (obj$data_type == "duration") {
    return(plogis(edge_samples))
  }
}

build_stan_model <- function(model_name) {
  model_filepath <- system.file("stan", paste0(model_name, ".stan"), package="bisonR")
  model <- cmdstanr::cmdstan_model(model_filepath, compile=FALSE)
  model$compile(dir=tempdir())
  return(model)
}
