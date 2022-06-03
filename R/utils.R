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
