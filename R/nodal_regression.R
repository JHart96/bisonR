#' Fit a nodel regression model
#'
#' @param formula
#' @param edgemodel
#' @param df
#' @param mc_cores
#'
#' @return
#' @export
#'
#' @examples
nodal_regression <- function(formula, edgemodel, df, mc_cores=4) {
  design_matrices <- build_design_matrix(formula, df)

  num_nodes <- edgemodel$num_nodes
  num_fixed <- ncol(design_matrices$X)

  metric_samples <- draw_node_metric_samples(edgemodel, "strength", 1000)
  metric_mu <- apply(metric_samples, 2, mean)
  metric_cov <- cov(metric_samples)

  model_data <- list(
    num_nodes=num_nodes,
    num_fixed=num_fixed,
    centrality_mu=metric_mu,
    centrality_cov=metric_cov,
    design_fixed=design_matrices$X
  )
  fit <- rstan::sampling(stanmodels$nodal_regression, model_data, cores=mc_cores)
  chain <- rstan::extract(fit)
  obj <- list()
  obj$formula <- formula
  obj$edgemodel <- edgemodel
  obj$fit <- fit
  obj$num_nodes <- num_nodes
  obj$model_data <- model_data
  obj$design_matrices <- design_matrices
  obj$chain <- chain
  obj$metric_samples <- metric_samples
  class(obj) <- "nodal_model"
  obj
}

#' Print information about a fitted nodal regression model
#'
#' @param obj
#'
#' @return
#' @export
#'
#' @examples
print.nodal_model <- function(obj) {
  coefficients <- t(apply(obj$chain[["beta_fixed"]], 2, function(x) quantile(x, probs=c(0.5, 0.05, 0.95))))
  rownames(coefficients) <- colnames(obj$design_matrices$X)
  coefficients <- round(coefficients, 3)
  cat(paste0(
    "=== Fitted nodal regression model ===\n",
    "Formula: ", format(obj$formula), "\n",
    "Number of nodes: ", obj$edgemodel$num_nodes, "\n",
    "=== Coefficient summary ==="
  ))
  print(coefficients)
}

#' Summary of a fitted nodal regression model
#'
#' @param obj
#'
#' @return
#' @export
#'
#' @examples
summary.nodal_model <- function(obj) {
  print(obj)
}

plot_trace.nodal_model <- function (obj, ...) {
  rstan::traceplot(obj$fit, ...)
}

plot_predictions.nodal_model <- function(obj, num_draws=20) {
  # Determine edge label
  xlab <- "Logit centrality"

  # Extract edge samples and predictions
  metric_samples <- obj$metric_samples
  metric_preds <- rstan::extract(obj$fit)$centrality_pred

  # Generate densities for edge sample and prediction
  sample_densities <- list()
  pred_densities <- list()
  for (i in 1:num_draws) {
    pred_densities[[i]] <- density(metric_samples[i, ])
    sample_densities[[i]] <- density(metric_preds[i, ])
  }

  # Set plot limits according to maximum density of samples
  xlim = c(
    min(sapply(sample_densities, function(x) min(x$x))),
    max(sapply(sample_densities, function(x) max(x$x))) * 1.1
  )
  ylim = c(
    min(sapply(sample_densities, function(x) min(x$y))),
    max(sapply(sample_densities, function(x) max(x$y))) * 1.1
  )

  # Plot densities for subsequent draws
  for (i in 1:num_draws) {
    if (i == 1) {
      plot(sample_densities[[i]], main="Observed vs predicted nodal metrics", xlab=xlab, col=rgb(0, 0, 0, 0.5), xlim=xlim, ylim=ylim)
    } else {
      lines(sample_densities[[i]], col=rgb(0, 0, 0, 0.5))
    }
    lines(pred_densities[[i]], col=rgb(0, 0, 1, 0.5))
  }
}
