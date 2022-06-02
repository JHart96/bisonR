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
  # obj$formula <- formula
  # obj$edgemodel <- edgemodel
  obj$fit <- fit
  # obj$model_data <- model_data
  # obj$design_matrices <- design_matrices
  # obj$chain <- chain
  # obj$edge_samples <- edge_samples
  # obj$dyad_ids <- dyad_ids
  class(obj) <- "nodal_model"
  obj
}
