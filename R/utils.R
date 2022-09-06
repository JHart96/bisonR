#' Traceplot of MCMC chains
#'
#' @param obj A fitted S3 object for edge weight, dyadic regression, or nodal regression models.
#' @param par_ids A vector of parameter IDs to be plotted.
#' @param ... Additional arguments to be passed to base R plots.
#'
#' @details Plots MCMC chains for a fitted model. If the model has too many parameters to plot, only the first 12
#' will be plotted. To change the parameters being plotted, pass the `par_ids` argument
#'
#' @export
plot_trace <- function(obj, par_ids=1:5, ...) {
  # Extract chains from objects
  if (is(obj, "bison_model")) {
    if (obj$data_type %in% c("binary_conjugate", "count_conjugate")) {
      stop("Trace plots not needed for conjugate models.")
    }
    chain <- obj$fit$draws("edge_weight")
  }
  # Set parameter IDs
  if (dim(chain)[3] < max(par_ids)) par_ids <- 1:dim(chain)[3]
  num_chains <- dim(chain)[2]
  for (par_id in par_ids) {
    for (chain_id in 1:num_chains) {
      if (chain_id == 1) {
        plot(chain[, chain_id, par_id], type="l", col=bison_colors[chain_id], ylab=paste0("Parameter ", par_id), ...)
      } else {
        lines(chain[, chain_id, par_id], type="l", col=bison_colors[chain_id], ...)
      }
    }
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
  if (is(obj, "bison_model")) {
    plot_predictions.bison_model(obj, ...)
  }
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
  model <- cmdstanr::cmdstan_model(model_filepath, compile=FALSE, stanc_options = list("O1"))
  model$compile(dir=tempdir())
  return(model)
}

bison_colors <- c("#2b9392", "#5aabaa", "#89c3c2", "#b8dbda")

col2rgba <- function(col, alpha) {
  rgba <- col2rgb(col, alpha=TRUE)
  rgba <- rgba/255
  rgba[4] <- alpha
  return(do.call(rgb, as.list(c(rgba))))
}
