#' Draw samples from the posterior node metric for a fitted edge weight model
#'
#' @param obj Fitted S3 edge weight model.
#' @param metric_name Character name of node-level metric function.
#' @param num_draws Number of posterior draws to use.
#' @param standardise TRUE/FALSE indicating whether to mean-center the metrics within sample.
#'
#' @return A matrix of metric samples where each column corresponds to a node, and each row corresponds to a posterior draw.
#' @export
draw_node_metric_samples <- function(obj, metric_name, num_draws=1000, standardise=FALSE) {
  edgelist_samples <- draw_edgelist_samples(obj, num_draws)
  net <- igraph::graph_from_edgelist(as.matrix(edgelist_samples[, 1:2]), directed = FALSE)
  metric_fn <- get_metric_fn(metric_name)

  metric_samples <- matrix(0, num_draws, obj$num_nodes)
  colnames(metric_samples) <- names(obj$node_to_idx)
  for (i in 1:num_draws) {
    if (obj$data_type == "binary") {
      igraph::E(net)$weight <- plogis(edgelist_samples[, 2 + i])
    }
    if (obj$data_type == "count") {
      igraph::E(net)$weight <- exp(edgelist_samples[, 2 + i])
    }
    if (obj$data_type == "duration") {
      igraph::E(net)$weight <- plogis(edgelist_samples[, 2 + i])
    }
    metric_sample <- metric_fn(net)
    if (standardise) {
      metric_samples[i, ] <- metric_sample - mean(metric_sample)
    } else {
      metric_samples[i, ] <- metric_sample
    }
  }
  metric_samples
}

#' Draw samples from posterior network-level metrics for a fitted edge weight model
#'
#' @param obj Fitted S3 edge weight model.
#' @param metric_name Character name of network-level metric function.
#' @param num_draws Number of posterior draws to use.
#' @param standardise TRUE/FALSE indicating whether to mean-center the metrics within sample.
#'
#' @return A matrix of metric samples where each column corresponds to a node, and each row corresponds to a posterior draw.
#' @export
draw_network_metric_samples <- function(obj, metric_name, num_draws=1000, standardise=FALSE) {
  edgelist_samples <- draw_edgelist_samples(obj, num_draws)
  net <- igraph::graph_from_edgelist(as.matrix(edgelist_samples[, 1:2]), directed = FALSE)
  metric_fn <- get_metric_fn(metric_name)

  metric_samples <- rep(0, num_draws)
  for (i in 1:num_draws) {
    if (obj$data_type == "binary") {
      igraph::E(net)$weight <- plogis(edgelist_samples[, 2 + i])
    }
    if (obj$data_type == "count") {
      igraph::E(net)$weight <- exp(edgelist_samples[, 2 + i])
    }
    if (obj$data_type == "duration") {
      igraph::E(net)$weight <- plogis(edgelist_samples[, 2 + i])
    }
    metric_samples[i] <- metric_fn(net)
  }
  if (standardise) {
    metric_samples <- metric_samples - mean(metric_samples)
  }
  obj <- list()
  class(obj) <- "network_metric_samples"
  obj$samples <- metric_samples
  obj$metric_name <- metric_name
  obj
}

"[.network_metric_samples" <- function(x, i) {
  return(unclass(x)[i])
}

#' Print a network metric samples object
#'
#' @param x An S3 network metric samples object
#' @param ci Credible interval to use for summary
#' @param ... Additional parameters to be passed to print
#'
#' @export
print.network_metric_samples <- function(x, ci=0.9, ...) {
  qt <- quantile(x$samples, probs=c(0.5, 0.5 * (1 - ci), ci + 0.5 * (1 - ci)))
  print(qt)
}

#' Plot network metric samples
#'
#' @param x Network metric samples object
#' @param ... Additional parameters to be passed to plot function
#'
#' @export
plot.network_metric_samples <- function(x, ...) {
  fit_logspline <- logspline::logspline(x$samples)
  xmin <- logspline::qlogspline(0.001, fit_logspline)
  xmax <- logspline::qlogspline(0.999, fit_logspline)
  curve(
    logspline::dlogspline(x, fit_logspline), xlim=c(xmin, xmax),
    lwd=2, col=bison_colors[1], xlab=x$metric_name, ylab="Probability density", main="", ...
  )
}

get_metric_fn <- function(metric_name) {
  if (metric_name == "strength") {
    return(igraph::strength)
  }
  if (metric_name == "eigenvector") {
    return(function(x) igraph::eigen_centrality(x)$vector)
  }
  if (metric_name == "betweenness") {
    return(function(x) igraph::betweenness(x, weights=1/igraph::E(x)$weight))
  }
  if (!is.na(stringr::str_match(metric_name, "^degree\\[.*\\]$"))) {
    threshold <- as.numeric(str_split(metric_name, "\\[|\\]")[[1]][2])
    return(function(x) {
      return(igraph::strength(x, weights=1 * (igraph::E(x)$weight < threshold)))
    })
  }
  if (metric_name == "weighted_density") {
    return(function(net) mean(igraph::E(net)$weight))
  }
  if (metric_name == "social_differentiation") {
    return(function(net) sd(igraph::E(net)$weight)/mean(igraph::E(net)$weight))
  }
  if (metric_name == "standard_deviation") {
    return(function(net) sd(igraph::E(net)$weight))
  }
  stop("Network metric does not exist.")
}
