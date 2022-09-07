#' Extract samples from posterior metrics for a fitted edge weight model
#'
#' @param obj Fitted S3 edge weight model.
#' @param metric_string Character name of metric function (e.g. node_strength).
#' @param num_draws Number of posterior draws to use.
#' @param standardise TRUE/FALSE indicating whether to mean-center the metrics within sample.
#'
#' @return A matrix of metric samples where each column corresponds to a node, and each row corresponds to a posterior draw.
#' @export
extract_metric <- function(obj, metric_string, num_draws=1000, standardise=FALSE) {
  if (metric_string == "edge_weight") {
    metric_samples <- matrix(as.numeric(obj$edge_samples), ncol=obj$num_dyads)
    if (standardise) {
      return(metric_samples - mean(metric_samples))
    } else {
      return(metric_samples)
    }
  }

  edgelist_samples <- draw_edgelist_samples(obj, num_draws)
  net <- igraph::graph_from_edgelist(as.matrix(edgelist_samples[, 1:2]), directed = FALSE)

  metric_fn <- get_metric_fn(metric_string)
  target_name <- stringr::str_split(metric_string, "_")[[1]][1]

  if (target_name == "edge") {
    num_targets <- obj$num_dyads
  }
  if (target_name == "node") {
    num_targets <- obj$num_nodes
  }
  if (target_name == "global") {
    num_targets <- 1
  }

  metric_samples <- matrix(0, num_draws, num_targets)

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
    metric_samples[i, ] <- metric_fn(net)
  }
  if (standardise) {
    metric_samples <- metric_samples - mean(metric_samples)
  }
  metric_samples
}

#' Plot metric samples
#'
#' @param x Network metric samples object
#' @param metric_name Metric name to use for plot axis label
#' @param ... Additional parameters to be passed to plot function
#'
#' @export
plot_metric <- function(x, metric_name="Network Metric", ...) {
  plot_name <- ""
  for (i in 1:ncol(x)) {
    if (ncol(x) > 1) {
      plot_name <- i
    }
    fit_logspline <- logspline::logspline(x[, i])
    xmin <- logspline::qlogspline(0.001, fit_logspline)
    xmax <- logspline::qlogspline(0.999, fit_logspline)
    curve(
      logspline::dlogspline(x, fit_logspline), xlim=c(xmin, xmax),
      lwd=2, col=bison_colors[1], xlab=paste0(metric_name, " ", plot_name), ylab="Probability density", main="",
      ...
    )
  }
}

get_metric_fn <- function(metric_string) {
  target_name <- stringr::str_split(metric_string, "_")[[1]][1]
  metric_name <- stringr::str_split(metric_string, "_")[[1]][2]
  if (target_name == "node") {
    if (metric_name == "strength") {
      return(igraph::strength)
    }
    if (metric_name == "eigen") {
      return(function(x) igraph::eigen_centrality(x)$vector)
    }
    if (metric_name == "betweenness") {
      return(function(x) igraph::betweenness(x, weights=1/igraph::E(x)$weight))
    }
    if (metric_name == "closeness") {
      return(function(x) igraph::closeness(x, weights=1/igraph::E(x)$weight))
    }
    if (!is.na(stringr::str_match(metric_name, "^degree\\[.*\\]$"))) {
      threshold <- as.numeric(str_split(metric_name, "\\[|\\]")[[1]][2])
      return(function(x) {
        return(igraph::strength(x, weights=1 * (igraph::E(x)$weight < threshold)))
      })
    }
  }
  if (target_name == "global") {
    if (metric_name == "density") {
      return(function(net) mean(igraph::E(net)$weight))
    }
    if (metric_name == "cv") {
      return(function(net) sd(igraph::E(net)$weight)/mean(igraph::E(net)$weight))
    }
    if (metric_name == "std") {
      return(function(net) sd(igraph::E(net)$weight))
    }
  }
  stop(paste0("Network metric ", metric_string, " does not exist."))
}
