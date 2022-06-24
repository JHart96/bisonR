#' Draw samples from the posterior node metric for a fitted edge weight model
#'
#' @param obj
#' @param metric_name
#'
#' @return
#' @export
#'
#' @examples
draw_node_metric_samples <- function(obj, metric_name, num_draws=1000, standardise=TRUE) {
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

get_metric_fn <- function(metric_name) {
  if (metric_name == "strength") {
    return(igraph::strength)
  }
  return(NULL)
}
