#' Extract posterior samples of network metrics
#'
#' Draws networks from the posterior of a fitted edge model and calculates a
#' network metric for each one, giving a posterior distribution of the metric
#' that carries the uncertainty in the edge weights.
#'
#' @param obj A fitted `bison_model` object.
#' @param metric_string Name of the metric, made of a target and a metric
#'   separated by an underscore, e.g. `"node_strength"` or `"global_cv"`. See
#'   *Available metrics*.
#' @param num_draws Number of posterior draws to use.
#' @param standardise If `TRUE`, metrics are mean-centred within each posterior
#'   draw (so each row has mean 0).
#' @param transform Only used for `"edge_weight"`. If `TRUE` (default), edge
#'   weights are on the natural scale (probabilities or rates); if `FALSE` they're
#'   on the link scale (logit or log). All other metrics are always calculated
#'   from natural-scale edge weights. In models with an intercept or other
#'   effects, natural-scale edge weights are for an observation where those
#'   effects are 0.
#'
#' @section Available metrics:
#' All metrics are calculated on the *weighted* network unless they take a
#' threshold `[t]`, in which case edges with weight at or above `t` count as
#' present and the rest as absent.
#'
#' * `edge_weight`: the edge weights themselves.
#' * `node_strength`: sum of the weights of a node's edges (in and out edges
#'   combined in directed networks).
#' * `node_degree[t]`: number of edges with weight at or above `t`.
#' * `node_eigen`: eigenvector centrality. Edge direction is ignored in directed
#'   networks.
#' * `node_betweenness`: betweenness, using 1/weight as the distance (directed
#'   paths in directed networks).
#' * `node_closeness`: closeness, using 1/weight as the distance (paths leaving
#'   the node in directed networks).
#' * `node_clustering[t]`: local clustering coefficient of the thresholded network.
#' * `global_density`: mean edge weight of the dyads in the model.
#' * `global_cv`: coefficient of variation of edge weights (social differentiation).
#' * `global_std`: standard deviation of edge weights.
#' * `global_diameter`: diameter of the network.
#' * `global_clustering[t]`: global clustering coefficient of the thresholded
#'   network (0 if there are no connected triples).
#'
#' @section Nodes and dyads without data:
#' Metrics are calculated from the dyads in the model. Dyads that aren't in the
#' data have no edge, so they don't contribute to strength or degree and act as
#' infinitely long paths for betweenness and closeness. Nodes with no dyads at
#' all (kept from factor levels) therefore have fixed values, such as 0 strength
#' or `NaN` closeness, and a warning is given.
#'
#' @return A matrix with one row per posterior draw and one column per target:
#'   per dyad (in the order of [get_edgelist()]) for edge metrics, per node (in
#'   the order of the node factor levels, with nodes as column names) for node
#'   metrics, and a single column for global metrics.
#'
#' @examples
#' sim_data <- simulate_bison_model("binary", aggregated = TRUE)
#' fit_edge <- bison_model(
#'   (event | effort) ~ dyad(node_1_id, node_2_id),
#'   data = sim_data$df_sim,
#'   model_type = "binary_conjugate",
#'   priors = get_default_priors("binary_conjugate")
#' )
#' strength <- extract_metric(fit_edge, "node_strength", num_draws = 100)
#' apply(strength, 2, quantile, probs = c(0.05, 0.5, 0.95))
#'
#' cv <- extract_metric(fit_edge, "global_cv", num_draws = 100)
#' plot_metric(cv, "Social differentiation")
#' @export
extract_metric <- function(obj, metric_string, num_draws = 1000, standardise = FALSE, transform = TRUE) {
  if (is.null(obj$dyad_to_idx)) {
    stop("This model has no edge weights (it has no dyad() term).", call. = FALSE)
  }
  draw_ids <- sample_draw_ids(nrow(obj$edge_samples), num_draws)

  if (metric_string == "edge_weight") {
    edge_samples <- if (transform) edge_transform(obj) else obj$edge_samples
    metric_samples <- edge_samples[draw_ids, , drop = FALSE]
    colnames(metric_samples) <- obj$dyad_names
  } else {
    metric_fn <- get_metric_fn(metric_string)
    target_name <- stringr::str_split(metric_string, "_")[[1]][1]
    num_targets <- switch(target_name, edge = obj$num_dyads, node = obj$num_nodes, global = 1)

    net <- bison_base_igraph(obj)
    check_metric_network(obj, net, metric_string, target_name)
    edge_samples <- edge_transform(obj)
    metric_samples <- matrix(0, num_draws, num_targets)
    for (i in seq_len(num_draws)) {
      igraph::E(net)$weight <- edge_samples[draw_ids[i], ]
      metric_samples[i, ] <- metric_fn(net)
    }
    if (target_name == "node") {
      colnames(metric_samples) <- names(obj$node_to_idx)
    }
  }

  if (standardise) {
    if (ncol(metric_samples) < 2) {
      stop(
        "`standardise = TRUE` centres a metric across nodes or edges within each draw, ",
        "so it isn't available for global metrics.",
        call. = FALSE
      )
    }
    metric_samples <- metric_samples - rowMeans(metric_samples)
  }
  metric_samples
}

#' Plot the posterior distribution of a network metric
#'
#' Continuous metrics are shown as smoothed densities; metrics that only take
#' whole-number values (such as degree) are shown as bar charts.
#'
#' @param x A matrix of metric samples from [extract_metric()].
#' @param metric_name Label for the x axis.
#' @param ... Additional arguments passed to the plotting function.
#'
#' @examples
#' sim_data <- simulate_bison_model("binary", aggregated = TRUE)
#' fit_edge <- bison_model(
#'   (event | effort) ~ dyad(node_1_id, node_2_id),
#'   data = sim_data$df_sim,
#'   model_type = "binary_conjugate",
#'   priors = get_default_priors("binary_conjugate")
#' )
#' plot_metric(extract_metric(fit_edge, "global_cv", num_draws = 200), "CV of edge weights")
#' @export
plot_metric <- function(x, metric_name = "Network Metric", ...) {
  x <- as.matrix(x)
  for (i in seq_len(ncol(x))) {
    label <- metric_name
    if (ncol(x) > 1) {
      label <- paste(metric_name, if (!is.null(colnames(x))) colnames(x)[i] else i)
    }
    values <- x[, i]
    values <- values[is.finite(values)]
    if (length(values) == 0) {
      next
    }
    if (all(values == round(values)) || length(unique(values)) < 3) {
      graphics::barplot(
        discrete_distribution(values), col = bison_colors[1], border = "white",
        xlab = label, ylab = "Probability", main = "", ...
      )
    } else {
      fit_logspline <- logspline::logspline(values)
      xmin <- logspline::qlogspline(0.001, fit_logspline)
      xmax <- logspline::qlogspline(0.999, fit_logspline)
      graphics::curve(
        logspline::dlogspline(x, fit_logspline), from = xmin, to = xmax,
        lwd = 2, col = bison_colors[1], xlab = label, ylab = "Probability density", main = "",
        ...
      )
    }
  }
  invisible(NULL)
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
      threshold <- parse_threshold(metric_string)
      return(function(x) {
        return(igraph::strength(x, weights = 1 * (igraph::E(x)$weight >= threshold)))
      })
    }
    if (!is.na(stringr::str_match(metric_name, "^clustering\\[.*\\]$"))) {
      threshold <- parse_threshold(metric_string)
      return(function(x) {
        x_binary <- x
        edges_to_remove <- which(igraph::E(x)$weight < threshold)
        x_binary <- igraph::delete_edges(x_binary, edges_to_remove)
        clustering <- igraph::transitivity(x_binary, type="local")
        clustering[is.nan(clustering)] <- 0
        return(clustering)
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
    if (metric_name == "diameter") {
      return(function(net) igraph::diameter(net, weights = 1 / igraph::E(net)$weight))
    }
    if (!is.na(stringr::str_match(metric_name, "^clustering\\[.*\\]$"))) {
      threshold <- parse_threshold(metric_string)
      return(function(x) {
        x_binary <- x
        edges_to_remove <- which(igraph::E(x)$weight < threshold)
        x_binary <- igraph::delete_edges(x_binary, edges_to_remove)
        clustering <- igraph::transitivity(x_binary, type="global")
        if (is.nan(clustering)) clustering <- 0
        return(clustering)
      })
    }
  }
  stop(paste0("Network metric ", metric_string, " does not exist."))
}

# Read the threshold from metric names such as "node_degree[0.2]".
parse_threshold <- function(metric_string) {
  value <- regmatches(metric_string, regexpr("\\[[^]]*\\]", metric_string))
  threshold <- suppressWarnings(as.numeric(gsub("\\[|\\]", "", value)))
  if (length(threshold) != 1 || is.na(threshold)) {
    stop(
      "Couldn't read the threshold in \"", metric_string, "\". Give a number in square brackets, ",
      "e.g. \"node_degree[0.2]\".",
      call. = FALSE
    )
  }
  threshold
}

# Warn about network structures that make some metrics hard to interpret.
check_metric_network <- function(obj, net, metric_string, target_name) {
  if (target_name != "node") {
    return(invisible(NULL))
  }
  isolated <- names(obj$node_to_idx)[igraph::degree(net) == 0]
  if (length(isolated) > 0) {
    warning(
      length(isolated), " node(s) have no dyads in the data (", paste(utils::head(isolated, 5), collapse = ", "),
      if (length(isolated) > 5) ", ..." else "", "), so their values of ", metric_string,
      " are fixed (e.g. 0 strength, or NaN closeness) rather than estimated.",
      call. = FALSE
    )
  }
  if (metric_string %in% c("node_eigen", "node_closeness") && igraph::components(net)$no > 1) {
    warning(
      "The network has more than one connected component. ",
      if (metric_string == "node_eigen") {
        "Eigenvector centrality is then concentrated in one component, which can change between draws."
      } else {
        "Closeness is calculated within each component, so values can't be compared across components."
      },
      call. = FALSE
    )
  }
  invisible(NULL)
}

# Proportion of draws taking each value of a discrete metric. Whole-number
# metrics get a bar for every value in their range; other discrete values
# (e.g. clustering of 0 or 1/3) get a bar for each value seen.
discrete_distribution <- function(values) {
  bar_levels <- if (all(values == round(values))) seq(min(values), max(values)) else sort(unique(values))
  counts <- table(factor(values, levels = bar_levels))
  counts / sum(counts)
}
