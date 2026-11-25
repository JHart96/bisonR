#' Trace plots of MCMC chains
#'
#' Plots the MCMC chains for some of the edge weights of a fitted model. Well-mixed
#' chains look like overlapping "fuzzy caterpillars" with no trends.
#'
#' @param obj A fitted `bison_model` object (fitted with Stan, not a conjugate
#'   model).
#' @param par_ids Indices of the edge weights to plot.
#' @param ... Additional arguments passed to [plot()].
#'
#' @examples
#' sim_data <- simulate_bison_model("binary", aggregated = TRUE)
#' fit_edge <- bison_model(
#'   (event | duration) ~ dyad(node_1_id, node_2_id),
#'   data = sim_data$df_sim, model_type = "binary",
#'   priors = get_default_priors("binary")
#' )
#' plot_trace(fit_edge, par_ids = 1:4)
#' @export
plot_trace <- function(obj, par_ids = 1:5, ...) {
  if (!inherits(obj, "bison_model")) {
    stop("plot_trace() needs a fitted bison_model object.", call. = FALSE)
  }
  if (isTRUE(obj$conjugate)) {
    stop(
      "Conjugate models are sampled exactly rather than by MCMC, so trace plots aren't needed.",
      call. = FALSE
    )
  }
  if (obj$num_dyads == 0) {
    stop("This model has no edge weights to plot.", call. = FALSE)
  }
  chain <- obj$fit$draws("edge_weight")
  par_ids <- par_ids[par_ids <= dim(chain)[3]]
  if (length(par_ids) == 0) {
    par_ids <- seq_len(min(5, dim(chain)[3]))
  }
  num_chains <- dim(chain)[2]

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  num_cols <- min(length(par_ids), 2)
  graphics::par(mfrow = c(ceiling(length(par_ids) / num_cols), num_cols), mar = c(4, 4, 2, 1))

  palette <- grDevices::colorRampPalette(c(bison_colors[1], "#e0a458", bison_highlight, "#5b5f97"))(max(num_chains, 4))
  for (par_id in par_ids) {
    values <- matrix(chain[, , par_id], ncol = num_chains)
    graphics::plot(
      NULL, xlim = c(1, nrow(values)), ylim = range(values),
      xlab = "Iteration", ylab = "Edge weight (link scale)",
      main = obj$dyad_names[par_id], ...
    )
    for (chain_id in seq_len(num_chains)) {
      graphics::lines(values[, chain_id], col = col2rgba(palette[chain_id], 0.7))
    }
  }
  invisible(NULL)
}

#' Posterior predictive checks
#'
#' Compares data simulated from a fitted model with the observed data. If the
#' model describes the data well, the observed data should look like a typical
#' simulated dataset.
#'
#' @param obj A fitted model.
#' @param ... Arguments passed to methods.
#'
#' @details
#' Three types of check are available, chosen with `type`:
#'
#' * `"density"`: the distribution of observed event *rates* per dyad (events
#'   divided by observation effort), compared with the same distribution in
#'   `num_draws` simulated datasets. Rates are used rather than raw event counts
#'   because raw counts mostly reflect how much each dyad was observed, which the
#'   model takes as given, so almost any model would match them.
#' * `"stats"`: summary statistics of the per-dyad rates in each simulated
#'   dataset (proportion of dyads with no events, mean, standard deviation and
#'   maximum), with the observed value shown as a vertical line. These are good at
#'   revealing too many zeroes or too much or too little variation between dyads.
#' * `"point"`: posterior medians and intervals of the edge weights plotted against
#'   simple point estimates (events / effort). This isn't a predictive check; it
#'   shows how BISoN estimates are pulled towards the prior when there are few
#'   data.
#'
#' No check can show that a model is correct, but these can reveal common ways
#' that it is wrong.
#'
#' @examples
#' sim_data <- simulate_bison_model("binary", aggregated = TRUE)
#' fit_edge <- bison_model(
#'   (event | duration) ~ dyad(node_1_id, node_2_id),
#'   data = sim_data$df_sim,
#'   model_type = "binary_conjugate",
#'   priors = get_default_priors("binary_conjugate")
#' )
#' plot_predictions(fit_edge, type = "density")
#' plot_predictions(fit_edge, type = "stats")
#' @export
plot_predictions <- function(obj, ...) {
  UseMethod("plot_predictions")
}

#' @rdname plot_predictions
#' @param num_draws Number of simulated datasets to compare against.
#' @param type One or more of `"density"`, `"stats"` and `"point"`.
#' @param draw_data Whether to show the observed data (set to `FALSE` for prior
#'   predictive checks).
#' @export
plot_predictions.bison_model <- function(obj, num_draws = 20, type = c("density", "point"),
                                         draw_data = TRUE, ...) {
  type <- match.arg(type, c("density", "stats", "point"), several.ok = TRUE)
  if ("point" %in% type && obj$num_dyads == 0) {
    type <- setdiff(type, "point")
  }

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))

  single_panels <- intersect(type, c("density", "point"))
  if ("stats" %in% type) {
    graphics::par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
  } else {
    graphics::par(mfrow = c(1, max(1, length(single_panels))), mar = c(4, 4, 2, 1))
  }

  # Group rows by dyad (or treat each row separately if there are no dyads).
  groups <- if (obj$num_dyads > 0) obj$model_info$row_dyad_ids else seq_len(obj$model_data$num_rows)
  divisor_by_group <- as.numeric(rowsum(as.numeric(obj$model_data$divisor), groups))
  keep <- divisor_by_group > 0
  group_rates <- function(events) {
    (as.numeric(rowsum(as.numeric(events), groups)) / divisor_by_group)[keep]
  }
  observed_rates <- group_rates(obj$model_data$event)
  draw_ids <- sample_draw_ids(nrow(obj$event_preds), num_draws)
  predicted_rates <- lapply(draw_ids, function(i) group_rates(obj$event_preds[i, ]))
  rate_label <- if (obj$model_type == "binary") "Proportion of sampling periods with an event" else "Event rate per unit effort"

  if ("density" %in% type) {
    plot_rate_densities(observed_rates, predicted_rates, draw_data, rate_label, bounded = obj$model_type == "binary")
  }
  if ("point" %in% type) {
    plot_point_comparison(obj, groups, keep)
  }
  if ("stats" %in% type) {
    plot_predictive_stats(observed_rates, predicted_rates, draw_data)
  }
  invisible(NULL)
}

plot_rate_densities <- function(observed, predicted, draw_data, xlab, bounded) {
  all_rates <- c(unlist(predicted), if (draw_data) observed)
  lower <- 0
  upper <- if (bounded) 1 else max(all_rates, 1e-6)
  # Use one bandwidth for every curve, so differences reflect the data and not the smoothing.
  bw <- stats::bw.nrd0(if (draw_data && length(unique(observed)) > 1) observed else unlist(predicted))
  bw <- max(bw, (upper - lower) / 200)
  dens <- function(x) stats::density(x, bw = bw, from = lower, to = upper, cut = 0)
  pred_dens <- lapply(predicted, dens)
  ymax <- max(sapply(pred_dens, function(d) max(d$y)))
  obs_dens <- NULL
  if (draw_data) {
    obs_dens <- dens(observed)
    ymax <- max(ymax, obs_dens$y)
  }
  graphics::plot(
    NULL, xlim = c(lower, upper), ylim = c(0, ymax * 1.1),
    main = if (draw_data) "Observed vs predicted rates" else "Predicted rates",
    xlab = xlab, ylab = "Density"
  )
  for (d in pred_dens) {
    graphics::lines(d, col = col2rgba(bison_colors[1], 0.35))
  }
  if (draw_data) {
    graphics::lines(obs_dens, lwd = 2.5, col = bison_highlight)
    graphics::legend("topright", legend = c("observed", "predicted"),
                     col = c(bison_highlight, bison_colors[1]), lwd = c(2.5, 1), bty = "n")
  } else {
    graphics::legend("topright", legend = "predicted", col = bison_colors[1], lwd = 1, bty = "n")
  }
}

plot_predictive_stats <- function(observed, predicted, draw_data) {
  stats_fns <- list(
    "Proportion of dyads with no events" = function(x) mean(x == 0),
    "Mean rate" = mean,
    "SD of rates" = function(x) if (length(x) > 1) stats::sd(x) else 0,
    "Maximum rate" = max
  )
  for (name in names(stats_fns)) {
    fn <- stats_fns[[name]]
    pred <- sapply(predicted, fn)
    obs <- if (draw_data) fn(observed) else NULL
    xlim <- range(c(pred, obs))
    if (diff(xlim) == 0) xlim <- xlim + c(-0.5, 0.5) * max(abs(xlim[1]), 1e-3)
    graphics::hist(
      pred, breaks = 15, xlim = xlim, main = name, xlab = "",
      col = col2rgba(bison_colors[2], 0.6), border = "white"
    )
    if (draw_data) {
      graphics::abline(v = obs, col = bison_highlight, lwd = 2.5)
    }
  }
  graphics::mtext(
    if (draw_data) "Histograms: simulated datasets. Red line: observed data." else "Histograms: simulated datasets.",
    side = 1, outer = TRUE, line = -1, cex = 0.8
  )
}

plot_point_comparison <- function(obj, groups, keep) {
  events <- as.numeric(rowsum(as.numeric(obj$model_data$event), groups))
  divisor <- as.numeric(rowsum(as.numeric(obj$model_data$divisor), groups))
  point_estimate <- (events / divisor)[keep]
  edgelist <- get_edgelist(obj, ci = 0.9)[keep, ]
  edgelist_inner <- get_edgelist(obj, ci = 0.5)[keep, ]
  graphics::plot(
    point_estimate, edgelist[, 3],
    ylim = range(c(edgelist[, 4], edgelist[, 5])),
    main = "Point vs BISoN estimates",
    xlab = "Point estimate (events / effort)", ylab = "BISoN estimate", type = "n"
  )
  graphics::abline(a = 0, b = 1, lty = 2, col = "grey50")
  graphics::segments(point_estimate, edgelist[, 4], point_estimate, edgelist[, 5], col = "grey40")
  graphics::segments(point_estimate, edgelist_inner[, 4], point_estimate, edgelist_inner[, 5],
                     lwd = 4, col = bison_colors[1])
  graphics::points(point_estimate, edgelist[, 3], pch = 19, cex = 0.6)
}

#' Plot the network with uncertainty
#'
#' Draws a sociogram where the uncertainty in each edge weight is shown by two
#' overlapping lines: a solid dark line whose width is the *lower* bound of the
#' credible interval (how weak the edge could plausibly be), and a wider,
#' translucent line whose width is the *upper* bound (how strong it could be).
#' Edges with wide translucent bands are uncertain.
#'
#' @param obj A fitted `bison_model` object.
#' @param ci Width of the credible interval, e.g. `0.9` for a 90% interval.
#' @param lwd Scaling factor for edge widths.
#' @param threshold Hide edges whose weight (see `threshold_on`) is below this
#'   value. This only affects the plot. `NULL` shows every edge.
#' @param threshold_on Which estimate to compare with `threshold`: the posterior
#'   `"median"` (default) or the `"lower"` bound of the credible interval (a
#'   stricter filter that only keeps edges that are clearly above the threshold).
#' @param show_all_nodes If `TRUE` (default), nodes with no remaining edges are
#'   still drawn.
#' @param layout Optional layout matrix (as returned invisibly by a previous call),
#'   so several networks can be drawn with nodes in the same positions.
#' @param ... Further arguments passed to [igraph::plot.igraph()], such as
#'   `vertex.size`.
#'
#' @return The layout matrix, invisibly.
#'
#' @examples
#' sim_data <- simulate_bison_model("binary", aggregated = TRUE)
#' fit_edge <- bison_model(
#'   (event | duration) ~ dyad(node_1_id, node_2_id),
#'   data = sim_data$df_sim,
#'   model_type = "binary_conjugate",
#'   priors = get_default_priors("binary_conjugate")
#' )
#' plot_network(fit_edge, lwd = 5)
#' plot_network(fit_edge, lwd = 5, threshold = 0.3)
#' @export
plot_network <- function(obj, ci = 0.9, lwd = 2, threshold = NULL,
                         threshold_on = c("median", "lower"), show_all_nodes = TRUE,
                         layout = NULL, ...) {
  threshold_on <- match.arg(threshold_on)
  edgelist <- get_edgelist(obj, ci = ci, transform = TRUE)
  keep <- rep(TRUE, nrow(edgelist))
  if (!is.null(threshold)) {
    keep <- edgelist[, if (threshold_on == "median") 3 else 4] >= threshold
  }

  net <- bison_base_igraph(obj, dyad_rows = which(keep))
  if (!show_all_nodes) {
    net <- igraph::delete_vertices(net, which(igraph::degree(net) == 0))
  }
  lb <- edgelist[keep, 4]
  ub <- edgelist[keep, 5]

  # Probabilities are already on a 0-1 scale; rates are rescaled so the widest edge has width lwd.
  scale <- if (obj$model_type == "binary" || length(ub) == 0) 1 else 1 / max(ub)

  if (is.null(layout)) {
    layout <- igraph::layout_nicely(net)
  }
  curved <- if (obj$directed) 0.15 else 0
  common <- list(
    net, layout = layout, vertex.label.color = "white", vertex.color = bison_colors[1],
    vertex.frame.color = "white", vertex.label.family = "sans", edge.curved = curved
  )
  do.call(igraph::plot.igraph, c(common, list(
    edge.width = ub * scale * lwd, edge.color = col2rgba(bison_colors[1], 0.3), edge.arrow.size = 0
  ), list(...)))
  do.call(igraph::plot.igraph, c(common, list(
    edge.width = lb * scale * lwd, edge.color = "grey20",
    edge.arrow.size = if (obj$directed) 0.4 else 0, add = TRUE
  ), list(...)))

  lb_pct <- format_percent(0.5 * (1 - ci))
  ub_pct <- format_percent(1 - 0.5 * (1 - ci))
  graphics::title(
    sub = paste0("Dark: ", lb_pct, "% bound of edge weight. Light: ", ub_pct, "% bound."),
    cex.sub = 0.8
  )
  invisible(layout)
}
