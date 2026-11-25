#' Fit an edge mixture model
#'
#' Looks for clusters (components) in the edge weights, which can be interpreted
#' as different types of social relationship. A Gaussian mixture model is fitted
#' to the edge weights (on the link scale) in each of `num_draws` posterior
#' draws, for every number of components from 1 to `num_components`. This gives
#' the posterior probability of each number of components, and of each edge
#' belonging to each component.
#'
#' Components are labelled in increasing order of their mean, so component 1 is
#' always the weakest type of edge.
#'
#' @param edgemodel A fitted `bison_model` object.
#' @param num_components Maximum number of components to consider.
#' @param criterion `"ICL"` or `"BIC"`: the information criterion used to weigh
#'   the different numbers of components.
#' @param verbose Whether to show a progress bar.
#' @param num_draws Number of posterior draws to fit mixtures to (at most the
#'   number available).
#'
#' @return A `bison_mixture` object.
#' @seealso [get_network_component_probabilities()],
#'   [get_edge_component_probabilities()], [get_component_means()]
#'
#' @examples
#' \donttest{
#' sim_data <- simulate_bison_model_mixture(
#'   "binary", num_components = 2, component_weights = c(0.5, 0.5)
#' )
#' fit_edge <- bison_model(
#'   (event | duration) ~ dyad(node_1_id, node_2_id),
#'   data = sim_data$df_sim,
#'   model_type = "binary_conjugate",
#'   priors = get_default_priors("binary_conjugate")
#' )
#' fit_mixture <- bison_mixture(fit_edge, num_components = 3, verbose = FALSE)
#' summary(fit_mixture)
#' }
#' @export
bison_mixture <- function(edgemodel, num_components = 5, criterion = c("ICL", "BIC"),
                          verbose = TRUE, num_draws = 200) {
  criterion <- match.arg(criterion)
  if (is.null(edgemodel$edge_samples)) {
    stop("This model has no edge weights (it has no dyad() term).", call. = FALSE)
  }
  num_edges <- ncol(edgemodel$edge_samples)
  num_samples <- min(num_draws, nrow(edgemodel$edge_samples))
  draw_ids <- sample.int(nrow(edgemodel$edge_samples), num_samples)
  component_range <- seq_len(num_components)

  component_probability_samples <- matrix(0, num_samples, num_components)
  edge_component_samples <- array(NA_integer_, c(num_samples, num_components, num_edges))
  component_mean_samples <- lapply(component_range, function(k) matrix(NA_real_, num_samples, k))

  if (verbose) {
    message("Fitting ", num_components, " mixture model(s) to edge weights...")
    pb <- utils::txtProgressBar(max = num_samples, style = 3)
  }

  for (i in seq_len(num_samples)) {
    if (verbose) utils::setTxtProgressBar(pb, i)
    edge_weights <- as.numeric(edgemodel$edge_samples[draw_ids[i], ])
    ics <- rep(-Inf, num_components)
    for (k in component_range) {
      fit <- mclust::Mclust(edge_weights, G = k, verbose = FALSE)
      if (is.null(fit)) {
        next
      }
      # Relabel components in order of their means so labels agree across draws.
      means <- as.numeric(fit$parameters$mean)
      ord <- order(means)
      component_mean_samples[[k]][i, ] <- means[ord]
      edge_component_samples[i, k, ] <- match(fit$classification, ord)
      ics[k] <- fit[[tolower(criterion)]]
    }

    # Model weights from information criterion differences (mclust ICs are "larger is better").
    weights <- exp(0.5 * (ics - max(ics)))
    component_probability_samples[i, ] <- weights / sum(weights)
  }
  if (verbose) {
    close(pb)
  }

  component_probabilities <- colMeans(component_probability_samples)

  edge_component_probabilities <- lapply(component_range, function(k) {
    sapply(seq_len(k), function(component) {
      apply(edge_component_samples[, k, , drop = FALSE], 3, function(x) mean(x == component, na.rm = TRUE))
    })
  })
  edge_component_probabilities <- lapply(edge_component_probabilities, function(x) matrix(x, nrow = num_edges))

  obj <- list(
    component_probabilities = component_probabilities,
    edge_component_probabilities = edge_component_probabilities,
    component_probability_samples = component_probability_samples,
    edge_component_samples = edge_component_samples,
    edgemodel = edgemodel,
    component_range = component_range,
    num_components = num_components,
    component_mean_samples = component_mean_samples,
    criterion = criterion
  )
  class(obj) <- "bison_mixture"
  obj
}

#' Summary of a edge mixture model
#'
#' @param object An S3 edge mixture model
#' @param ... Additional arguments to be passed to summary calculations.
#'
#' @return Returns a summary object of an edge mixture model.
#' @examples
#' \donttest{
#' sim_data <- simulate_bison_model_mixture(
#'   "count", num_components = 2, component_weights = c(0.5, 0.5)
#' )
#' fit_edge <- bison_model(
#'   (event | duration) ~ dyad(node_1_id, node_2_id),
#'   data = sim_data$df_sim,
#'   model_type = "count_conjugate",
#'   priors = get_default_priors("count_conjugate")
#' )
#' fit_mixture <- bison_mixture(fit_edge, num_components = 3, verbose = FALSE)
#' summary(fit_mixture)
#' }
#' @export
summary.bison_mixture <- function(object, ...) {
  summary_obj <- list()
  summary_obj$description <- paste0(
    "=== Fitted dyadic mixture model ===\n",
    "Number of dyads: ", object$edgemodel$num_dyads, "\n",
    "Number of components: ", object$num_components, "\n"
  )
  summary_obj$component_probabilities <- object$component_probabilities
  names(summary_obj$component_probabilities) <- sapply(object$component_range, function(x) paste0("K = ", x))
  summary_obj$edge_component_probabilities <- object$edge_component_probabilities
  summary_obj$bison_mixture_obj <- object
  summary_obj$num_components <- object$num_components
  class(summary_obj) <- "summary.bison_mixture"
  summary_obj
}

#' Print information about an edge mixture model
#'
#' @param x An S3 edge mixture model summary.
#' @param digits Number of digits to round summary values to.
#' @param ... Additional arguments to be passed to print functions.
#'
#' @seealso [summary.bison_mixture()]
#' @keywords internal
#' @export
print.summary.bison_mixture <- function(x, digits=3, ...) {
  best_model <- which.max(x$component_probabilities)

  cat(paste0(
    "=== Fitted dyadic mixture model ===\n",
    "Maximum number of components: ", x$num_components, "\n",
    "Best model: ", best_model, " components\n",
    "Probability of best model: ", round(100 * x$component_probabilities[best_model], 1), "%\n"
  ))
  cat("=== Component probabilities ===\n")
  cp <- matrix(0, ncol=length(x$component_probabilities), nrow=1)
  rownames(cp) <- "P(K=k)"
  colnames(cp) <- 1:length(x$component_probabilities)
  cp[1, ] <- round(x$component_probabilities, digits)
  print(cp)

  component_means <- apply(
    mixture_transform(x$bison_mixture_obj, x$bison_mixture_obj$component_mean_samples[[best_model]]),
    2, stats::median, na.rm = TRUE
  )
  cat(paste0("=== Component means for best model (K = ", best_model, ") ===\n"))
  cm <- matrix(0, ncol=length(component_means), nrow=1)
  rownames(cm) <- "median"
  colnames(cm) <- 1:best_model
  cm[1, ] <- round(component_means, digits)
  print(cm)

  cat(paste0("=== Edge component probabilities for best model (K = ", best_model, ") ===\n"))
  ecp <- x$edge_component_probabilities[[best_model]]

  rownames(ecp) <- x$bison_mixture_obj$edgemodel$dyad_names
  colnames(ecp) <- 1:best_model

  print(head(round(ecp, digits)))
  if (nrow(ecp) > 6) cat("... (", nrow(ecp) - 6, " more rows; see get_edge_component_probabilities())\n", sep = "")
  invisible(x)
}

# Transform mixture component means from the link scale to the natural scale.
mixture_transform <- function(object, samples) {
  if (object$edgemodel$model_type == "binary") {
    return(stats::plogis(samples))
  }
  exp(samples)
}

#' Get probabilities of component membership for each edge
#'
#' @param object An S3 edge mixture model object
#' @param num_components The number of components in the mixture model
#'
#' @return Dataframe summarising edge component probabilities for each edge
#' @examples
#' \donttest{
#' sim_data <- simulate_bison_model_mixture(
#'   "count", num_components = 2, component_weights = c(0.5, 0.5)
#' )
#' fit_edge <- bison_model(
#'   (event | duration) ~ dyad(node_1_id, node_2_id),
#'   data = sim_data$df_sim,
#'   model_type = "count_conjugate",
#'   priors = get_default_priors("count_conjugate")
#' )
#' fit_mixture <- bison_mixture(fit_edge, num_components = 3, verbose = FALSE)
#' head(get_edge_component_probabilities(fit_mixture, 2))
#' }
#' @export
get_edge_component_probabilities <- function(object, num_components) {
  node_names <- names(object$edgemodel$node_to_idx)

  if (num_components <= length(object$edge_component_probabilities)) {
    df_names <- data.frame(
      node_1 = node_names[object$edgemodel$dyad_to_idx[, 1]],
      node_2 = node_names[object$edgemodel$dyad_to_idx[, 2]]
    )
    df_components <- data.frame(object$edge_component_probabilities[[num_components]])
    df <- cbind(df_names, df_components)
    colnames(df)[1:2] <- c("node_1", "node_2")
    colnames(df)[3:dim(df)[2]] <- sapply(1:num_components, function(i) paste0("P(K = ", i, ")"))
    return(df)
  }
  stop("This component does not exist.")
}

#' Get probabilities of component membership for the entire network
#'
#' @param object An S3 edge mixture model object
#'
#' @return A dataframe of probabilities of numbers of components over the entire network
#' @examples
#' \donttest{
#' sim_data <- simulate_bison_model_mixture(
#'   "count", num_components = 2, component_weights = c(0.5, 0.5)
#' )
#' fit_edge <- bison_model(
#'   (event | duration) ~ dyad(node_1_id, node_2_id),
#'   data = sim_data$df_sim,
#'   model_type = "count_conjugate",
#'   priors = get_default_priors("count_conjugate")
#' )
#' fit_mixture <- bison_mixture(fit_edge, num_components = 3, verbose = FALSE)
#' get_network_component_probabilities(fit_mixture)
#' }
#' @export
get_network_component_probabilities <- function(object) {
  num_components <- length(object$component_probabilities)
  df <- data.frame(num_components=1:num_components, probability=object$component_probabilities)
  return(df)
}

#' Get component means from edge mixture model
#'
#' @param object An S3 edge mixture model object
#' @param num_components The number of components in the mixture model
#' @param ci Credible interval for component mean estimates
#'
#' @return A dataframe describing the posteriors of the component means
#' @examples
#' \donttest{
#' sim_data <- simulate_bison_model_mixture(
#'   "count", num_components = 2, component_weights = c(0.5, 0.5)
#' )
#' fit_edge <- bison_model(
#'   (event | duration) ~ dyad(node_1_id, node_2_id),
#'   data = sim_data$df_sim,
#'   model_type = "count_conjugate",
#'   priors = get_default_priors("count_conjugate")
#' )
#' fit_mixture <- bison_mixture(fit_edge, num_components = 3, verbose = FALSE)
#' get_component_means(fit_mixture, 2)
#' }
#' @export
get_component_means <- function(object, num_components, ci=0.90) {
  lb_prob <- 0.5 * (1 - ci)
  ub_prob <- 1 - 0.5 * (1 - ci)
  component_mean_samples <- mixture_transform(object, object$component_mean_samples[[num_components]])
  mu <- apply(component_mean_samples, 2, stats::median, na.rm = TRUE)
  lb <- apply(component_mean_samples, 2, stats::quantile, probs = lb_prob, na.rm = TRUE)
  ub <- apply(component_mean_samples, 2, stats::quantile, probs = ub_prob, na.rm = TRUE)
  summary_table <- matrix(0, num_components, 3)
  rownames(summary_table) <- sapply(1:num_components, function(i) paste0("K = ", i))
  colnames(summary_table) <- c("50%", paste0(lb_prob * 100, "%"), paste0(ub_prob * 100, "%"))
  summary_table[, 1] <- mu
  summary_table[, 2] <- lb
  summary_table[, 3] <- ub
  return(as.data.frame(summary_table))
}
