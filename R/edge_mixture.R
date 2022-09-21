#' Fit a dyadic mixture model
#'
#' @param edgemodel A fitted edge weight model.
#' @param num_components Maximum number of components to fit.
#' @param verbose TRUE/FALSE indicating whether to output messages and progress when fitting.
#'
#' @return A bison_mixture object
#'
#' @export
bison_mixture <- function(edgemodel, num_components=5, verbose=TRUE) {
  # Prepare for fitting mixtures over the posterior networks
  num_samples <- dim(edgemodel$edge_samples)[1]
  num_edges <- dim(edgemodel$edge_samples)[2]
  component_range <- 1:num_components
  num_samples <- 200

  mclustBIC <- mclust::mclustBIC

  component_probability_samples <- matrix(0, num_samples, length(component_range))
  edge_component_samples <- array(0, c(num_samples, length(component_range), num_edges))

  if (verbose) {
    message(paste0("Fitting ", length(component_range), " mixture model(s) to edge weights..."))
    pb <- txtProgressBar(max=num_samples, style=3)
  }

  # Fit mixture to each posterior draw of the networks
  for (i in 1:num_samples) {
    if (verbose) setTxtProgressBar(pb, i)
    # Fit a model for each number of components
    fit_mixtures <- list()
    for (k in component_range) {
      fit_mixtures[[k]] <- mclust::Mclust(edgemodel$edge_samples[i, ], G=k, verbose=FALSE)
    }

    # Calculate model weights from BIC differences
    bics <- sapply(fit_mixtures, function(x) x$BIC[1])
    bics_diff <- max(bics) - bics
    half_exp_diffs <- exp(-0.5 * bics_diff)
    component_probability_samples[i, ] <- half_exp_diffs/sum(half_exp_diffs)

    # Draw samples of component memberships per posterior draw
    edge_component_samples[i, , ] <- t(sapply(fit_mixtures, function(x) x$classification))
  }
  if (verbose) {
    close(pb)
  }

  # Calculate network-level component probabilities from samples
  component_probabilities <- apply(component_probability_samples, 2, mean)

  # Calculate edge-level component probabilities for each model from samples
  edge_component_probabilities <- lapply(component_range, function(k) {
    cbind(sapply(1:k, function(i) {apply(edge_component_samples[, k, ], 2, function(x) mean(x == i))}))
  })

  obj <- list()
  obj$component_probabilities <- component_probabilities
  obj$edge_component_probabilities <- edge_component_probabilities
  obj$component_probability_samples <- component_probability_samples
  obj$edge_component_samples <- edge_component_samples
  obj$edgemodel <- edgemodel
  obj$component_range <- component_range
  obj$num_components <- length(component_range)
  class(obj) <- "bison_mixture"
  obj
}

#' Summary of a edge mixture model
#'
#' @param object An S3 edge mixture model
#' @param ... Additional arguments to be passed to summary calculations.
#'
#' @return Returns a summary object of an edge mixture model.
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
  class(summary_obj) <- "summary.bison_mixture"
  summary_obj
}

#' Print information about an edge mixture model
#'
#' @param x An S3 edge mixture model summary.
#' @param digits Number of digits to round summary values to.
#' @param ... Additional arguments to be passed to print functions.
#'
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
  cp <- matrix(round(x$component_probabilities, digits), nrow=1)
  rownames(cp) <- "P(K=k)"
  colnames(cp) <- 1:length(x$component_probabilities)
  print(cp)
  cat(paste0("=== Edge component probabilities for best model (K = ", best_model, ") ===\n"))
  ecp <- x$edge_component_probabilities[[best_model]]

  rownames(ecp) <- do.call(paste, c(get_edgelist(x$bison_mixture_obj$edgemodel)[, 1:2], sep=" <-> "))
  colnames(ecp) <- 1:best_model

  print(head(ecp))
  cat("...")
}

#' Get probabilities of component membership for each edge
#'
#' @param object An S3 edge mixture model object
#' @param num_components The number of components in the mixture model
#'
#' @return Dataframe summarising edge component probabilities for each edge
#' @export
get_edge_component_probabilities <- function(object, num_components) {

  node_names <- sapply(
    1:nrow(object$edgemodel$dyad_to_idx),
    function(x) names(object$edgemodel$node_to_idx)[object$edgemodel$dyad_to_idx[x, ]]
  )

  if (num_components < length(object$edge_component_probabilities)) {
    df_names <- data.frame(t(node_names))
    df_components <- data.frame(object$edge_component_probabilities[[num_components]])
    df <- cbind(df_names, df_components)
    colnames(df)[1:2] <- c("node_1", "node_2")
    colnames(df)[3:dim(df)[2]] <- sapply(1:num_components, function(i) paste0("P(K = ", i, ")"))
    return(df)
  }
}

#' Get probabilities of component membership for the entire network
#'
#' @param object An S3 edge mixture model object
#'
#' @return A dataframe of probabilities of numbers of components over the entire network
#' @export
get_network_component_probabilities <- function(object) {
  num_components <- length(object$component_probabilities)
  df <- data.frame(num_components=1:num_components, probability=object$component_probabilities)
  return(df)
}
