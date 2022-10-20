#' Fit a dyadic mixture model
#'
#' @param edgemodel A fitted edge weight model.
#' @param num_components Maximum number of components to fit.
#' @param criterion "ICL" or "BIC" describing the criterion to use for model comparison.
#' @param verbose TRUE/FALSE indicating whether to output messages and progress when fitting.
#'
#' @return A bison_mixture object
#'
#' @export
bison_mixture <- function(edgemodel, num_components=5, criterion="ICL", verbose=TRUE) {
  # Prepare for fitting mixtures over the posterior networks
  num_samples <- dim(edgemodel$edge_samples)[1]
  num_edges <- dim(edgemodel$edge_samples)[2]
  component_range <- 1:num_components
  num_samples <- 200

  mclustBIC <- mclust::mclustBIC
  mclustICL <- mclust::mclustICL

  component_probability_samples <- matrix(0, num_samples, length(component_range))
  edge_component_samples <- array(0, c(num_samples, length(component_range), num_edges))

  if (verbose) {
    message(paste0("Fitting ", length(component_range), " mixture model(s) to edge weights..."))
    pb <- txtProgressBar(max=num_samples, style=3)
  }

  component_mean_samples <- list()
  # Fit mixture to each posterior draw of the networks
  for (i in 1:num_samples) {
    if (verbose) setTxtProgressBar(pb, i)
    # Fit a model for each number of components
    fit_mixtures <- list()
    for (k in component_range) {
      fit_mixtures[[k]] <- mclust::Mclust(edgemodel$edge_samples[i, ], G=k, verbose=FALSE)
      if (i == 1) {
        component_mean_samples[[k]] <- matrix(0, num_samples, k)
      }
      component_mean_samples[[k]][i, ] <- fit_mixtures[[k]]$parameters$mean
    }

    # Calculate model weights from IC differences
    ics <- sapply(fit_mixtures, function(x) x[[tolower(criterion)]])
    ics_diff <- max(ics) - ics
    half_exp_diffs <- exp(-0.5 * ics_diff)
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
  obj$fit_mixtures <- fit_mixtures
  obj$component_mean_samples <- component_mean_samples
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
  cp <- matrix(0, ncol=length(x$component_probabilities), nrow=1)
  rownames(cp) <- "P(K=k)"
  colnames(cp) <- 1:length(x$component_probabilities)
  cp[1, ] <- round(x$component_probabilities, digits)
  print(cp)

  component_means <- apply(x$bison_mixture_obj$component_mean_samples[[best_model]], 2, mean)
  cat(paste0("=== Component means for best model (K = ", best_model, ") ===\n"))
  cm <- matrix(0, ncol=length(component_means), nrow=1)
  rownames(cm) <- "mean"
  colnames(cm) <- 1:best_model
  component_means <- round(component_means, digits)
  if (x$bison_mixture_obj$edgemodel$model_type == "binary") {
    component_means <- plogis(component_means)
  } else if (x$bison_mixture_obj$edgemodel$model_type == "count") {
    component_means <- exp(component_means)
  }
  cm[1, ] <- component_means
  print(cm)

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

  if (num_components <= length(object$edge_component_probabilities)) {
    df_names <- data.frame(t(node_names))
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
#' @export
get_component_means <- function(object, num_components, ci=0.90) {
  lb_prob <- 0.5 * (1 - ci)
  ub_prob <- 1 - 0.5 * (1 - ci)
  component_mean_samples <- object$component_mean_samples[[num_components]]
  if (object$edgemodel$model_type == "binary") {
    component_mean_samples <- plogis(component_mean_samples)
  }
  if (object$edgemodel$model_type == "count") {
    component_mean_samples <- exp(component_mean_samples)
  }
  mu <- apply(component_mean_samples, 2, median)
  lb <- apply(component_mean_samples, 2, function(x) quantile(x, probs=c(lb_prob)))
  ub <- apply(component_mean_samples, 2, function(x) quantile(x, probs=c(ub_prob)))
  summary_table <- matrix(0, num_components, 3)
  rownames(summary_table) <- sapply(1:num_components, function(i) paste0("K = ", i))
  colnames(summary_table) <- c("50%", paste0(lb_prob * 100, "%"), paste0(ub_prob * 100, "%"))
  summary_table[, 1] <- mu
  summary_table[, 2] <- lb
  summary_table[, 3] <- ub
  return(as.data.frame(summary_table))
}
