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
    if (obj$model_type %in% c("binary_conjugate", "count_conjugate")) {
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
  if (obj$model_type == "binary") {
    return(plogis(edge_samples))
  }
  if (obj$model_type == "count") {
    return(exp(edge_samples))
  }
}

#' Convert bisonR edgemodel posterior to igraph network objects
#'
#' @param obj bisonR edge model
#' @param num_draws Number of draws to sample from network posterior
#'
#' @return List of igraph objects representing posterior distribution of networks
#' @export
bison_to_igraph <- function(obj, num_draws=1000) {
  edgelist_samples <- draw_edgelist_samples(obj, num_draws)
  base_net <- igraph::graph_from_edgelist(as.matrix(edgelist_samples[, 1:2]), directed = obj$directed)
  nets <- lapply(1:num_draws, function(i) {
    if (obj$model_type == "binary") {
      igraph::E(base_net)$weight <- plogis(edgelist_samples[, 2 + i])
    }
    if (obj$model_type == "count") {
      igraph::E(base_net)$weight <- exp(edgelist_samples[, 2 + i])
    }
    base_net
  })
  return(nets)
}

#' Convert duration data to binary
#' @description Uses a naive method for converting duration data to binary data. Use with caution. Ensure that all time periods are in the same units.
#'
#' @param df Original dataframe, ready for use with `bison_model`.
#' @param event_name The variable name in the dataframe corresponding to event durations.
#' @param divisor_name The variable name in the dataframe corresponding to divisor durations (e.g. observation time).
#' @param max_event_length The maximum length of time a single social event could take place for.
#'
#' @return A dataframe the same as the input but with `event_name` and `divisor_name` columns replaced with values compatible with binary models.
#' @export
convert_duration_to_binary <- function(df, event_name, divisor_name, max_event_length) {
  event_data <- dplyr::pull(df, event_name)
  divisor_data <- dplyr::pull(df, divisor_name)
  event_samples <- round(event_data/max_event_length)
  divisor_samples <- round(divisor_data/max_event_length)
  new_df <- df
  new_df[event_name] <- event_samples
  new_df[divisor_name] <- divisor_samples
  new_df
}

#' Convert group-by-individual matrix to bisonR format
#'
#' @param gbi Group-by-individual (R x C) matrix corresponding to R groups of C individuals.
#' @param group_properties A vector of length R describing properties of the groups to be added to the dataframe.
#' @param individual_properties A vector of length C describing properties of the individuals to be added to the dataframe.
#' @param individual_constraints A vector of length C describing constraints of the individuals.
#'
#' @return A dataframe for use in bisonR, where each row corresponds to a possible social event between two individuals.
#' @export
convert_gbi_to_bison <- function(gbi, group_properties=NULL, individual_properties=NULL, individual_constraints=NULL) {
  # Define empty bisonR dataframe
  df = data.frame(node_1=numeric(), node_2=numeric(), event=numeric(), group_id=numeric())
  if (!is.null(group_properties)) {
    df$group_property = numeric()
  }
  if (!is.null(individual_properties)) {
    df$node_1_property = numeric()
    df$node_2_property = numeric()
  }

  # If there are no individual constraints, set the constraints vector to the same value.
  if (is.null(individual_constraints)) {
    individual_constraints <- rep(0, ncol(gbi))
  }

  # If there are no node names on the GBI, define them numerically.
  node_names <- colnames(gbi)
  if (is.null(names(gbi))) {
    node_names <- 1:ncol(gbi)
  }

  # For each group in the GBI matrix, calculate all possible dyadic events.
  for (row_idx in 1:nrow(gbi)) {
    for (i in which(gbi[row_idx, ] == 1)) {
      for (j in 1:ncol(gbi)) {
        if (i != j && individual_constraints[i] == individual_constraints[j]) {
          event = gbi[row_idx, i] * gbi[row_idx, j]
          new_row <- list(
            node_1=node_names[i],
            node_2=node_names[j],
            event=event,
            group_id=row_idx
          )
          if (!is.null(group_properties)) {
            new_row$group_property=group_properties[row_idx]
          }
          if (!is.null(individual_properties)) {
            new_row$node_1_property=individual_properties[i]
            new_row$node_2_property=individual_properties[j]
          }
          df[nrow(df) + 1, ] <- new_row
        }
      }
    }
  }

  df
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
