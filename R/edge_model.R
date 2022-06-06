require(igraph)
require(bridgesampling)
require(dplyr)
require(stringr)
require(bayesplot)

#' Fit an edge model to data
#'
#' @param formula
#' @param data
#' @param data_type
#' @param directed
#' @param method
#' @param verbose
#'
#' @return
#'
#' @export
edge_model <- function(formula, data, data_type=c("binary", "count", "duration"), directed=FALSE, method="mcmc", verbose=FALSE, mc_cores=4) {
  # If verbose, print out MCMC chains.
  if (verbose) {
    refresh <- 500
  } else {
    refresh <- 0
  }

  # Set fitting method
  if (length(method) > 1) {
    method = "mcmc"
  }


  # If using duration data, require data to be in list with aggregated data as well as raw observations.
  if (data_type == "duration") {
    data_agg <- data$obs_agg
    data <- data$obs
  }

  # Extract nodes from data and assign indices to them.
  obs <- get_all_vars(formula, data)
  if (is.factor(obs[, 2])) {
    node_list <- levels(obs[, 2])
  } else {
    node_list <- sort(unique(c(obs[, 2], obs[, 3])))
  }
  node_to_idx <- sapply(node_list, function(x) which(node_list == x))

  model_data <- prepare_data(formula, data, directed, node_to_idx, node_list, data_type)
  dyad_mapping <- model_data$dyad_id
  num_dyads <- max(dyad_mapping)
  dyad_names <- sapply(
    1:num_dyads,
    function(x) paste0(names(node_to_idx)[which(dyad_mapping == x, arr.ind=TRUE)[1, 2:1]], collapse=" <-> ")
  )

  # Set up model data depending on data type.
  if (data_type %in% c("binary", "count")) {
    model_data <- prepare_data(formula, data, directed, node_to_idx, node_list, data_type)
  } else if (data_type == "duration") {
    model_data <- prepare_data(formula, list(obs=data, obs_agg=data_agg), directed, node_to_idx, node_list, data_type)
  }

  # Build model
  model <- build_stan_model(data_type)

  # Fit model
  fit <- model$sample(data=model_data, refresh=refresh, chains=4, parallel_chains=mc_cores)

  # Extract edge weights from fitted edge model.
  chain <- fit$draws("beta_fixed", format="matrix")
  colnames(chain) <- colnames(data$X)

  # Prepare output object.
  obj <- list()
  obj$chain <- chain
  obj$num_nodes <- model_data$num_nodes
  obj$num_dyads <- num_dyads
  obj$dyad_mapping <- dyad_mapping
  obj$dyad_names <- dyad_names
  obj$directed <- directed
  obj$fit <- fit
  obj$formula <- formula
  # obj$sri <- model_data$sri
  obj$data_type <- data_type
  obj$model_data <- model_data
  obj$data = data
  obj$node_to_idx <- node_to_idx
  obj$stan_model = model
  obj$fit_method = method
  class(obj) <- "edge_model"
  obj
}

#' Summarises a fitted edge model object
#'
#' @param obj
#'
#' @return
#' @export
#'
#' @examples
summary.edge_model <- function(obj) {
  print(obj)
}

#' Prints out details of a fitted edge model object
#'
#' @param obj
#' @param ci
#'
#' @export
print.edge_model <- function(obj, ci=0.90) {
  cat(paste0(
    "=== Fitted BISoN edge model ===",
    "\nData type: ", obj$data_type,
    "\nFormula: ", format(obj$formula),
    "\nNumber of nodes: ", obj$num_nodes,
    "\nNumber of dyads: ", obj$num_dyads,
    "\nDirected: ", obj$directed,
    "\n=== Edge list summary ===\n"
  ))

  edgelist <- get_edgelist(obj, ci=ci)
  dyad_names <- do.call(paste, c(get_edgelist(fit_edge)[, 1:2], sep=" <-> "))
  summary_matrix <- as.matrix(edgelist[, 3:5])
  rownames(summary_matrix) <- dyad_names
  print(summary_matrix)
}

#' Retrieves an edgelist with uncertainty for a fitted edge weight model object
#'
#' @param obj
#'
#' @return
#' @export
#'
#' @examples
get_edgelist <- function (obj, ci=0.9) {
  node_names <- sapply(
    1:max(obj$dyad_mapping),
    function(x) names(obj$node_to_idx)[which(obj$dyad_mapping == x, arr.ind=TRUE)[1, 2:1]]
  )
  lb <- 0.5 * (1 - ci)
  ub <- 1 - lb
  dyad_start <- 1
  model_spec <- get_edge_model_spec(obj$formula)
  if(model_spec$intercept) {
    dyad_start <- 2
  }
  dyad_samples <- obj$chain[, dyad_start:(obj$num_dyads + dyad_start - 1)]
  edge_lower <- apply(dyad_samples, 2, function(x) quantile(x, probs=lb))
  edge_upper <- apply(dyad_samples, 2, function(x) quantile(x, probs=ub))
  edge_median <- apply(dyad_samples, 2, function(x) quantile(x, probs=0.5))
  edgelist <- data.frame(
    node_1 = node_names[1, ],
    node_2 = node_names[2, ],
    median = round(edge_median, 3),
    lb = round(edge_lower, 3),
    ub = round(edge_upper, 3)
  )
  colnames(edgelist)[4] <- paste0(as.character(lb * 100), "%")
  colnames(edgelist)[5] <- paste0(as.character(ub * 100), "%")
  edgelist
}

#' Draw samples from edgelist posterior
#'
#' @param obj
#' @param num_draws
#'
#' @return
#' @export
#'
#' @examples
draw_edgelist_samples <- function (obj, num_draws) {
  node_names <- sapply(
    1:max(obj$dyad_mapping),
    function(x) names(obj$node_to_idx)[which(obj$dyad_mapping == x, arr.ind=TRUE)[1, 2:1]]
  )
  edgelist_samples <- data.frame(
    node_1 = factor(node_names[1, ], levels=names(obj$node_to_idx)),
    node_2 = factor(node_names[2, ], levels=names(obj$node_to_idx)),
    draw = t(obj$chain[sample(1:dim(obj$chain)[1], size=num_draws, replace=TRUE), ])
  )
  edgelist_samples
}

plot_predictions.edge_model <- function(obj, num_draws=20) {
  y_preds <- obj$fit$draws("y_pred", format="matrix")
  df_draw <- data.frame(y=obj$model_data$y, dyad_ids=obj$model_data$dyad_ids)
  df_summed <- aggregate(y ~ as.factor(dyad_ids), df_draw, sum)
  pred_density <- density(df_summed$y)
  plot(pred_density, main="Observed vs predicted social events", xlab="Social events", ylim=c(0, max(pred_density$y) * 1.1))
  for (i in 1:num_draws) {
    df_draw$y <- as.vector(y_preds[i, ])
    df_summed <- aggregate(y ~ as.factor(dyad_ids), df_draw, sum)
    lines(density(df_summed$y), col=rgb(0, 0, 1, 0.5))
  }
}

plot_trace.edge_model <- function(obj, par_ids=1:12, ...) {
  if (obj$fit_method %in% c("mcmc")) {
    if (dim(obj$chain)[2] < 12) {
      par_ids <- 1:dim(obj$chain[2])
    }
    bayesplot::mcmc_trace(obj$fit$draws("beta_fixed")[,,par_ids])
  } else {
    message("plot_trace function not applicable to this fitting method")
  }
}

#' Sociogram plot with uncertainty of a fitted edge weight model object
#'
#' @param obj
#'
#' @return
#' @export
#'
#' @examples
plot_network <- function(obj, ci=0.9, lwd=1, ciwd=10) {
  edgelist <- get_edgelist(obj, ci=ci)
  net <- igraph::graph_from_edgelist(as.matrix(edgelist[, 1:2]), directed=FALSE)
  edge_weights <- plogis(edgelist[, 3])
  weights <- (edge_weights - min(edge_weights))/(max(edge_weights) - min(edge_weights))
  coords <- igraph::layout_nicely(net)
  igraph::plot.igraph(net, edge.width=weights * lwd, layout=coords)
  ci_widths <- plogis(edgelist[, 5]) - plogis(edgelist[, 4])
  weights <- (ci_widths - min(edge_weights))/(max(edge_weights) - min(edge_weights))
  igraph::plot.igraph(net, edge.width=weights * lwd * ciwd, layout=coords, edge.color=rgb(0, 0, 0, 0.25), add=TRUE)
}

prepare_data <- function(formula, observations, directed, node_to_idx, node_list, data_type) {
  if (data_type == "duration") {
    observations_agg <- observations$obs_agg
    observations <- observations$obs
  }

  model_spec <- get_edge_model_spec(formula)

  # Get number of nodes
  n <- length(node_to_idx)

  # Build matrix of dyad IDs
  dyad_id <- matrix(0, n, n)
  if (directed == FALSE) {
    dyad_id[upper.tri(dyad_id)] <- 1:(0.5 * n * (n - 1))
    dyad_id <- dyad_id + t(dyad_id)
  } else {
    dyad_id[upper.tri(dyad_id)] <- 1:(0.5 * n * (n - 1))
    dyad_id[lower.tri(dyad_id)] <- (0.5 * n * (n - 1) + 1):n
  }

  # Get responses
  y <- observations[, model_spec$event_var_name]

  # Create empty design matrices for fixed (X) and random (Z) effects
  X <- data.frame(empty_col = rep(0, nrow(observations)))
  Z <- data.frame(empty_col = rep(0, nrow(observations)))

  # Add intercept term if needed
  if (model_spec$intercept) {
    X[, "intercept"] <- 1
  }

  # If using dyad-level weights, populate design matrix
  if (!is.null(model_spec$node_1_name)) {
    # Get dyad IDs in the correct order
    dyad_ids=as.factor(dyad_id[cbind(node_to_idx[observations[, model_spec$node_1_name]], node_to_idx[observations[, model_spec$node_2_name]])])

    # Populate design matrix
    term_levels <- levels(dyad_ids)
    for (i in 1:length(term_levels)) {
      term_level <- term_levels[i]
      new_term_name <-  paste0("dyad_", term_level)
      X[, new_term_name] <- 1 * (dyad_ids == term_level)
    }
  }

  # Variable grouping for random effects
  G <- c()

  # Get additional fixed effects
  if (!is.null(model_spec$fixed)) {
    for (term_name in model_spec$fixed) {
      # If it's a factor, create a column for each level.
      if (is.factor(observations[, term_name])) {
        var_group <- paste0("fixed_", term_name)
        term_levels <- levels(observations[, term_name])
        for (term_level in term_levels) {
          new_term_name <-  paste0("fixed_", term_name, term_level)
          X[, new_term_name] <- 1 * (observations[, term_name] == term_level)
        }
      } else {
        # Otherwise, create a single column:
        new_term_name <- paste0(c("fixed_", term_name), collapse="")
        X[, new_term_name] <- observations[, term_name]
      }
    }
  }

  # Get additional random effects
  if (!is.null(model_spec$random)) {
    for (term_name in model_spec$random) {
      var_group <- paste0("random_", term_name)
      term_levels <- levels(as.factor(observations[, term_name]))
      for (term_level in term_levels) {
        new_term_name <-  paste0("random_", term_name, term_level)
        Z[, new_term_name] <- 1 * (as.factor(observations[, term_name]) == term_level)
        G[length(G) + 1] <- var_group
      }
    }
  }

  # Check if divisor is a real
  if (!is.na(str_match(model_spec$divisor_var_name, "\\d+"))) {
    divisor <- rep(as.numeric(model_spec$divisor_var_name), nrow(observations))
  } else {
    divisor <- observations[, model_spec$divisor_var_name]
  }

  # If divisor is a single value, convert it to a list.
  if (length(divisor) == 1) {
    divisor <- rep(divisor, nrow(obs))
  }

  data <- list(
    N=nrow(X),
    M=length(levels(dyad_ids)),
    y=y,
    X=data.matrix(X[, -1]),
    Z=data.matrix(Z[, -1]),
    K_fixed=ncol(X[, -1]),
    K_random=ncol(Z[, -1]),
    R=length(unique(G)),
    G=as.integer(as.factor(G)),
    node_to_idx=node_to_idx,
    dyad_id=dyad_id,
    dyad_ids=as.integer(dyad_ids),
    num_nodes=n,
    divisor=divisor
    # sri=sri
  )

  if (data_type == "duration") {
    node_1_name <- colnames(observations)[2]
    node_2_name <- colnames(observations)[3]
    observations_agg$dyad_id <- dyad_id[cbind(node_to_idx[dplyr::pull(observations_agg[, node_1_name])], node_to_idx[dplyr::pull(observations_agg[, node_2_name])])]
    observations_agg <- observations_agg[order(observations_agg$dyad_id), ]
    data$k <- observations_agg$k
    data$d <- observations_agg$d
    data$observations_agg <- observations_agg
  }

  data
}

get_edge_model_spec <- function(formula) {
  model_spec <- list()

  x <- str_split(deparse1(formula), "~")[[1]]
  lhs <- x[1]
  rhs <- x[2]

  # Process left hand side
  lhs_split <- str_split(lhs, "\\|")[[1]]
  event_var_name <- lhs_split[1]
  event_var_name <- str_replace_all(event_var_name, "\\(", "")
  event_var_name <- str_replace_all(event_var_name, " ", "")
  model_spec$event_var_name <- event_var_name

  divisor_var_name <- lhs_split[2]
  divisor_var_name <- str_replace_all(divisor_var_name, "\\)", "")
  divisor_var_name <- str_replace_all(divisor_var_name, " ", "")
  model_spec$divisor_var_name <- divisor_var_name

  # Set intercept to false by default
  model_spec$intercept <- FALSE

  model_spec$fixed <- c()
  model_spec$random <- c()

  rhs_split <- str_split(rhs, "\\+")[[1]]
  for (term in rhs_split) {
    term <- str_replace_all(term, " ", "")
    # Is it an intercept, a dyad, a fixed effect, or a random effect?
    if (!is.na(str_match(term, "^0|1$"))) {
      # Intercept (or lack thereof)
      if (term == "0") {
        model_spec$intercept <- FALSE
      } else {
        model_spec$intercept <- TRUE
      }
    } else if (!is.na(str_match(term, "^dyad\\(.*,.*\\)$")[[1]])) {
      # dyad(,) term
      node_names <- str_split(term, "\\(|\\)")[[1]][2]
      node_names <- str_replace_all(node_names, " ", "")
      node_names_split <- str_split(node_names, ",")[[1]]
      model_spec$node_1_name <- node_names_split[1]
      model_spec$node_2_name <- node_names_split[2]
    } else if (is.na(str_match(term, "[^a-zA-Z0-9]"))) {
      # No non-alphanumeric characters, and it can't be an intercept, so it's a fixed effect
      model_spec$fixed[length(model_spec$fixed) + 1] <- term
    } else if (!is.na(str_match(term, "^\\(1\\|.*\\)$"))) {
      # Contains a (1 | *) structure, so it's a basic random effect
      term_name <- str_split(term, "\\(|\\||\\)")[[1]][3]
      model_spec$random[length(model_spec$random) + 1] <- term_name
    } else {
      warning(paste0("Formula term \"", term, "\" not supported by bisonR. Check the formula is correctly specified."))
    }
  }

  model_spec
}

# Code for INLA and conjugate fitting. Not in use for now.
# prepare_data_inla <- function (formula, observations, directed, node_to_idx, node_list, data_type) {
#   obs <- get_all_vars(formula, observations)
#
#   # If nodes are in factors, convert to idx.
#   obs[, 2] <- as.integer(obs[, 2])
#   obs[, 3] <- as.integer(obs[, 3])
#
#   # Get number of nodes
#   n <- length(node_to_idx)
#
#   dyad_id <- matrix(0, n, n)
#   if (directed == FALSE) {
#     dyad_id[upper.tri(dyad_id)] <- 1:(0.5 * n * (n - 1))
#     dyad_id <- dyad_id + t(dyad_id)
#   } else {
#     dyad_id[upper.tri(dyad_id)] <- 1:(0.5 * n * (n - 1))
#     dyad_id[lower.tri(dyad_id)] <- (0.5 * n * (n - 1) + 1):n
#   }
#
#   # Calculate SRI for future reference
#   obs_temp <- data.frame(
#     association=obs[, 1],
#     id_1=obs[, 2],
#     id_2=obs[, 3]
#   )
#   sri <- matrix(0, n, n)
#   for (i in 1:n) {
#     for (j in 1:n) {
#       if (i < j) {
#         sri[i, j] <- mean(subset(obs_temp, (id_1 == i & id_2 == j) | (id_1 == j & id_2 == i))$association)
#       }
#     }
#   }
#
#   obs <- get_all_vars(formula, df)
#   obs$dyad_id <- fit_edge$dyad_mapping[cbind(obs[, 2], obs[, 3])]
#   obs$dyad_id <- factor(obs$dyad_id)
#
#
#   term_labels <- c('0', 'dyad_id')
#
#   if (length(labels(terms(formula))) > 1) {
#     additional_effects <- labels(terms(formula))
#     additional_effects <- additional_effects[2:length(additional_effects)]
#     for (i in 1:length(additional_effects)) {
#       if (any(grep("\\|", additional_effects[i]))) {
#         # Random effect
#         term_vars <- sapply(colnames(obs), function(x) grepl(x, additional_effects[i]))
#         term_var <- term_vars[term_vars == TRUE][1]
#
#         term_name <- names(term_var)
#         term_labels <- c(term_labels, paste0("f(", term_name, ", model='iid', hyper=prior.random)"))
#
#       } else {
#         # Fixed effect
#         term_vars <- sapply(colnames(obs), function(x) grepl(x, additional_effects[i]))
#         term_var <- term_vars[term_vars == TRUE][1]
#
#         term_name <- names(term_var)
#         term_labels <- c(term_labels, term_name)
#       }
#     }
#   }
#
#   formula_inla <- reformulate(termlabels = term_labels, response = all.vars(formula)[[1]])
#   data <- list(
#     formula = formula_inla,
#     df = obs
#   )
# }
# EDGE MODEL ADDITION
# else if (method == "inla") {
#   model_data_inla <- prepare_data_inla(formula, data, directed, node_to_idx, node_list, data_type)
#   prior.fixed <- list(mean=0, prec=1)
#   prior.random <- list(prec=list(prior="normal", param=c(0, 1)))
#   model = NULL
#   # Fit the INLA model
#   if (length(divisor) == 1) {
#     divisor <- rep(divisor, nrow(model_data_inla$df))
#   }
#   if (data_type == "binary") {
#     fit <- INLA::inla(model_data_inla$formula,
#                       family="binomial",
#                       data=model_data_inla$df,
#                       Ntrials=divisor,
#                       control.fixed=prior.fixed,
#                       control.compute=list(config = TRUE)
#     )
#   } else if (data_type == "count") {
#     model_data$df$divisor <- divisor
#     inla_formula <- reformulate(termlabels = c(labels(terms(formula)), "offset(log(divisor))"), response=all.vars(formula)[1])
#     fit <- INLA::inla(model_data$formula, # Add offset(log(duration)) to this and Stan model. Figure out syntax for formula.
#                       family="poisson",
#                       data=model_data$df,
#                       control.fixed=prior.fixed,
#                       control.compute=list(config = TRUE)
#     )
#   } else if (data_type == "duration") {
#     message("This model type is not supported by INLA")
#   }
#
#   # Extract edge weights from fitted edge model.
#   num_samples <- 1000
#   inla_samples <- inla.posterior.sample(num_samples, fit)
#   chain <- matrix(0, length(inla_samples), num_dyads)
#   for (i in 1:length(inla_samples)) {
#     chain[i, ] <- tail(inla_samples[[i]]$latent, num_dyads)
#   }
# } else if (method == "conjugate") {
#   fit <- NA
#   model <- NA
#   num_samples <- 1000
#   if (length(labels(terms(fit_edge$formula))) == 1) {
#     chain <- matrix(0, num_samples, num_dyads)
#     if (data_type == "binary") {
#       for (i in 1:num_dyads) {
#         chain[, i] <- qlogis(rbeta(num_samples, 1 + data[i, lhs_components[1]], 1 + divisor[i]))
#       }
#     } else if (data_type == "count") {
#       for (i in 1:num_dyads) {
#         chain[, i] <- log(rgamma(num_samples, 0.001 + data[i, lhs_components[1]], 0.001 + 1)/divisor[i])
#       }
#     }
#   } else {
#     message("Additional effects cannot be fitted using the conjugate method")
#   }
# }
