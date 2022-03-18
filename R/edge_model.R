require(rstan)
require(igraph)
require(bridgesampling)
require(dplyr)

#' Title
#'
#' @param formula
#' @param data
#' @param data_type
#' @param directed
#' @param method
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
edge_model <- function(formula, data, data_type=c("binary", "count", "duration"), directed=FALSE, method=c("mcmc", "vb", "inla"), verbose=FALSE, mc_cores=1) {
  # If verbose, print out MCMC chains.
  if (verbose) {
    refresh <- 200
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

  # Set up model and model data depending on data type.
  if (data_type == "binary") {
    model_data <- prepare_data(formula, data, directed, node_to_idx, node_list, data_type)
    if (dim(model_data$Z)[2] > 0) {
      # Mixed effects model
      model <- stanmodels$binary_mixed
    } else {
      # Fixed effects model
      model <- stanmodels$binary_fixed
    }
  }
  if (data_type == "count") {
    model_data <- prepare_data(formula, data, directed, node_to_idx, node_list, data_type)
    if (dim(model_data$Z)[2] > 0) {
      # Mixed effects model
      model <- stanmodels$count_mixed
    } else {
      # Fixed effects model
      model <- stanmodels$count_fixed
    }
  }
  if (data_type == "duration") {
    model_data <- prepare_data(formula, list(obs=data, obs_agg=data_agg), directed, node_to_idx, node_list, data_type)
    if (dim(model_data$Z)[2] > 0) {
      # Mixed effects model
      model <- stanmodels$duration_mixed
    } else {
      # Fixed effects model
      model <- stanmodels$duration_fixed
    }
  }

  # Fit model
  if (method == "mcmc"){
    fit <- rstan::sampling(model, data=model_data, refresh=refresh, cores=mc_cores)
  }
  if (method == "vb") {
    fit <- rstan::vb(model, data=model_data, output_samples=2000, iter=1e5, tol_rel_obj=1e-3, importance_resampling=TRUE)
  }

  # Extract edge weights from fitted edge model.
  chain <- rstan::extract(fit)$beta_fixed
  colnames(chain) <- colnames(data$X)

  # Prepare output object.
  obj <- list()
  obj$chain <- chain
  obj$num_nodes <- model_data$num_nodes
  obj$dyad_mapping <- model_data$dyad_id
  obj$directed <- directed
  obj$fit <- fit
  obj$formula <- formula
  obj$sri <- model_data$sri
  obj$data_type <- data_type
  obj$model_data <- model_data
  obj$data = data
  obj$stan_model = model
  class(obj) <- "edge_model"
  obj
}

#' Print edge model object
#'
#' @param obj
#'
#' @export
print.edge_model <- function(obj) {
  cat(paste0(
    "### Fitted BISoN edge model ###",
    "\nData type: ", obj$data_type,
    "\nFormula: ", format(obj$formula),
    "\nNumber of nodes: ", obj$num_nodes
  ))
}

#' Title
#'
#' @param obj
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot.edge_model <- function(obj, ...) {
  rstan::traceplot(obj$fit, ...)
}

prepare_data <- function(formula, observations, directed, node_to_idx, node_list, data_type) {
  if (data_type == "duration") {
    observations_agg <- observations$obs_agg
    observations <- observations$obs
  }

  obs <- get_all_vars(formula, observations)

  # If nodes are in factors, convert to idx.
  obs[, 2] <- as.integer(obs[, 2])
  obs[, 3] <- as.integer(obs[, 3])

  # Get number of nodes
  n <- length(node_to_idx)

  dyad_id <- matrix(0, n, n)
  if (directed == FALSE) {
    dyad_id[upper.tri(dyad_id)] <- 1:(0.5 * n * (n - 1))
    dyad_id <- dyad_id + t(dyad_id)
  } else {
    dyad_id[upper.tri(dyad_id)] <- 1:(0.5 * n * (n - 1))
    dyad_id[lower.tri(dyad_id)] <- (0.5 * n * (n - 1) + 1):n
  }

  # Calculate SRI for future reference
  obs_temp <- data.frame(
    association=obs[, 1],
    id_1=obs[, 2],
    id_2=obs[, 3]
  )
  sri <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      if (i < j) {
        sri[i, j] <- mean(subset(obs_temp, (id_1 == i & id_2 == j) | (id_1 == j & id_2 == i))$association)
      }
    }
  }

  y <- obs[, 1]

  dyad_ids=as.factor(dyad_id[cbind(node_to_idx[obs[, 2]], node_to_idx[obs[, 3]])])

  term_levels <- levels(dyad_ids)

  X <- data.frame(empty_col = rep(0, nrow(obs)))
  Z <- data.frame(empty_col = rep(0, nrow(obs)))

  for (i in 1:length(term_levels)) {
    term_level <- term_levels[i]
    new_term_name <-  paste0("dyad_", term_level)
    X[, new_term_name] <- 1 * (dyad_ids == term_level)
  }

  K_fixed <- length(term_levels) # Number of fixed effects
  K_random <- 0 # Number of random effects
  G <- c() # Variable grouping for random effects

  # Get additional effects.
  if (length(labels(terms(formula))) > 1) {
    additional_effects <- labels(terms(formula))
    additional_effects <- additional_effects[2:length(additional_effects)]
    for (i in 1:length(additional_effects)) {
      if (any(grep("\\|", additional_effects[i]))) {
        # Random effect
        term_vars <- sapply(colnames(obs), function(x) grepl(x, additional_effects[i]))
        term_var <- term_vars[term_vars == TRUE][1]

        term_name <- names(term_var)

        var_group <- paste0("random_", term_name)
        term_levels <- levels(as.factor(obs[, term_name]))
        for (term_level in term_levels) {
          new_term_name <-  paste0("random_", term_name, term_level)
          Z[, new_term_name] <- 1 * (as.factor(obs[, term_name]) == term_level)
          G[length(G) + 1] <- var_group
        }

      } else {
        # Fixed effect
        term_vars <- sapply(colnames(obs), function(x) grepl(x, additional_effects[i]))
        term_var <- term_vars[term_vars == TRUE][1]

        term_name <- names(term_var)

        # If it's a factor, create a column for each level.
        if (is.factor(obs[, term_name])) {
          var_group <- paste0("fixed_", term_name)
          term_levels <- levels(obs[, term_name])
          for (term_level in term_levels) {
            new_term_name <-  paste0("fixed_", term_name, term_level)
            X[, new_term_name] <- 1 * (obs[, term_name] == term_level)
          }
        } else {
          # Otherwise, create a single column:
          new_term_name <- paste0(c("fixed_", term_name), collapse="")
          X[, new_term_name] <- obs[, term_name]
        }
      }
    }
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
    node_names=node_list,
    num_nodes=n,
    sri=sri
  )

  if (data_type == "duration") {
    node_1_name <- colnames(obs)[2]
    node_2_name <- colnames(obs)[3]
    observations_agg$dyad_id <- dyad_id[cbind(node_to_idx[dplyr::pull(observations_agg[, node_1_name])], node_to_idx[dplyr::pull(observations_agg[, node_2_name])])]
    observations_agg <- observations_agg[order(observations_agg$dyad_id), ]
    print(observations_agg)
    data$k <- observations_agg$k
    data$d <- observations_agg$d
    data$observations_agg <- observations_agg
  }

  data
}
