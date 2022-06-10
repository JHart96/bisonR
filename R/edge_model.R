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
edge_model <- function(formula, data, data_type=c("binary", "count"), directed=FALSE, priors=NULL, refresh=0, mc_cores=4) {
  if (data_type == "duration") {
    stop("Duration model not yet supported")
  }

  # If user-specified priors haven't been set, use the defaults
  if (is.null(priors)) priors <- get_default_priors(data_type)

  # Set up model data depending on data type.
  model_info <- prepare_data(formula, data, directed, data_type)

  # Set the priors in model data
  prior_parameters <- extract_prior_parameters(priors)
  model_data <- c(model_info$model_data, prior_parameters)

  # Build model
  model <- build_stan_model(data_type)

  # Fit model
  fit <- model$sample(data=model_data, refresh=refresh, chains=4, parallel_chains=mc_cores, step_size=0.1)

  # Extract edge weights from fitted edge model.
  chain <- fit$draws("beta_fixed", format="matrix")
  colnames(chain) <- colnames(model_data$design_fixed)

  # Prepare output object.
  obj <- list(
    chain = chain,
    num_nodes = model_info$num_nodes,
    num_dyads = model_info$num_dyads,
    node_to_idx = model_info$node_to_idx,
    dyad_to_idx = model_info$dyad_to_idx,
    dyad_names = model_info$dyad_names,
    model_info = model_info,
    directed = directed,
    fit = fit,
    formula = formula,
    data_type = data_type,
    model_data = model_data,
    input_data = data,
    stan_model = model
  )
  class(obj) <- "edge_model"
  return(obj)
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
    1:max(obj$dyad_to_idx),
    function(x) names(obj$node_to_idx)[which(obj$dyad_to_idx == x, arr.ind=TRUE)[1, 2:1]]
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
  event_preds <- obj$fit$draws("event_pred", format="matrix")
  df_draw <- data.frame(event=obj$model_data$event, dyad_id=obj$model_info$row_dyad_ids)
  df_summed <- aggregate(event ~ as.factor(dyad_id), df_draw, sum)
  pred_density <- density(df_summed$event)
  plot(pred_density, main="Observed vs predicted social events", xlab="Social events", ylim=c(0, max(pred_density$y) * 1.1))
  for (i in 1:num_draws) {
    df_draw$event <- as.vector(event_preds[i, ])
    df_summed <- aggregate(event ~ as.factor(dyad_id), df_draw, sum)
    lines(density(df_summed$event), col=rgb(0, 0, 1, 0.5))
  }
}

plot_trace.edge_model <- function(obj, par_ids=1:12, ...) {
  if (dim(obj$chain)[2] < 12) {
    par_ids <- 1:dim(obj$chain[2])
  }
  bayesplot::mcmc_trace(obj$fit$draws("beta_fixed")[,,par_ids])
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

prepare_data <- function(formula, observations, directed, data_type) {
  if (data_type == "duration") {
    observations_agg <- observations[[2]]
    observations <- observations[[1]]
  }

  # Get model specification from formula
  model_spec <- get_edge_model_spec(formula)

  # Get node to index mapping
  node_to_idx <- 1:length(levels(dplyr::pull(observations, model_spec$node_1_name)))
  names(node_to_idx) <- levels(dplyr::pull(observations, model_spec$node_1_name))

  # Get number of nodes
  num_nodes <- length(node_to_idx)

  # Get dyad to index mapping
  dyad_to_idx <- matrix(0, num_nodes, num_nodes)
  if (directed == FALSE) {
    num_dyads <- (0.5 * num_nodes * (num_nodes - 1))
    dyad_to_idx[upper.tri(dyad_to_idx)] <- 1:(0.5 * num_nodes * (num_nodes - 1))
    dyad_to_idx <- dyad_to_idx + t(dyad_to_idx)
  } else {
    num_dyads <- num_nodes * (num_nodes - 1)
    dyad_to_idx[upper.tri(dyad_to_idx)] <- 1:(0.5 * num_nodes * (num_nodes - 1))
    dyad_to_idx[lower.tri(dyad_to_idx)] <- (0.5 * num_nodes * (num_nodes - 1) + 1):num_nodes
  }

  # Get dyad names
  dyad_names <- sapply(
    1:num_dyads,
    function(x) paste0(names(node_to_idx)[which(dyad_to_idx == x, arr.ind=TRUE)[1, 2:1]], collapse=" <-> ")
  )

  # Get events
  event <- dplyr::pull(observations, model_spec$event_var_name)

  # Create empty design matrices for fixed (design_fixed) and random (design_random) effects
  design_fixed <- data.frame(empty_col = rep(0, nrow(observations)))
  design_random <- data.frame(empty_col = rep(0, nrow(observations)))

  # Add intercept term if needed
  if (model_spec$intercept) {
    design_fixed[, "intercept"] <- 1
  }

  # If using dyad-level weights, populate design matrix
  if (!is.null(model_spec$node_1_name)) {
    # Get dyad IDs in the correct order
    node_1_names <- dplyr::pull(observations, model_spec$node_1_name)
    node_2_names <- dplyr::pull(observations, model_spec$node_2_name)
    dyad_ids=as.factor(dyad_to_idx[cbind(node_to_idx[node_1_names], node_to_idx[node_2_names])])

    # Populate design matrix
    term_levels <- levels(dyad_ids)
    for (i in 1:length(term_levels)) {
      term_level <- term_levels[i]
      new_term_name <-  paste0("dyad_", term_level)
      design_fixed[, new_term_name] <- 1 * (dyad_ids == term_level)
    }
  }

  # Variable grouping for random effects
  random_group_index <- c()

  # Get additional fixed effects
  if (!is.null(model_spec$fixed)) {
    for (term_name in model_spec$fixed) {
      # If it's a factor, create a column for each level.
      if (is.factor(dplyr::pull(observations, term_name))) {
        var_group <- paste0("fixed_", term_name)
        term_levels <- levels(dplyr::pull(observations, term_name))
        for (term_level in term_levels) {
          new_term_name <-  paste0("fixed_", term_name, term_level)
          design_fixed[, new_term_name] <- 1 * (dplyr::pull(observations, term_name) == term_level)
        }
      } else {
        # Otherwise, create a single column:
        new_term_name <- paste0(c("fixed_", term_name), collapse="")
        design_fixed[, new_term_name] <- dplyr::pull(observations, term_name)
      }
    }
  }

  # Get additional random effects
  if (!is.null(model_spec$random)) {
    for (term_name in model_spec$random) {
      var_group <- paste0("random_", term_name)
      term_levels <- levels(as.factor(dplyr::pull(observations, term_name)))
      for (term_level in term_levels) {
        new_term_name <-  paste0("random_", term_name, term_level)
        design_random[, new_term_name] <- 1 * (as.factor(dplyr::pull(observations, term_name)) == term_level)
        random_group_index[length(random_group_index) + 1] <- var_group
      }
    }
  }

  if (data_type %in% c("binary", "count")) {
    # Check if divisor is a real
    if (!is.na(str_match(model_spec$divisor_var_name, "\\d+"))) {
      divisor <- rep(as.numeric(model_spec$divisor_var_name), nrow(observations))
    } else {
      divisor <- dplyr::pull(observations, model_spec$divisor_var_name)
    }
  }

  # If divisor is a single value, convert it to a list.
  if (length(divisor) == 1) {
    divisor <- rep(divisor, nrow(observations))
  }

  model_data <- list(
    num_rows=nrow(design_fixed),
    event=event,
    divisor=divisor,
    design_fixed=data.matrix(design_fixed[, -1]),
    design_random=data.matrix(design_random[, -1]),
    num_fixed=ncol(design_fixed[, -1]),
    num_random=ncol(design_random[, -1]),
    num_random_groups=length(unique(random_group_index)),
    random_group_index=as.integer(as.factor(random_group_index))
  )

  # ,
  #   node_to_idx=node_to_idx,
  #   dyad_id=dyad_id,
  #   dyad_ids=as.integer(dyad_ids),
  #   num_nodes=n,
  #
  #   # sri=sri
  # )

  if (data_type == "duration") {
    model_data$k <- dplyr::pull(observations_agg, model_spec$divisor_var_name)
    # data$d <- observations_agg$d
    model_data$observations_agg <- observations_agg
  }

  obj <- list(
    model_data=model_data,
    node_to_idx=node_to_idx,
    dyad_to_idx=dyad_to_idx,
    dyad_names=dyad_names,
    num_nodes=num_nodes,
    num_dyads=num_dyads,
    row_dyad_ids=dyad_ids
  )

  return(obj)
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
      # Intercept
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
    } else if (is.na(str_match(term, "[^a-zA-Z0-9_]"))) {
      # Fixed effect
      model_spec$fixed[length(model_spec$fixed) + 1] <- term
    } else if (!is.na(str_match(term, "^\\(1\\|.*\\)$"))) {
      # Random intercept
      term_name <- str_split(term, "\\(|\\||\\)")[[1]][3]
      model_spec$random[length(model_spec$random) + 1] <- term_name
    } else {
      warning(paste0("Formula term \"", term, "\" not supported by bisonR. Check the formula is correctly specified."))
    }
  }

  model_spec
}
