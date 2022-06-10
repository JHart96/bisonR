prepare_edge_model_data <- function(formula, observations, directed, data_type) {
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
    dyad_to_idx[upper.tri(dyad_to_idx)] <- 1:(0.5 * num_nodes * (num_nodes - 1))
    dyad_to_idx <- dyad_to_idx + t(dyad_to_idx)
  } else {
    dyad_to_idx[upper.tri(dyad_to_idx)] <- 1:(0.5 * num_nodes * (num_nodes - 1))
    dyad_to_idx[lower.tri(dyad_to_idx)] <- (0.5 * num_nodes * (num_nodes - 1) + 1):num_nodes
  }

  # Get number of dyads
  num_dyads <- length(dyad_to_idx)

  # Get dyad names
  dyad_names <- sapply(
    1:num_dyads,
    function(x) paste0(names(node_to_idx)[which(dyad_to_idx == x, arr.ind=TRUE)[1, 2:1]], collapse=" <-> ")
  )

  # Get responses
  response <- dplyr::pull(observations, model_spec$event_var_name)

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
    dyad_ids=as.factor(dyad_id[cbind(node_to_idx[dplyr::pull(observations, model_spec$node_1_name)], node_to_idx[dplyr::pull(observations, model_spec$node_2_name)])])

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

  if (data_type %in% c("binary", ))
    # Check if divisor is a real
    if (!is.na(str_match(model_spec$divisor_var_name, "\\d+"))) {
      divisor <- rep(as.numeric(model_spec$divisor_var_name), nrow(observations))
    } else {
      divisor <- dplyr::pull(observations, model_spec$divisor_var_name)
    }

  # If divisor is a single value, convert it to a list.
  if (length(divisor) == 1) {
    divisor <- rep(divisor, nrow(observations))
  }

  data <- list(
    N=nrow(design_fixed),
    M=length(levels(dyad_ids)),
    y=y,
    design_fixed=data.matrix(design_fixed[, -1]),
    design_random=data.matrix(design_random[, -1]),
    K_fixed=ncol(design_fixed[, -1]),
    K_random=ncol(design_random[, -1]),
    R=length(unique(random_group_index)),
    random_group_index=as.integer(as.factor(random_group_index)),
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
    observations_agg$dyad_id <- dyad_id[cbind(node_to_idx[observations_agg[, node_1_name]], node_to_idx[observations_agg[, node_2_name]])]
    observations_agg <- observations_agg[order(observations_agg$dyad_id), ]
    data$k <- dplyr::pull(observations_agg, model_spec$divisor_var_name)
    data$d <- observations_agg$d
    data$observations_agg <- observations_agg
  }

  data
}
