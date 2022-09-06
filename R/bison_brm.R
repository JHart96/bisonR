#' Fit a brms model on BISoN network data
#'
#' @param formula A formula compatible with a brms model
#' @param edgemodel_list A bisonR edge model (or a list of them)
#' @param data_list A dataframe of regression variables compatible with a brms model (or a list of them)
#' @param num_draws Number of draws from the network posterior to use when fitting model
#' @param ...
#'
#' @return A fitted brms model
#' @export
bison_brm <- function(formula, edgemodel_list, data_list, num_draws=100, ...) {
  # Parse formula
  parsed_formula <- parse_bison_brms_formula(formula)
  brms_formula <- parsed_formula$brms_formula

  # Generate mice object
  mice_obj <- bison_mice(
    edgemodel_list,
    data_list,
    parsed_formula$param_names,
    parsed_formula$target_name,
    parsed_formula$metric_name,
    num_draws
  )

  # Run brms imputation.
  brms::brm_multiple(brms_formula, mice_obj, backend="cmdstanr", ...)
}

#' Generate a mice mids object
#'
#' @param edgemodel_list A bisonR edge model (or a list of them)
#' @param data_list A dataframe of regression variables compatible with a brms model (or a list of them)
#' @param num_draws Number of draws from the network posterior
#' @param param_names List of column names for network features (e.g. node_1_id, node_2_id or just node_id)
#' @param target_name Name of network feature (edge, node, global)
#' @param metric_name Name of network metric to calculate
#'
#' @return A mice mids object
#' @export
bison_mice <- function(edgemodel_list, data_list, param_names, target_name, metric_name, num_draws=100) {
  # If the data aren't lists yet, convert them
  if (class(edgemodel_list) != "list") edgemodel_list <- list(edgemodel_list)
  if (class(data_list) != "list") data_list <- list(data_list)

  new_bison_term <- paste0("bison_", target_name, "_", metric_name)

  posterior_samples_list <- list()

  for (i in 1:length(edgemodel_list)) {
    # Depending
    if (target_name == "edge") {
      node_1_ids <- sapply(
        dplyr::pull(data_list[[i]], param_names[1]),
        function(x) which(names(edgemodel_list[[i]]$node_to_idx) == x)
      )
      node_2_ids <- sapply(
        dplyr::pull(data_list[[i]], param_names[2]),
        function(x) which(names(edgemodel_list[[i]]$node_to_idx) == x)
      )
      dyad_ids <- get_dyad_ids(
        node_1_ids,
        node_2_ids,
        edgemodel_list[[i]]$dyad_to_idx,
        edgemodel_list[[i]]$directed
      )

      if (metric_name == "weight") {
        random_idxs <- sample(1:num_draws)
        posterior_samples <- extract_metric(edgemodel_list[[i]], paste0(target_name, "_", metric_name), num_draws)[random_idxs, dyad_ids]
      }
    } else if (target_name == "node") {
      node_ids <- sapply(
        dplyr::pull(data_list[[i]], param_names[1]),
        function(x) which(names(edgemodel_list[[i]]$node_to_idx) == x)
      )
      posterior_samples <- extract_metric(edgemodel_list[[i]], paste0(target_name, "_", metric_name), num_draws)[, node_ids]
    } else {
      stop("Network feature not supported.")
    }

    # Make sure posterior samples are in matrix format (for num_draws=1).
    posterior_samples_list[[i]] <- matrix(posterior_samples, nrow=num_draws)
  }

  # Generate list of dataframes.
  imputed_dataframes <- lapply(1:num_draws, function(i) {
    # For each posterior draw, combine dataframes from different edge models and data
    new_data_list <- lapply(1:length(edgemodel_list), function(j) {
      new_data <- data_list[[j]]
      new_data[new_bison_term] <- posterior_samples_list[[j]][i, ]
      new_data["bison_model"] <- j
      new_data
    })
    new_data <- dplyr::bind_rows(new_data_list)
    new_data$bison_model <- factor(new_data$bison_model)
    new_data
  })

  imputed_dataframes <- dplyr::bind_rows(imputed_dataframes, .id=".imp")
  imputed_dataframes <- dplyr::mutate(imputed_dataframes, .imp=as.integer(.imp))

  original_dataframes <- lapply(1:length(edgemodel_list), function(j) {
    new_data <- data_list[[j]]
    new_data["bison_model"] <- j
    new_data$bison_model <- factor(new_data$bison_model)
    new_data
  })
  original_dataframes <- bind_rows(original_dataframes)
  original_dataframes[".imp"] <- 0

  combined_dataframes <- bind_rows(list(original_dataframes, imputed_dataframes))

  mice::as.mids(as.data.frame(combined_dataframes))
}

parse_bison_brms_formula <- function(formula) {
  # Parse formula string
  formula_string <- deparse1(formula)[1]

  # Locate bison(...) term
  locs <- stringr::str_locate(formula_string, "bison\\([a-zA-Z0-9_,\\(\\)\\s]*\\)")
  formula_string_pre <- substr(formula_string, 1, locs[1] - 1)
  formula_string_post <- substr(formula_string, locs[2] + 1, nchar(formula_string))

  # Extract bison(...) term
  bison_term <- stringr::str_extract(formula_string, "bison\\([a-zA-Z0-9_,\\(\\)\\s]*\\)")
  bison_term <- stringr::str_replace(bison_term, " ", "")

  # Extract function name and arguments
  bison_term_split <- stringr::str_split(bison_term, "bison\\(|\\(|\\)")[[1]]
  bison_term_split <- bison_term_split[nchar(bison_term_split) > 0]
  metric_name <- bison_term_split[1]
  params <- bison_term_split[2]

  # Extract parameter names
  param_names <- stringr::str_split(params, ",")[[1]]
  target_name <- stringr::str_split(metric_name, "_")[[1]][1]
  metric_name <- stringr::str_split(metric_name, "_")[[1]][2]

  # Generate bison term for brms model.
  new_bison_term <- paste0("bison_", target_name, "_", metric_name)
  brms_formula <- paste0(formula_string_pre, new_bison_term, formula_string_post)

  return(list(
    param_names=param_names,
    target_name=target_name,
    metric_name=metric_name,
    new_bison_term=new_bison_term,
    brms_formula=brms_formula
  ))
}
