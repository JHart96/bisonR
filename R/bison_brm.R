#' Fit a brms model on BISoN network data
#'
#' @param formula A formula compatible with a brms model
#' @param edgemodel A fitted bisonR edge model
#' @param data A dataframe of regression variables compatible with a brms model
#' @param family Regression family compatible with a brms model
#' @param num_draws Number of draws from the network posterior to use when fitting model
#' @param ...
#'
#' @return A fitted brms model
#' @export
bison_brm <- function(formula, edgemodel, data, num_draws=100, ...) {
  # Parse formula
  parsed_formula <- parse_bison_brms_formula(formula)
  param_names <- parsed_formula$param_names
  target_name <- parsed_formula$target_name
  metric_name <- parsed_formula$metric_name
  new_bison_term <- parsed_formula$new_bison_term
  brms_formula <- parsed_formula$brms_formula

  # Depending
  if (target_name == "edge") {
    node_1_ids <- sapply(
      dplyr::pull(data, param_names[1]),
      function(x) which(names(fit_edge$node_to_idx) == x)
    )
    node_2_ids <- sapply(
      dplyr::pull(data, param_names[2]),
      function(x) which(names(fit_edge$node_to_idx) == x)
    )
    dyad_ids <- get_dyad_ids(node_1_ids, node_2_ids, edgemodel$dyad_to_idx, edgemodel$directed)

    if (metric_name == "weight") {
      random_idxs <- sample(1:edgemodel$num_dyads, size=num_draws, replace=FALSE)
      posterior_samples <- extract_metric(edgemodel, paste0(target_name, "_", metric_name), num_draws)[random_idxs, dyad_ids]
    }
  } else if (target_name == "node") {
    node_ids <- sapply(
      dplyr::pull(data, param_names[1]),
      function(x) which(names(fit_edge$node_to_idx) == x)
    )
    posterior_samples <- extract_metric(edgemodel, paste0(target_name, "_", metric_name), num_draws)[, node_ids]
  } else {
    stop("Network feature not supported.")
  }

  # Make sure posterior samples are in matrix format (for num_draws=1).
  posterior_samples <- matrix(posterior_samples, nrow=num_draws)

  # Generate list of dataframes.
  data_list <- lapply(1:num_draws, function(i) {
    new_data <- data
    new_data[new_bison_term] <- posterior_samples[i, ]
    new_data
  })

  # Run brms imputation.
  brms::brm_multiple(brms_formula, data_list, backend="cmdstanr", ...)
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
