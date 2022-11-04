require(brms)

utils::globalVariables(c(".imp"))

#' Fit a brms model on BISoN network data
#'
#' @param formula A formula compatible with a brms model
#' @param edgemodel_list A bisonR edge model (or a list of them)
#' @param data_list A dataframe of regression variables compatible with a brms model (or a list of them)
#' @param num_draws Number of draws from the network posterior to use when fitting model
#' @param z_score Whether to Z-score bison variable or not.
#' @param cores Number of computational cores to use.
#' @param chains Number of chains to run per model.
#' @param ... Additional arguments to be passed to brm(), such as family or priors
#'
#' @return A fitted brms model
#' @export
bison_brm <- function(formula, edgemodel_list, data_list, num_draws=100, z_score=FALSE, cores=4, chains=2, ...) {
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
    num_draws,
    z_score
  )

  # Run brms imputation.
  brms::brm_multiple(brms_formula, mice_obj, backend="cmdstanr", cores=cores, chains=chains, ...)
}

#' Get priors for BISoN brm model
#'
#' Wrapper for brms::get_prior() for use with bison_brm.
#'
#' @param formula A formula compatible with a brms model
#' @param edgemodel_list A bisonR edge model (or a list of them)
#' @param data_list A dataframe of regression variables compatible with a brms model (or a list of them)
#' @param ... Additional parameters to be passed to brms::get_prior()
#'
#' @return The output from brms::get_prior()
#'
#' @export
bison_brm_get_prior <- function(formula, edgemodel_list, data_list, ...) {
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
    1
  )

  data <- mice::complete(mice_obj)

  brms::get_prior(brms_formula, data, ...)
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
