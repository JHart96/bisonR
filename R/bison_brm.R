
utils::globalVariables(c(".imp"))

#' Fit a brms model on BISoN network data
#'
#' @param formula A formula compatible with a brms model
#' @param edgemodel_list A bisonR edge model (or a list of them)
#' @param data_list A dataframe of regression variables compatible with a brms model (or a list of them)
#' @param num_draws Number of draws from the network posterior to use when fitting model
#' @param z_score Whether to standardise the network variable within each
#'   posterior draw. Equivalent to `scale = "z"`.
#' @param cores Number of computational cores to use.
#' @param chains Number of chains to run per model.
#' @param scale How to transform the network variable within each posterior
#'   draw: `"none"` (default), `"centre"` (subtract the mean) or `"z"` (subtract
#'   the mean and divide by the standard deviation). Transformations are applied
#'   across the rows used in the model, separately for each network.
#' @param subset An optional logical expression, evaluated in each data frame,
#'   selecting the rows to use, e.g. `subset = sex == "F"`. Network metrics are
#'   still calculated from the whole network; only the rows of the regression
#'   are selected. Not available for global metrics.
#' @param ... Additional arguments to be passed to brm(), such as family or priors
#'
#' @return A fitted brms model
#' @examples
#' sim_data <- simulate_bison_model("binary", aggregated = TRUE)
#' df <- sim_data$df_sim
#' fit_edge <- bison_model(
#'   (event | effort) ~ dyad(node_1_id, node_2_id),
#'   data = df,
#'   model_type = "binary_conjugate",
#'   priors = get_default_priors("binary_conjugate")
#' )
#' # Does age difference predict edge weight? The multi-membership term accounts
#' # for edges that share a node not being independent.
#' fit_dyadic <- bison_brm(
#'   bison(edge_weight(node_1_id, node_2_id)) ~ age_diff + (1 | mm(node_1_id, node_2_id)),
#'   fit_edge,
#'   df,
#'   num_draws = 5 # Use more (e.g. 100) in a real analysis
#' )
#' summary(fit_dyadic)
#' @export
bison_brm <- function(formula, edgemodel_list, data_list, num_draws=100, z_score=FALSE, cores=4, chains=2,
                      scale = c("none", "centre", "z"), subset = NULL, ...) {
  check_suggested("brms", "bison_brm")
  check_suggested("mice", "bison_brm")
  scale <- match.arg(scale)
  # Parse formula
  parsed_formula <- parse_bison_brms_formula(formula)
  brms_formula <- parsed_formula$brms_formula
  data_list <- subset_bison_data(data_list, substitute(subset), parent.frame(), parsed_formula$target_name)

  if (parsed_formula$target_name == "edge" && !grepl("mm\\(", brms_formula)) {
    warning(
      "Dyadic regressions should normally account for the non-independence of edges that share ",
      "a node, for example by adding a multi-membership term: + (1 | mm(",
      parsed_formula$param_names[1], ", ", parsed_formula$param_names[2], ")).",
      call. = FALSE
    )
  }

  # Generate mice object
  mice_obj <- bison_mice(
    edgemodel_list,
    data_list,
    parsed_formula$param_names,
    parsed_formula$target_name,
    parsed_formula$metric_name,
    num_draws,
    z_score,
    scale = scale
  )

  # Run brms imputation. brm_multiple() warns that E-BFMI can't be computed for
  # the combined fit; that warning isn't informative here, so it is muffled.
  withCallingHandlers(
    brms::brm_multiple(brms_formula, mice_obj, backend = "cmdstanr", cores = cores, chains = chains, ...),
    warning = function(w) {
      if (grepl("E-BFMI not computed", conditionMessage(w), fixed = TRUE)) invokeRestart("muffleWarning")
    }
  )
}

#' Get priors for BISoN brm model
#'
#' Wrapper for brms::get_prior() for use with bison_brm.
#'
#' @param formula A formula compatible with a brms model
#' @param edgemodel_list A bisonR edge model (or a list of them)
#' @param data_list A dataframe of regression variables compatible with a brms model (or a list of them)
#' @param scale,subset As for [bison_brm()]. The scale of the network variable
#'   affects the default priors, so use the same value as in [bison_brm()].
#' @param ... Additional parameters to be passed to brms::get_prior()
#'
#' @return The output from brms::get_prior()
#'
#' @examples
#' sim_data <- simulate_bison_model("binary", aggregated = TRUE)
#' df <- sim_data$df_sim
#' fit_edge <- bison_model(
#'   (event | effort) ~ dyad(node_1_id, node_2_id),
#'   data = df,
#'   model_type = "binary_conjugate",
#'   priors = get_default_priors("binary_conjugate")
#' )
#' bison_brm_get_prior(
#'   bison(edge_weight(node_1_id, node_2_id)) ~ age_diff + (1 | mm(node_1_id, node_2_id)),
#'   fit_edge,
#'   df
#' )
#' @export
bison_brm_get_prior <- function(formula, edgemodel_list, data_list, scale = c("none", "centre", "z"),
                                subset = NULL, ...) {
  scale <- match.arg(scale)
  check_suggested("brms", "bison_brm_get_prior")
  check_suggested("mice", "bison_brm_get_prior")
  # Parse formula
  parsed_formula <- parse_bison_brms_formula(formula)
  brms_formula <- parsed_formula$brms_formula
  data_list <- subset_bison_data(data_list, substitute(subset), parent.frame(), parsed_formula$target_name)

  # Generate mice object
  mice_obj <- bison_mice(
    edgemodel_list,
    data_list,
    parsed_formula$param_names,
    parsed_formula$target_name,
    parsed_formula$metric_name,
    num_draws=2,
    scale = scale
  )

  data <- mice::complete(mice_obj)

  brms::get_prior(brms_formula, data, ...)
}

parse_bison_brms_formula <- function(formula) {
  # Parse formula string
  formula_string <- deparse1(formula)[1]

  # Locate bison(...) term
  locs <- stringr::str_locate(formula_string, "bison\\([a-zA-Z0-9_.,\\[\\]\\(\\)\\s]*\\)")
  if (is.na(locs[1])) {
    stop("The formula must contain a bison(...) term, e.g. bison(node_strength(node_id)).", call. = FALSE)
  }
  formula_string_pre <- substr(formula_string, 1, locs[1] - 1)
  formula_string_post <- substr(formula_string, locs[2] + 1, nchar(formula_string))

  # Extract bison(...) term
  bison_term <- stringr::str_extract(formula_string, "bison\\([a-zA-Z0-9_.,\\[\\]\\(\\)\\s]*\\)")
  bison_term <- stringr::str_replace_all(bison_term, "\\s", "")

  # Extract function name and arguments
  bison_term_split <- stringr::str_split(bison_term, "bison\\(|\\(|\\)")[[1]]
  bison_term_split <- bison_term_split[nchar(bison_term_split) > 0]
  metric_name <- bison_term_split[1]
  params <- bison_term_split[2]

  # Extract parameter names
  param_names <- stringr::str_split(params, ",")[[1]]
  metric_parts <- stringr::str_split(metric_name, "_", n = 2)[[1]]
  if (length(metric_parts) != 2 || !metric_parts[1] %in% c("edge", "node", "global")) {
    stop(
      "Network metrics in bison() must start with edge_, node_ or global_, e.g. ",
      "bison(node_strength(node_id)), but got \"", metric_name, "\".",
      call. = FALSE
    )
  }
  target_name <- metric_parts[1]
  metric_name <- metric_parts[2]

  # Generate bison term for brms model.
  new_bison_term <- bison_term_name(target_name, metric_name)
  brms_formula <- paste0(formula_string_pre, new_bison_term, formula_string_post)

  return(list(
    param_names=param_names,
    target_name=target_name,
    metric_name=metric_name,
    new_bison_term=new_bison_term,
    brms_formula=brms_formula
  ))
}

# Name of the data column holding a network metric in bison_brm() models, e.g.
# "bison_node_strength" or "bison_node_degree_0_5" for node_degree[0.5].
bison_term_name <- function(target_name, metric_name) {
  name <- paste0("bison_", target_name, "_", metric_name)
  name <- gsub("[^A-Za-z0-9_]+", "_", name)
  sub("_+$", "", name)
}

# Select rows of the regression data with a logical expression evaluated in
# each data frame.
subset_bison_data <- function(data_list, subset_expr, env, target_name) {
  if (is.null(subset_expr)) {
    return(data_list)
  }
  if (target_name == "global") {
    stop("`subset` can't be used with global metrics; select the edge models and rows directly.", call. = FALSE)
  }
  select_rows <- function(df) {
    keep <- eval(subset_expr, df, env)
    if (!is.logical(keep) || length(keep) != nrow(df)) {
      stop("`subset` must be a logical expression with one value per row of the data.", call. = FALSE)
    }
    keep[is.na(keep)] <- FALSE
    if (!any(keep)) {
      stop("`subset` doesn't select any rows.", call. = FALSE)
    }
    df[keep, , drop = FALSE]
  }
  if (is.data.frame(data_list)) {
    return(select_rows(data_list))
  }
  lapply(data_list, select_rows)
}
