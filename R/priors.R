#' Return default priors for a given model
#'
#' @param model_type
#'
#' @return
#' @export
#'
#' @examples
get_default_priors <- function(model_type) {
  if (model_type == "binary") {
    priors <- list(
      edge="normal(0, 2.5)",
      fixed="normal(0, 2.5)",
      random_mean="normal(0, 1)",
      random_std="half-normal(1)"
    )
    return(priors)
  }
  if (model_type == "count") {
    priors <- list(
      edge="normal(0, 2.5)",
      fixed="normal(0, 2.5)",
      random_mean="normal(0, 1)",
      random_std="half-normal(1)"
    )
    return(priors)
  }
  if (model_type == "duration") {
    priors <- list(
      edge="normal(0, 2.5)",
      fixed="normal(0, 2.5)",
      random_mean="normal(0, 1)",
      random_std="half-normal(1)"
    )
    return(priors)
  }
  if (model_type == "dyadic_regression") {
    priors <- list(
      fixed="normal(0, 2.5)",
      random_mean="normal(0, 1)",
      random_std="half-normal(1)",
      multimembership="half-normal(1)",
      error="half-normal(1)"
    )
    return(priors)
  }
  if (model_type == "nodal_regression") {
    priors <- list(
      fixed="normal(0, 2.5)",
      random_mean="normal(0, 1)",
      random_std="half-normal(1)",
      error="half-normal(1)"
    )
    return(priors)
  }
  warning("Model type", model_type, "not supported.")
}

extract_prior_parameters <- function(priors) {
  prior_parameters <- list()
  for (parameter_name in names(priors)) {
    parameter_string <- stringr::str_replace_all(priors[parameter_name], " ", "")
    parameter_split <- stringr::str_split(parameter_string, "\\(|\\)|,")[[1]]
    distribution_name <- parameter_split[1]
    parameter_values <- parameter_split[2:(length(parameter_split) - 1)]
    parameter_values <- as.numeric(parameter_values)
    # prior_parameters[[parameter_name]] <- list(distribution=distribution_name, parameter_values=parameter_values)
    if (distribution_name == "normal") {
      prior_parameters[[paste("prior", parameter_name, "mu", sep="_")]] <- parameter_values[1]
      prior_parameters[[paste("prior", parameter_name, "sigma", sep="_")]] <- parameter_values[2]
    }
    if (distribution_name == "half-normal") {
      prior_parameters[[paste("prior", parameter_name, "sigma", sep="_")]] <- parameter_values[1]
    }
  }
  prior_parameters
}
