#' Return default priors for a given model
#'
#' @param model_type A character specifying the type of model to retrieve priors for ("binary", "count", "duration", etc)
#'
#' @return A list of default priors that can be used in the `edge_model()` function.
#' @export
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
  if (model_type == "binary_conjugate") {
    priors <- list(
      edge="beta(1, 1)"
    )
    return(priors)
  }
  if (model_type == "count_conjugate") {
    priors <- list(
      edge="gamma(1, 1)"
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

#' Prior checks
#'
#' @param priors List of priors for a model, can be retrieved using `get_default_prior()`.
#' @param model_type Type of model the priors will be used for (same as the argument for `get_default_prior()`).
#' @param type Type of prior check to run, `"value"` or `"prediction"`. Details below.
#'
#' @details The parameter `type` determines what type of prior check to run. `type="value"` will
#' plot the prior probability over the parameter value on the original scale. `type="predictive"`
#' will run a prior predictive plot, where predictions from the model are generated using only
#' prior probabilities (the model when not updated from the data).
#'
#' @export
prior_check <- function(priors, model_type, type="density") {
  if (type == "density") {
    num_cols = 2
    num_rows = ceiling(length(priors)/num_cols)
    par(mfrow=c(num_rows, num_cols))
    for (parameter_name in names(priors)) {
      prior_distribution <- extract_distribution(priors[parameter_name])
      distribution_name <- prior_distribution$distribution_name
      parameter_values <- prior_distribution$parameter_values

      density_fn <- get_density_fn(parameter_name, model_type, distribution_name)
      xlim <- get_xlim(function(x) do.call(density_fn, as.list(c(x, parameter_values)))) # This needs improving
      x_values <- seq(xlim[1], xlim[2], (xlim[2] - xlim[1])/100)
      y_values <- sapply(
        x_values,
        function(x) do.call(density_fn, as.list(c(x, parameter_values)))
      )
      plot(x_values, y_values, type="l", col=bison_colors[1], lwd=2, xlab=parameter_name, ylab="Prior probability")
    }
    par(mfrow=c(1, 1))
  }
}

# Get reasonable limits of density function. Could be improved a lot.
get_xlim <- function(density_fn) {
  start <- 0.5
  while(is.nan(density_fn(start)) || density_fn(start) <= 10^-3) {
    if (density_fn(start + 0.01) >= density_fn(start)) {
      start <- start + 0.01
    } else {
      start <- start - 0.01
    }
  }
  xmin <- start
  while(!is.nan(suppressWarnings(density_fn(xmin))) && density_fn(xmin) > 10^-3) {
    xmin <- xmin - 0.01
  }
  xmax <- start
  while(!is.nan(suppressWarnings(density_fn(xmax))) && density_fn(xmax) > 10^-3) {
    xmax <- xmax + 0.01
  }
  return(c(xmin + 0.01, xmax - 0.01))
}

get_density_fn <- function(parameter_name, model_type, distribution_name) {
  if (parameter_name == "edge") {
    if (model_type == "binary") {
      if (distribution_name == "normal") {
        density_fn <- function(x, ...) dnorm(qlogis(x), ...)
      }
    }
    if (model_type == "count") {
      if (distribution_name == "normal") {
        density_fn <- function(x, ...) dlnorm(x, ...)
      }
    }
  } else {
    if (distribution_name == "normal") {
      density_fn <- function(x, ...) dnorm(x, ...)
    }
    if (distribution_name == "half-normal") {
      density_fn <- function(x, ...) extraDistr::dhnorm(x, ...)
    }
  }

  return(density_fn)
}

extract_prior_parameters <- function(priors) {
  prior_parameters <- list()
  for (parameter_name in names(priors)) {
    prior_distribution <- extract_distribution(priors[parameter_name])
    distribution_name <- prior_distribution$distribution_name
    parameter_values <- prior_distribution$parameter_values
    if (distribution_name == "normal") {
      prior_parameters[[paste("prior", parameter_name, "mu", sep="_")]] <- parameter_values[1]
      prior_parameters[[paste("prior", parameter_name, "sigma", sep="_")]] <- parameter_values[2]
    }
    if (distribution_name == "half-normal") {
      prior_parameters[[paste("prior", parameter_name, "sigma", sep="_")]] <- parameter_values[1]
    }
    if (distribution_name == "beta") {
      prior_parameters[[paste("prior", parameter_name, "alpha", sep="_")]] <- parameter_values[1]
      prior_parameters[[paste("prior", parameter_name, "beta", sep="_")]] <- parameter_values[2]
    }
    if (distribution_name == "gamma") {
      prior_parameters[[paste("prior", parameter_name, "alpha", sep="_")]] <- parameter_values[1]
      prior_parameters[[paste("prior", parameter_name, "beta", sep="_")]] <- parameter_values[2]
    }
  }
  # print(prior_parameters)
  prior_parameters
}

# Extracts distribution name and parameter values for a prior string (e.g. from "normal(0, 1)")
extract_distribution <- function(prior_string) {
  parameter_string <- stringr::str_replace_all(prior_string, " ", "")
  parameter_split <- stringr::str_split(parameter_string, "\\(|\\)|,")[[1]]
  distribution_name <- parameter_split[1]
  parameter_values <- parameter_split[2:(length(parameter_split) - 1)]
  parameter_values <- as.numeric(parameter_values)
  return(list(distribution_name=distribution_name, parameter_values=parameter_values))
}
