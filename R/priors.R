#' Get the default priors for an edge model
#'
#' Returns the default priors for a type of edge model, as a starting point for
#' choosing priors. A prior describes which values of a parameter are plausible
#' before seeing the data. The defaults may not suit every study, so check
#' them with [prior_check()] and [prior_predictive_check()] and change them
#' where needed before passing them to [bison_model()].
#'
#' @param model_type The type of model: `"binary"`, `"count"`, `"binary_conjugate"` or `"count_conjugate"`.
#'
#' @details Each prior is a string naming a distribution and its parameters,
#' such as `"normal(0, 2.5)"`. The entries are:
#' * `edge`: the edge weights. For `"binary"` and `"count"` models this is a
#'   normal distribution on the link scale (log-odds for binary models, log
#'   rates for count models). With `partial_pooling = TRUE`, the second value is
#'   instead the scale of a half-normal prior on the spread of the edge weights.
#'   For `"binary_conjugate"` it is a beta distribution on the probability, and
#'   for `"count_conjugate"` a gamma distribution (shape, rate) on the rate.
#' * `fixed`: fixed effects, as a normal distribution on the link scale.
#' * `random_std`: the standard deviation of random effects, as a half-normal
#'   distribution with the given scale.
#' * `zero_prob`: in zero-inflated models, the probability that a row of data
#'   is a structural zero (a zero regardless of the edge weight), as a beta
#'   distribution.
#'
#' Conjugate models only have an `edge` prior. Priors for terms that a model
#' doesn't include, such as `fixed` in a model without fixed effects, have no
#' effect.
#'
#' @return A named list of priors, in the format used by the `priors` argument
#'   of [bison_model()].
#'
#' @seealso [prior_check()], [prior_predictive_check()], [bison_model()]
#' @examples
#' priors <- get_default_priors("binary")
#' priors
#' # Priors are strings describing distributions, and can be changed:
#' priors$edge <- "normal(-1, 1.5)"
#' @export
get_default_priors <- function(model_type) {
  if (model_type == "binary") {
    priors <- list(
      edge="normal(0, 2.5)",
      fixed="normal(0, 2.5)",
      random_std="half-normal(1)",
      zero_prob="beta(1, 1)"
    )
    return(priors)
  }

  if (model_type == "count") {
    priors <- list(
      edge="normal(0, 1)",
      fixed="normal(0, 2.5)",
      random_std="half-normal(1)",
      zero_prob="beta(1, 1)"
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
  stop("Model type \"", model_type, "\" not supported. Use \"binary\", \"count\", \"binary_conjugate\" or \"count_conjugate\".", call. = FALSE)
}

#' Plot prior distributions
#'
#' Plots the density of each prior in a list, to check that the priors give
#' reasonable weight to plausible values. Edge weight priors are shown on the
#' natural scale (probabilities for binary models, rates for count models); other
#' priors are shown on the scale on which they are defined. To see the data
#' implied by the priors together, use [prior_predictive_check()].
#'
#' @param priors List of priors for a model, as returned by [get_default_priors()].
#' @param model_type Type of model the priors will be used for, as in
#'   [get_default_priors()].
#' @param type Type of plot. Only `"density"` is currently available.
#'
#' @details With `partial_pooling = TRUE` in [bison_model()], the second value of
#'   the edge prior is the scale of a half-normal prior on the spread of the edge
#'   weights, so the edge prior shown here is not the implied prior on each edge.
#'
#' @return `NULL`. The function is called for its plots.
#'
#' @seealso [get_default_priors()], [prior_predictive_check()]
#'
#' @examples
#' priors <- get_default_priors("binary")
#' prior_check(priors, "binary")
#' @export
prior_check <- function(priors, model_type, type = "density") {
  type <- match.arg(type, "density")
  if (type == "density") {
    if (length(priors) > 1) {
      num_cols = 2
    } else {
      num_cols = 1
    }
    num_rows = ceiling(length(priors)/num_cols)
    old_par <- par(mfrow = c(num_rows, num_cols), mar = c(4, 4, 1, 1))
    on.exit(par(old_par))
    for (parameter_name in names(priors)) {
      prior_distribution <- extract_distribution(priors[parameter_name])
      distribution_name <- prior_distribution$distribution_name
      parameter_values <- prior_distribution$parameter_values

      quantile_fn <- get_quantile_fn(parameter_name, model_type, distribution_name)
      density_fn <- get_density_fn(parameter_name, model_type, distribution_name)

      if (is.null(quantile_fn) || is.null(density_fn)) {
        stop("Distribution not available for this model type. Try using a different distribution in the priors.")
      }

      # Get limits on the X-axis
      x_lower <- do.call(quantile_fn, as.list(c(0.001, parameter_values)))
      x_upper <- do.call(quantile_fn, as.list(c(0.999, parameter_values)))

      x_values <- seq(x_lower, x_upper, (x_upper - x_lower)/100)
      y_values <- sapply(
        x_values,
        function(x) do.call(density_fn, as.list(c(x, parameter_values)))
      )
      plot(x_values, y_values, type="l", col=bison_colors[1], lwd=2, xlab=parameter_name, ylab="Prior probability")
    }
  }
}

get_density_fn <- function(parameter_name, model_type, distribution_name) {
  if (parameter_name == "edge") {
    if (model_type == "binary") {
      if (distribution_name == "normal") {
        # Density of a logit-normal distribution on the probability scale.
        return(function(x, ...) dnorm(qlogis(x), ...) / (x * (1 - x)))
      }
    }
    if (model_type == "count") {
      if (distribution_name == "normal") {
        return(function(x, ...) dlnorm(x, ...))
      }
    }
    if (model_type == "binary_conjugate") {
      if (distribution_name == "beta") {
        return(function(x, ...) dbeta(x, ...))
      }
    }
    if (model_type == "count_conjugate") {
      if (distribution_name == "gamma") {
        return(function(x, ...) dgamma(x, ...))
      }
    }
  } else {
    if (distribution_name == "normal") {
      return(function(x, ...) dnorm(x, ...))
    }
    if (distribution_name == "half-normal") {
      return(function(x, ...) extraDistr::dhnorm(x, ...))
    }
    if (distribution_name == "beta") {
      return(function(x, ...) dbeta(x, ...))
    }
  }
  return(NULL)
}

get_quantile_fn <- function(parameter_name, model_type, distribution_name) {
  if (parameter_name == "edge") {
    if (model_type == "binary") {
      if (distribution_name == "normal") {
        return(function(x, ...) plogis(qnorm(x, ...)))
      }

    }
    if (model_type == "count") {
      if (distribution_name == "normal") {
        return(function(x, ...) qlnorm(x, ...))
      }
    }
    if (model_type == "binary_conjugate") {
      if (distribution_name == "beta") {
        return(function(x, ...) qbeta(x, ...))
      }
    }
    if (model_type == "count_conjugate") {
      if (distribution_name == "gamma") {
        return(function(x, ...) qgamma(x, ...))
      }
    }
  } else {
    if (distribution_name == "normal") {
      return(function(x, ...) qnorm(x, ...))
    }
    if (distribution_name == "half-normal") {
      return(function(x, ...) extraDistr::qhnorm(x, ...))
    }
    if (distribution_name == "beta") {
      return(function(x, ...) qbeta(x, ...))
    }
  }

  return(NULL)
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

# Check that a list of priors is suitable for a model type, filling in any
# missing entries from the defaults.
validate_priors <- function(priors, model_type) {
  if (!is.list(priors) || is.null(names(priors)) || any(names(priors) == "")) {
    stop("`priors` must be a named list, like the output of get_default_priors().", call. = FALSE)
  }
  defaults <- get_default_priors(model_type)
  if ("random_mean" %in% names(priors)) {
    warning(
      "The `random_mean` prior is no longer used: random effects are centred on 0. ",
      "Add an intercept (`1`) to the formula if the model needs an overall level.",
      call. = FALSE
    )
    priors$random_mean <- NULL
  }
  unknown <- setdiff(names(priors), names(defaults))
  if (length(unknown) > 0) {
    stop(
      "Prior(s) ", paste0("`", unknown, "`", collapse = ", "), " aren't used by \"", model_type,
      "\" models. Expected: ", paste(names(defaults), collapse = ", "),
      ". Use get_default_priors(\"", model_type, "\") as a starting point.",
      call. = FALSE
    )
  }
  for (name in setdiff(names(defaults), names(priors))) {
    priors[[name]] <- defaults[[name]]
  }

  allowed <- list(
    edge = switch(model_type, binary_conjugate = "beta", count_conjugate = "gamma", "normal"),
    fixed = "normal",
    random_std = "half-normal",
    zero_prob = "beta"
  )
  for (name in names(priors)) {
    dist <- extract_distribution(priors[[name]])
    if (!identical(dist$distribution_name, allowed[[name]])) {
      stop(
        "The `", name, "` prior for \"", model_type, "\" models must be a ", allowed[[name]],
        " distribution, e.g. \"", defaults[[name]], "\", but got \"", priors[[name]], "\".",
        call. = FALSE
      )
    }
    expected_params <- if (allowed[[name]] == "half-normal") 1 else 2
    values <- dist$parameter_values
    if (length(values) != expected_params || anyNA(values)) {
      stop("Couldn't read the parameters of the `", name, "` prior \"", priors[[name]], "\".", call. = FALSE)
    }
    # Scale and shape parameters must be positive (the normal mean can be anything).
    positive <- if (allowed[[name]] == "normal") values[2] else values
    if (any(positive <= 0)) {
      stop(
        "The `", name, "` prior \"", priors[[name]], "\" has a parameter that must be positive ",
        "(a standard deviation, scale or shape).",
        call. = FALSE
      )
    }
  }
  priors[names(defaults)]
}
