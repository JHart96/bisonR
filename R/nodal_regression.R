#' Fit a nodel regression model
#'
#' @param formula Formula of the model
#' @param edgemodel Fitted edge weight model
#' @param df Dataframe for regression
#' @param mc_cores Number of cores to use for the MCMC sampler
#' @param refresh Frequency of print-outs from MCMC sampler
#' @param priors List of priors in the format supplied by `get_default_priors()`.
#'
#' @return An S3 nodal model object containing chain samples and processed data.
#' @export
#'
#' @examples
nodal_regression <- function(formula, edgemodel, df, mc_cores=4, refresh=500, priors=NULL) {
  # If user-specified priors haven't been set, use the defaults
  if (is.null(priors)) {
    priors <- get_default_priors("nodal_regression")
  }

  model_info <- get_nodal_regression_model_data(formula, edgemodel, df)
  model_data <- model_info$model_data

  # Set the priors in model data
  prior_parameters <- extract_prior_parameters(priors)
  model_data <- c(model_data, prior_parameters)

  model <- build_stan_model("nodal_regression")
  fit <- model$sample(data=model_data, chains=4, parallel_chains=mc_cores, refresh=refresh, step_size=0.1)
  chain <- fit$draws("beta_fixed", format="matrix")

  obj <- list(
    formula = formula,
    edgemodel = edgemodel,
    fit = fit,
    model_data = model_data,
    chain = chain,
    metric_samples = model_info$metric_samples
  )

  class(obj) <- "nodal_model"
  obj
}

get_nodal_regression_model_data <- function(formula, edgemodel, data) {
  design_fixed <- data.frame(empty_col = rep(0, nrow(data)))
  design_random <- data.frame(empty_col = rep(0, nrow(data)))

  # Get model specification
  model_spec <- get_nodal_regression_spec(formula)
  # If there is an intercept
  if (model_spec$intercept) {
    design_fixed[, "intercept"] <- 1
  }

  # Variable grouping for random effects
  random_group_index <- c()

  # Get additional fixed effects
  if (!is.null(model_spec$fixed)) {
    for (term_name in model_spec$fixed) {
      # If it's a factor, create a column for each level.
      if (is.factor(data[, term_name])) {
        var_group <- paste0("fixed_", term_name)
        term_levels <- levels(data[, term_name])
        for (term_level in term_levels) {
          new_term_name <-  paste0("fixed_", term_name, term_level)
          design_fixed[, new_term_name] <- 1 * (data[, term_name] == term_level)
        }
      } else {
        # Otherwise, create a single column:
        new_term_name <- paste0(c("fixed_", term_name), collapse="")
        design_fixed[, new_term_name] <- data[, term_name]
      }
    }
  }

  # Get additional random effects
  if (!is.null(model_spec$random)) {
    for (term_name in model_spec$random) {
      var_group <- paste0("random_", term_name)
      term_levels <- levels(as.factor(data[, term_name]))
      for (term_level in term_levels) {
        new_term_name <-  paste0("random_", term_name, term_level)
        design_random[, new_term_name] <- 1 * (as.factor(data[, term_name]) == term_level)
        random_group_index[length(random_group_index) + 1] <- var_group
      }
    }
  }

  num_nodes <- edgemodel$num_nodes
  num_dyads <- edgemodel$num_dyads
  num_fixed <- ncol(design_fixed) - 1
  num_random <- ncol(design_random) - 1
  num_random_groups <- length(unique(random_group_index))

  metric_samples <- draw_node_metric_samples(edgemodel, model_spec$metric_name, 1000, standardise=TRUE)
  metric_mu <- apply(metric_samples, 2, mean)
  metric_cov <- cov(metric_samples)

  model_data <- list(
    num_rows = edgemodel$num_dyads,
    num_nodes = num_nodes,
    num_fixed = num_fixed,
    num_random = num_random,
    num_random_groups = num_random_groups,
    metric_mu = metric_mu,
    metric_cov = metric_cov,
    design_fixed = as.matrix(design_fixed[, -1]),
    design_random = as.matrix(design_random[, -1]),
    random_group_index = random_group_index
  )

  model_info <- list(
    model_data = model_data,
    metric_samples = metric_samples
  )

  return(model_info)
}

#' Print information about a fitted nodal regression model
#'
#' @param x An S3 nodal regression model
#' @param digits Number of digits for rounding coefficients.
#' @param ... Additional arguments to be passed to summary calculations.
#'
#' @export
#'
print.summary.nodal_model <- function(x, digits=3, ...) {
  cat(x$description)
  coefficients <- round(x$coefficients, digits)
  print(coefficients)
}

#' Summary of a fitted nodal regression model
#'
#' @param object An S3 dyadic regression model
#' @param ci Credible interval to use in summary, based on quantiles.
#' @param ... Additional arguments to be passed to summary calculations.
#'
#' @export
#'
summary.nodal_model <- function(object, ci=0.90, ...) {
  summary_obj <- list()
  coefficients <- t(apply(object$chain, 2, function(object) quantile(object, probs=c(0.5, 0.5 * (1 - ci), ci + 0.5 * (1 - ci)))))
  rownames(coefficients) <- colnames(object$model_data$design_fixed)

  summary_obj$coefficients <- coefficients
  summary_obj$description <- paste0(
    "=== Fitted dyadic regression model ===\n",
    "Formula: ", format(object$formula), "\n",
    "Number of dyads: ", object$edgemodel$num_dyads, "\n",
    "=== Coefficient summary ==="
  )

  class(summary_obj) <- "summary.nodal_model"

  summary_obj
}

plot_trace.nodal_model <- function (obj, par_ids=1:12, ...) {
  if (dim(obj$chain)[2] < 12) {
    par_ids <- 1:dim(obj$chain)[2]
  }
  bayesplot::mcmc_trace(obj$fit$draws("beta_fixed")[,,par_ids])
}

plot_predictions.nodal_model <- function(obj, num_draws=20) {
  # Determine edge label
  xlab <- "Logit centrality"

  # Extract edge samples and predictions
  metric_samples <- obj$metric_samples
  metric_preds <- obj$fit$draws("metric_pred", format="matrix")

  # Generate densities for edge sample and prediction
  sample_densities <- list()
  pred_densities <- list()
  for (i in 1:num_draws) {
    sample_densities[[i]] <- density(metric_samples[i, ])
    pred_densities[[i]] <- density(as.vector(metric_preds[i, ]))
  }

  # Set plot limits according to maximum density of samples
  xlim = c(
    min(sapply(sample_densities, function(x) min(x$x))),
    max(sapply(sample_densities, function(x) max(x$x))) * 1.1
  )
  ylim = c(
    min(sapply(sample_densities, function(x) min(x$y))),
    max(sapply(sample_densities, function(x) max(x$y))) * 1.1
  )

  # Plot densities for subsequent draws
  for (i in 1:num_draws) {
    if (i == 1) {
      plot(sample_densities[[i]], main="Observed vs predicted nodal metrics", xlab=xlab, col=rgb(0, 0, 0, 0.5), xlim=xlim, ylim=ylim)
    } else {
      lines(sample_densities[[i]], col=rgb(0, 0, 0, 0.5))
    }
    lines(pred_densities[[i]], col=col2rgba(bison_colors[1], 0.5))
  }
  legend("topright", legend=c("observed", "predicted"), fill=c("black", bison_colors[1]))
}

extract_metric_name <- function(term) {
  metric_name <- stringr::str_split(term, "\\(")[[1]][1]
  metric_name <- stringr::str_replace_all(metric_name, " ", "")
  if(!is.null(get_metric_fn(metric_name))) {
    return(metric_name)
  }
  return(NULL)
}

get_nodal_regression_spec <- function(formula) {
  model_spec <- list()

  x <- str_split(deparse1(formula), "~")[[1]]
  lhs <- x[1]
  rhs <- x[2]

  # Remove whitespace
  lhs <- stringr::str_replace_all(lhs, " ", "")
  rhs <- stringr::str_replace_all(rhs, " ", "")

  # If lhs is a metric term
  if (!is.na(str_match(lhs, "[a-z]\\(.*\\)")[[1]])) {
    metric_name <- extract_metric_name(lhs)
    if (is.null(metric_name)) {
      stop(paste0("Invalid node metric name. Note that the use of functions within a formula is not supported, make sure any transformations are applied to the dataframe *before* fitting the model."))
    }
    model_spec$dyad_side <- "lhs"
    node_name <- stringr::str_split(lhs, "\\(|\\)")[[1]][2]
    model_spec$node_name <- node_name
    model_spec$metric_name <- metric_name
  }

  # Set intercept to true by default
  model_spec$intercept <- TRUE

  model_spec$fixed <- c()
  model_spec$random <- c()

  rhs_split <- str_split(rhs, "\\+")[[1]]
  for (term in rhs_split) {
    term <- str_replace_all(term, " ", "")
    # Is it an intercept, a dyad, a fixed effect, or a random effect?
    if (!is.na(str_match(term, "^0|1$"))) {
      # Intercept (or lack thereof)
      if (term == "0") {
        model_spec$intercept <- FALSE
      } else {
        model_spec$intercept <- TRUE
      }
    } else if (!is.na(str_match(term, "[a-z]\\(.*\\)")[[1]])) {
      metric_name <- extract_metric_name(term)
      if (is.null(metric_name)) {
        stop(paste0("Invalid node metric name. Note that the use of functions within a formula is not supported, make sure any transformations are applied to the dataframe *before* fitting the model."))
      }
      model_spec$dyad_side <- "rhs"
      node_name <- stringr::str_split(term, "\\(|\\)")[[1]][2]
      model_spec$metric_name <- metric_name
      model_spec$node_name <- node_name
    } else if (is.na(str_match(term, "[^a-zA-Z0-9_]"))) {
      # No non-alphanumeric characters, and it can't be an intercept, so it's a fixed effect
      model_spec$fixed[length(model_spec$fixed) + 1] <- term
    } else if (!is.na(str_match(term, "^\\(1\\|.*\\)$"))) {
      # Contains a (1 | *) structure, so it's a basic random effect
      term_name <- str_split(term, "\\(|\\||\\)")[[1]][3]
      model_spec$random[length(model_spec$random) + 1] <- term_name
    } else {
      warning(paste0("Formula term \"", term, "\" not supported by bisonR. Check the formula is correctly specified."))
    }
  }

  model_spec
}


