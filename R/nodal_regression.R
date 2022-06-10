require(bayesplot)
#' Fit a nodel regression model
#'
#' @param formula
#' @param edgemodel
#' @param df
#' @param mc_cores
#'
#' @return
#' @export
#'
#' @examples
nodal_regression <- function(formula, edgemodel, df, mc_cores=4, refresh=500, priors=NULL) {
  # If user-specified priors haven't been set, use the defaults
  if (is.null(priors)) {
    priors <- get_default_priors("nodal_regression")
  }

  design_matrices <- build_design_matrix_dyadic(formula, df)

  num_nodes <- edgemodel$num_nodes
  num_fixed <- ncol(design_matrices$X)
  num_random <- ncol(design_matrices$Z)
  num_random_groups <- length(unique(design_matrices$G))

  model_spec <- get_nodal_regression_spec(formula)
  metric_samples <- draw_node_metric_samples(edgemodel, model_spec$metric_name, 1000)
  metric_mu <- apply(metric_samples, 2, mean)
  metric_cov <- cov(metric_samples)

  model_data <- list(
    num_nodes=num_nodes,
    num_fixed=num_fixed,
    num_random=num_random,
    num_random_groups=num_random_groups,
    centrality_mu=metric_mu,
    centrality_cov=metric_cov,
    design_fixed=design_matrices$X,
    design_random=design_matrices$Z,
    random_group_index=design_matrices$G
  )

  # Set the priors in model data
  prior_parameters <- extract_prior_parameters(priors)
  model_data <- c(model_data, prior_parameters)

  model <- build_stan_model("nodal_regression")
  fit <- model$sample(data=model_data, chains=4, parallel_chains=mc_cores, refresh=refresh, step_size=0.1)
  chain <- fit$draws("beta_fixed", format="matrix")

  obj <- list()
  obj$formula <- formula
  obj$edgemodel <- edgemodel
  obj$fit <- fit
  obj$num_nodes <- num_nodes
  obj$model_data <- model_data
  obj$design_matrices <- design_matrices
  obj$chain <- chain
  obj$metric_samples <- metric_samples
  class(obj) <- "nodal_model"
  obj
}

#' Print information about a fitted nodal regression model
#'
#' @param obj
#'
#' @return
#' @export
#'
#' @examples
print.nodal_model <- function(obj) {
  coefficients <- t(apply(obj$chain, 2, function(x) quantile(x, probs=c(0.5, 0.05, 0.95))))
  rownames(coefficients) <- colnames(obj$design_matrices$X)
  coefficients <- round(coefficients, 3)
  cat(paste0(
    "=== Fitted nodal regression model ===\n",
    "Formula: ", format(obj$formula), "\n",
    "Number of nodes: ", obj$edgemodel$num_nodes, "\n",
    "=== Coefficient summary ==="
  ))
  print(coefficients)
}

#' Summary of a fitted nodal regression model
#'
#' @param obj
#'
#' @return
#' @export
#'
#' @examples
summary.nodal_model <- function(obj) {
  print(obj)
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
  metric_preds <- obj$fit$draws("centrality_pred", format="matrix")

  # Generate densities for edge sample and prediction
  sample_densities <- list()
  pred_densities <- list()
  for (i in 1:num_draws) {
    pred_densities[[i]] <- density(metric_samples[i, ])
    sample_densities[[i]] <- density(as.vector(metric_preds[i, ]))
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
    lines(pred_densities[[i]], col=rgb(0, 0, 1, 0.5))
  }
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


