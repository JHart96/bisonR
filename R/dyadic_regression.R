#' Fit dyadic regression model
#'
#' @param formula Formula of the model
#' @param edgemodel Fitted edge weight model
#' @param df Dataframe for regression
#' @param mc_cores Number of cores to use for the MCMC sampler
#' @param refresh Frequency of print-outs from MCMC sampler
#' @param mm TRUE/FALSE indicating whether to include multi-membership effects in the regression
#' @param priors List of priors in the format supplied by `get_default_priors()`.
#'
#' @return An S3 dyadic model object containing chain samples and processed data.
#' @export
#'
#' @details
#' Fits a dyadic regression mixed model of the form where edge weight (with uncertainty) is either a response or predictor.
dyadic_regression <- function(formula, edgemodel, df, mc_cores=4, refresh=500, mm=TRUE, priors=NULL) {
  # If user-specified priors haven't been set, use the defaults
  if (is.null(priors)) {
    priors <- get_default_priors("dyadic_regression")
  }

  # design_matrices <- build_design_matrix(formula, df)
  model_info <- get_dyadic_regression_model_data(formula, edgemodel, df)
  model_data <- model_info$model_data
  model_data$include_multimembership = as.integer(mm)

  # Set the priors in model data
  prior_parameters <- extract_prior_parameters(priors)
  model_data <- c(model_data, prior_parameters)

  model <- build_stan_model("dyadic_regression")
  fit <- model$sample(data=model_data, chains=4, parallel_chains=mc_cores, refresh=refresh, step_size=0.1)
  chain <- fit$draws("beta_fixed", format="matrix")

  obj <- list(
    formula = formula,
    edgemodel = edgemodel,
    fit = fit,
    model_data = model_data,
    chain = chain,
    dyad_ids = model_info$dyad_ids
  )

  class(obj) <- "dyadic_model"
  obj
}

#' Print information about a dyadic regression model
#'
#' @param x An S3 dyadic regression model
#' @param digits Number of digits for rounding coefficients.
#' @param ... Additional arguments to be passed to summary calculations.
#'
#' @export
print.summary.dyadic_model <- function(x, digits=3, ...) {
  cat(x$description)
  coefficients <- round(x$coefficients, digits)
  print(coefficients)
}

#' Summary of a fitted dyadic regression model
#'
#' @param object An S3 dyadic regression model
#' @param ci Credible interval to use in summary, based on quantiles.
#' @param ... Additional arguments to be passed to summary calculations.
#'
#' @return Returns a summary object of the fitted dyadic regression model.
#'
#' @export
#'
summary.dyadic_model <- function(object, ci=0.90, ...) {
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
  class(summary_obj) <- "summary.dyadic_model"
  summary_obj
}

plot_trace.dyadic_model <- function (obj, par_ids=1:12, ...) {
  if (dim(obj$chain)[2] < 12) {
    par_ids <- 1:dim(obj$chain)[2]
  }
  bayesplot::mcmc_trace(obj$fit$draws("beta_fixed")[,,par_ids])
}

plot_predictions.dyadic_model <- function(obj, num_draws=20) {
  # Determine edge label
  if (obj$edgemodel$data_type == "binary") {
    xlab <- "Logit edge weight"
  } else if (obj$edgemodel$data_type == "binary") {
    xlab <- "Log edge weight"
  }

  # Extract edge samples and predictions
  edge_samples <- obj$edgemodel$edge_samples
  edge_preds <- obj$fit$draws("edge_pred", format="matrix")

  # Generate densities for edge sample and prediction
  sample_densities <- list()
  pred_densities <- list()
  for (i in 1:num_draws) {
    df_draw <- data.frame(y=as.vector(edge_samples[i, ]), dyad_id=obj$dyad_ids)
    df_summed <- aggregate(y ~ as.factor(dyad_id), df_draw, sum)
    sample_densities[[i]] <- density(df_summed$y)
    df_draw$y <- as.vector(edge_preds[i, ])
    df_summed <- aggregate(y ~ as.factor(dyad_id), df_draw, sum)
    pred_densities[[i]] <- density(df_summed$y)
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
      plot(sample_densities[[i]], main="Observed vs predicted edge weight", xlab=xlab, col=rgb(0, 0, 0, 0.5), xlim=xlim, ylim=ylim)
    } else {
      lines(sample_densities[[i]], col=rgb(0, 0, 0, 0.5))
    }
    lines(pred_densities[[i]], col=col2rgba(bison_colors[1], 0.5))
  }
  legend("topright", legend=c("observed", "predicted"), fill=c("black", bison_colors[1]))
}

get_dyadic_regression_model_data <- function(formula, edgemodel, data) {
  design_fixed <- data.frame(empty_col = rep(0, nrow(data)))
  design_random <- data.frame(empty_col = rep(0, nrow(data)))

  # Get model specification
  model_spec <- get_dyadic_regression_spec(formula)
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

  # Should come from model spec
  node_ids_1 <- edgemodel$node_to_idx[dplyr::pull(data, model_spec$node_1_name)]
  node_ids_2 <- edgemodel$node_to_idx[dplyr::pull(data, model_spec$node_2_name)]

  dyad_ids <- edgemodel$dyad_to_idx[cbind(node_ids_1, node_ids_2)]
  edge_samples <- edgemodel$chain[, dyad_ids]
  edge_mu <- apply(edge_samples, 2, mean)
  edge_cov <- cov(edge_samples)

  model_data <- list(
    num_rows = edgemodel$num_dyads,
    num_nodes = num_nodes,
    num_fixed = num_fixed,
    num_random = num_random,
    num_random_groups = num_random_groups,
    edge_mu = edge_mu,
    edge_cov = edge_cov,
    design_fixed = as.matrix(design_fixed[, -1]),
    design_random = as.matrix(design_random[, -1]),
    node_ids_1 = node_ids_1,
    node_ids_2 = node_ids_2,
    random_group_index = random_group_index
  )

  model_info <- list(
    model_data = model_data,
    dyad_ids = dyad_ids
  )

  return(model_info)
}

get_dyadic_regression_spec <- function(formula) {
  model_spec <- list()

  x <- str_split(deparse1(formula), "~")[[1]]
  lhs <- x[1]
  rhs <- x[2]

  # Remove whitespace
  lhs <- stringr::str_replace_all(lhs, " ", "")
  rhs <- stringr::str_replace_all(rhs, " ", "")

  # If lhs is a dyad term
  if (!is.na(str_match(lhs, "^dyad\\(.*,.*\\)$")[[1]])) {
    # dyad(,) term
    model_spec$dyad_side <- "lhs"
    node_names <- str_split(lhs, "\\(|\\)")[[1]][2]
    node_names <- str_replace_all(node_names, " ", "")
    node_names_split <- str_split(node_names, ",")[[1]]
    model_spec$node_1_name <- node_names_split[1]
    model_spec$node_2_name <- node_names_split[2]
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
    } else if (!is.na(str_match(term, "^dyad\\(.*,.*\\)$")[[1]])) {
      # dyad(,) term
      model_spec$dyad_side <- "rhs"
      node_names <- str_split(term, "\\(|\\)")[[1]][2]
      node_names <- str_replace_all(node_names, " ", "")
      node_names_split <- str_split(node_names, ",")[[1]]
      model_spec$node_1_name <- node_names_split[1]
      model_spec$node_2_name <- node_names_split[2]
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
