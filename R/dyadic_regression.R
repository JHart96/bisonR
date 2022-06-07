require(bayesplot)

#' Fit dyadic regression model
#'
#' @param formula Formula of the model
#' @param edgemodel Fitted edge weight model
#' @param df Dataframe for regression
#' @param mc_cores Number of cores to use for the MCMC sampler
#' @param refresh Frequency of print-outs from MCMC sampler
#' @param mm TRUE/FALSE indicating whether to include multi-membership effects in the regression
#'
#' @return Fitted dyadic regression model object
#' @export
#'
#' @examples
dyadic_regression <- function(formula, edgemodel, df, mc_cores=4, refresh=500, mm=TRUE, priors=NULL) {
  # If user-specified priors haven't been set, use the defaults
  if (is.null(priors)) {
    priors <- get_default_priors("dyadic_regression")
  }

  design_matrices <- build_design_matrix(formula, df)

  num_nodes <- edgemodel$num_nodes
  N <- edgemodel$num_dyads
  K_fixed <- ncol(design_matrices$X)
  K_random <- ncol(design_matrices$Z)
  R <- length(unique(design_matrices$G))

  node_ids_1 <- edgemodel$node_to_idx[as.vector(df[all.vars(formula)[1]][, 1])]
  node_ids_2 <- edgemodel$node_to_idx[as.vector(df[all.vars(formula)[2]][, 1])]

  dyad_ids <- edgemodel$dyad_mapping[cbind(node_ids_1, node_ids_2)]
  edge_samples <- edgemodel$chain[, dyad_ids]
  edge_mu <- apply(edge_samples, 2, mean)
  edge_cov <- cov(edge_samples)

  model_data <- list(
    N=N,
    num_nodes=num_nodes,
    K_fixed=K_fixed,
    K_random=K_random,
    edge_mu=edge_mu,
    edge_cov=edge_cov,
    X=design_matrices$X,
    Z=design_matrices$Z,
    R=R,
    G=design_matrices$G,
    node_ids_1=node_ids_1,
    node_ids_2=node_ids_2,
    include_multimembership=as.integer(mm)
  )

  # Set the priors in model data
  prior_parameters <- extract_prior_parameters(priors)
  model_data <- c(model_data, prior_parameters)

  model <- build_stan_model("dyadic_regression")
  fit <- model$sample(data=model_data, chains=4, parallel_chains=mc_cores, refresh=refresh, step_size=0.1)
  chain <- fit$draws("beta_fixed", format="matrix")
  obj <- list()
  obj$formula <- formula
  obj$edgemodel <- edgemodel
  obj$fit <- fit
  obj$model_data <- model_data
  obj$design_matrices <- design_matrices
  obj$chain <- chain
  obj$edge_samples <- edge_samples
  obj$dyad_ids <- dyad_ids
  class(obj) <- "dyadic_model"
  obj
}

#' Print information about a fitted dyadic regression model
#'
#' @param obj
#'
#' @return
#' @export
#'
#' @examples
print.dyadic_model <- function(obj) {
  coefficients <- t(apply(obj$chain, 2, function(x) quantile(x, probs=c(0.5, 0.05, 0.95))))
  rownames(coefficients) <- colnames(obj$design_matrices$X)
  coefficients <- round(coefficients, 3)
  cat(paste0(
    "=== Fitted dyadic regression model ===\n",
    "Formula: ", obj$formula, "\n",
    "Number of dyads: ", obj$edgemodel$num_dyads, "\n",
    "=== Coefficient summary ==="
  ))
  print(coefficients)
}

#' Summary of a fitted dyadic regression model
#'
#' @param obj
#'
#' @return
#' @export
#'
#' @examples
summary.dyadic_model <- function(obj) {
  print(obj)
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
  edge_samples <- obj$edge_samples
  edge_preds <- obj$fit$draws("edge_pred", format="matrix")

  # Generate densities for edge sample and prediction
  sample_densities <- list()
  pred_densities <- list()
  for (i in 1:num_draws) {
    df_draw <- data.frame(y=as.vector(edge_samples[i, ]), dyad_ids=obj$dyad_ids)
    df_summed <- aggregate(y ~ as.factor(dyad_ids), df_draw, sum)
    pred_densities[[i]] <- density(df_summed$y)
    df_draw$y <- as.vector(edge_preds[i, ])
    df_summed <- aggregate(y ~ as.factor(dyad_ids), df_draw, sum)
    sample_densities[[i]] <- density(df_summed$y)
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
    lines(pred_densities[[i]], col=rgb(0, 0, 1, 0.5))
  }

}

build_design_matrix <- function(formula, data) {
  # Builds two design matrices for fixed and random effects.
  # Also builds a group vector for the random effects.
  # Don't forget intercept.

  X <- data.frame(empty_col = rep(0, nrow(data)))
  Z <- data.frame(empty_col = rep(0, nrow(data)))

  lpar <- list()
  hyper <- c() # Vector of hyperparameters.

  # Get model specification
  model_spec <- get_dyadic_regression_spec(formula)

  # If there is an intercept
  if (model_spec$intercept) {
    X[, "intercept"] <- 1
  }

  # print(model_terms)

  # Variable grouping for random effects
  G <- c()

  # Get additional fixed effects
  if (!is.null(model_spec$fixed)) {
    for (term_name in model_spec$fixed) {
      # If it's a factor, create a column for each level.
      if (is.factor(data[, term_name])) {
        var_group <- paste0("fixed_", term_name)
        term_levels <- levels(data[, term_name])
        for (term_level in term_levels) {
          new_term_name <-  paste0("fixed_", term_name, term_level)
          X[, new_term_name] <- 1 * (data[, term_name] == term_level)
        }
      } else {
        # Otherwise, create a single column:
        new_term_name <- paste0(c("fixed_", term_name), collapse="")
        X[, new_term_name] <- data[, term_name]
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
        Z[, new_term_name] <- 1 * (as.factor(data[, term_name]) == term_level)
        G[length(G) + 1] <- var_group
      }
    }
  }

  lpar$hyper <- hyper

  return(list(X=as.matrix(X[, -1]), Z=as.matrix(Z[, -1]), G=G))
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
