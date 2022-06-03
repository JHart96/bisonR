#' Fit dyadic regression model
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
dyadic_regression <- function(formula, edgemodel, df, mc_cores=4, mm=TRUE) {
  design_matrices <- build_design_matrix(formula, df)

  num_nodes <- edgemodel$num_nodes
  N <- edgemodel$num_dyads
  K_fixed <- ncol(design_matrices$X)

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
    edge_mu=edge_mu,
    edge_cov=edge_cov,
    X=design_matrices$X,
    node_ids_1=node_ids_1,
    node_ids_2=node_ids_2,
    include_multimembership=as.integer(mm)
  )
  fit <- rstan::sampling(stanmodels$dyadic_regression, model_data, cores=mc_cores)
  chain <- rstan::extract(fit)
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
  coefficients <- t(apply(obj$chain[["beta_fixed"]], 2, function(x) quantile(x, probs=c(0.5, 0.05, 0.95))))
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

plot_trace.dyadic_model <- function (obj, ...) {
  rstan::traceplot(obj$fit, ...)
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
  edge_preds <- rstan::extract(obj$fit)$edge_pred

  # Generate densities for edge sample and prediction
  sample_densities <- list()
  pred_densities <- list()
  for (i in 1:num_draws) {
    df_draw <- data.frame(y=edge_samples[i, ], dyad_ids=obj$dyad_ids)
    df_summed <- aggregate(y ~ as.factor(dyad_ids), df_draw, sum)
    pred_densities[[i]] <- density(df_summed$y)
    df_draw$y <- edge_preds[i, ]
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

  # If there is an intercept
  if (attr(terms(formula), "intercept") == 1) {
    lpar$intercept = 0
    X[, "intercept"] <- 1
  }

  #

  model_terms <- labels(terms(formula))

  # print(model_terms)

  if (length(model_terms) > 0) {
    for (i in 1:length(model_terms)) {
      if (any(grep("dyad\\(.*,.*\\)", model_terms[i]))) {
        # Node multimembership effect
        var_group <- "node"
        lpar[[var_group]] <- c()

        # Get ID names in a different way. Can't rely on them being first in model terms.
        id_vars <- sapply(colnames(data), function(x) grepl(x, model_terms[i]))
        id_names <- names(id_vars[id_vars == TRUE])
        term_levels <- levels(as.factor(unique(c(data[, id_names[1]], data[, id_names[2]]))))

        for (term_level in term_levels) {
          new_term_name <-  paste0("node_", term_level)
          Z[, new_term_name] <- 1 * (as.factor(data[, id_names[1]]) == term_level)
          Z[, new_term_name] <- Z[, new_term_name] + 1 * (as.factor(data[, id_names[2]]) == term_level)

          lpar[[var_group]][[new_term_name]] <- 0
          hyper[paste0(var_group, "_mu")] <- 0
          hyper[paste0(var_group, "_sigma")] <- 0 # Will get exp.
        }

      } else if (any(grep("\\l", model_terms[i]))) {
        # Random effect

        # Extract term name
        term_vars <- sapply(colnames(data), function(x) grepl(x, model_terms[i]))
        term_var <- term_vars[term_vars == TRUE][1]
        term_name <- names(term_var)

        # Create variable grouping
        var_group <- paste0("random_", term_name)
        lpar[[var_group]] <- c()
        term_levels <- levels(as.factor(data[, term_name]))
        for (term_level in term_levels) {
          new_term_name <-  paste0("random_", term_name, term_level)
          Z[, new_term_name] <- 1 * (as.factor(data[, term_name]) == term_level)
          lpar[[var_group]][[new_term_name]] <- 0
          hyper[paste0(var_group, "_mu")] <- 0
          hyper[paste0(var_group, "_sigma")] <- 0 # Will get exp.
        }
      } else {
        # Fixed effect
        term_vars <- sapply(colnames(data), function(x) grepl(model_terms[i], x)) ## Change how this gets processed
        term_var <- term_vars[term_vars == TRUE][1]
        term_name <- names(term_var)
        # If it's a factor, create a column for each level.
        if (is.factor(data[, term_name])) {
          var_group <- paste0("fixed_", term_name)
          lpar[[var_group]] <- c()
          term_levels <- levels(data[, term_name])
          for (term_level in term_levels) {
            new_term_name <-  paste0("fixed_", term_name, term_level)
            X[, new_term_name] <- 1 * (data[, term_name] == term_level)
            lpar[[var_group]][[new_term_name]] <- 0
          }
        } else {
          # Otherwise, create a single column:
          new_term_name <- paste0(c("fixed_", term_name), collapse="")
          X[, new_term_name] <- data[, term_name]
          lpar[[new_term_name]] <- 0
        }
      }
    }
  }

  lpar$hyper <- hyper

  return(list(lpar=lpar, X=as.matrix(X[, -1]), Z=as.matrix(Z[, -1])))
}

# Processes formula and returns information about the model configuration
read_formula <- function(formula) {

}
