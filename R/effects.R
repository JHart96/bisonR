#' Get the estimates of other effects in an edge model
#'
#' Edge weight models can include effects other than the edge weights: fixed
#' effects, random effects, partial pooling and zero inflation (see
#' [bison_model()]). This function returns their posterior medians and credible
#' intervals in a data frame. A fixed effect is a single estimated effect of a
#' variable, such as age difference. A random effect estimates a separate value
#' for each level of a grouping variable, such as location, with the values
#' assumed to come from a shared normal distribution.
#'
#' @param obj A fitted `bison_model` object, fitted with Stan (conjugate models
#'   only estimate edge weights).
#' @param ci Width of the credible interval, e.g. `0.9` for a 90% interval.
#'
#' @details The `type` column says what each row is:
#' * `"fixed"`: a fixed effect, such as `age` for a numeric column, or
#'   `sex = F` for one level of a factor. These are on the same scale as the
#'   edge weights before they are transformed: log-odds for binary models and
#'   log rates for count models.
#' * `"random"`: the random effect for one level of a group, such as
#'   `location = 2`, on the same scale as fixed effects.
#' * `"random sd"`: the standard deviation of the random effects in a group.
#' * `"pooling sd"`: with `partial_pooling = TRUE`, the standard deviation of
#'   the edge weights.
#' * `"zero inflation"`: with `zero_inflated = TRUE`, the probability that an
#'   observation is a structural zero.
#'
#' @return A data frame with columns `effect`, `type`, `median`, and the lower
#'   and upper bounds of the credible interval (named after their quantiles,
#'   e.g. `5%` and `95%`). It has no rows if the model has no effects other than
#'   edge weights.
#'
#' @seealso [get_edgelist()] for the edge weights, and [summary.bison_model()],
#'   which shows both.
#'
#' @examplesIf !is.null(tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL))
#' sim_data <- simulate_bison_model("binary", aggregated = FALSE)
#' df <- sim_data$df_sim
#' df$location <- factor(df$location)
#' fit_edge <- bison_model(
#'   (event | effort) ~ dyad(node_1_id, node_2_id) + age_diff + (1 | location),
#'   data = df,
#'   model_type = "binary",
#'   priors = get_default_priors("binary")
#' )
#' get_effects(fit_edge)
#' @export
get_effects <- function(obj, ci = 0.9) {
  check_ci(ci)
  if (isTRUE(obj$conjugate) || is.null(obj$fit)) {
    stop("Conjugate models only estimate edge weights, so there are no other effects to report.", call. = FALSE)
  }
  md <- obj$model_data
  info <- obj$model_info
  # Fall back to the design matrix names for models fitted with older versions.
  fixed_labels <- if (length(info$fixed_labels) == md$num_fixed) info$fixed_labels else colnames(md$design_fixed)
  random_labels <- if (length(info$random_labels) == md$num_random) info$random_labels else colnames(md$design_random)
  group_labels <- if (length(info$random_group_labels) == md$num_random_groups) {
    info$random_group_labels
  } else {
    paste("group", seq_len(md$num_random_groups))
  }

  pieces <- list()
  add <- function(variable, type, labels) {
    pieces[[length(pieces) + 1]] <<- list(variable = variable, type = type, labels = labels)
  }
  if (md$num_fixed > 0) add("beta_fixed", "fixed", fixed_labels)
  if (md$num_random > 0) {
    add("beta_random", "random", random_labels)
    add("random_group_sigma", "random sd", paste0("sd(", group_labels, ")"))
  }
  if (isTRUE(md$partial_pooling == 1) && md$num_edges > 0) add("edge_sigma", "pooling sd", "sd(edge weights)")
  if (isTRUE(md$zero_inflated == 1)) add("zero_prob", "zero inflation", "zero_prob")

  lb <- 0.5 * (1 - ci)
  bound_names <- paste0(format_percent(c(lb, 1 - lb)), "%")
  if (length(pieces) == 0) {
    out <- data.frame(effect = character(0), type = character(0), median = numeric(0),
                      lower = numeric(0), upper = numeric(0))
    colnames(out)[4:5] <- bound_names
    return(out)
  }

  variables <- vapply(pieces, function(p) p$variable, character(1))
  draws <- as_plain_matrix(obj$fit$draws(variables, format = "matrix"))
  rows <- lapply(pieces, function(p) {
    columns <- grep(paste0("^", p$variable, "\\["), colnames(draws))
    quantiles <- edge_quantiles(draws[, columns, drop = FALSE], ci)
    data.frame(
      effect = p$labels,
      type = p$type,
      median = round(quantiles[1, ], 3),
      lower = round(quantiles[2, ], 3),
      upper = round(quantiles[3, ], 3)
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  colnames(out)[4:5] <- bound_names
  out
}
