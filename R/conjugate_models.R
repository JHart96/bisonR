# Conjugate edge models. These give exact posteriors for edge weights without
# MCMC, but only support dyad() terms.
#
# The returned chain is on the link scale (logit for binary, log for count) so
# it matches the Stan models. Predictions and log-likelihoods are per data row,
# so plot_predictions() and model_comparison() work as for Stan models.
fit_conjugate_model <- function(model_type, model_data, num_samples, priors_only = FALSE) {
  if (model_type == "binary") {
    return(fit_conjugate_model_binary(model_data, num_samples, priors_only))
  }
  if (model_type == "count") {
    return(fit_conjugate_model_count(model_data, num_samples, priors_only))
  }
  stop("Conjugate models are only available for binary and count data.", call. = FALSE)
}

# Total events and effort for each dyad (rows are summed over observations).
conjugate_totals <- function(model_data) {
  dyads <- factor(model_data$dyad_ids, levels = seq_len(model_data$num_edges))
  list(
    events = as.numeric(tapply(model_data$event, dyads, sum, default = 0)),
    divisor = as.numeric(tapply(model_data$divisor, dyads, sum, default = 0))
  )
}

fit_conjugate_model_binary <- function(model_data, num_samples, priors_only) {
  totals <- conjugate_totals(model_data)
  if (priors_only) {
    totals$events[] <- 0
    totals$divisor[] <- 0
  }
  alpha <- model_data$prior_edge_alpha + totals$events
  beta <- model_data$prior_edge_beta + totals$divisor - totals$events
  num_edges <- model_data$num_edges

  # num_samples x num_edges matrix of edge probabilities
  p <- matrix(
    stats::rbeta(num_samples * num_edges, rep(alpha, each = num_samples), rep(beta, each = num_samples)),
    num_samples, num_edges
  )

  # Row-level predictions and log-likelihoods
  p_rows <- p[, model_data$dyad_ids, drop = FALSE]
  divisor_rows <- matrix(model_data$divisor, num_samples, model_data$num_rows, byrow = TRUE)
  event_rows <- matrix(model_data$event, num_samples, model_data$num_rows, byrow = TRUE)
  event_preds <- matrix(stats::rbinom(length(p_rows), divisor_rows, p_rows), num_samples, model_data$num_rows)
  log_lik <- matrix(stats::dbinom(event_rows, divisor_rows, p_rows, log = TRUE), num_samples, model_data$num_rows)

  # Keep draws away from exactly 0 or 1, which can occur with very small prior
  # shapes, so the link-scale values are finite.
  p <- pmin(pmax(p, .Machine$double.eps), 1 - .Machine$double.eps)
  list(chain = stats::qlogis(p), event_preds = event_preds, log_lik = log_lik)
}

fit_conjugate_model_count <- function(model_data, num_samples, priors_only) {
  totals <- conjugate_totals(model_data)
  if (priors_only) {
    totals$events[] <- 0
    totals$divisor[] <- 0
  }
  shape <- model_data$prior_edge_alpha + totals$events
  rate <- model_data$prior_edge_beta + totals$divisor
  num_edges <- model_data$num_edges

  lambda <- matrix(
    stats::rgamma(num_samples * num_edges, shape = rep(shape, each = num_samples), rate = rep(rate, each = num_samples)),
    num_samples, num_edges
  )

  lambda_rows <- lambda[, model_data$dyad_ids, drop = FALSE]
  divisor_rows <- matrix(model_data$divisor, num_samples, model_data$num_rows, byrow = TRUE)
  event_rows <- matrix(model_data$event, num_samples, model_data$num_rows, byrow = TRUE)
  mu <- lambda_rows * divisor_rows
  event_preds <- matrix(stats::rpois(length(mu), mu), num_samples, model_data$num_rows)
  log_lik <- matrix(stats::dpois(event_rows, mu, log = TRUE), num_samples, model_data$num_rows)

  list(chain = log(pmax(lambda, .Machine$double.xmin)), event_preds = event_preds, log_lik = log_lik)
}
