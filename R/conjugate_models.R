fit_conjugate_model <- function(data_type, model_data, num_samples=10000) {
  if (data_type == "binary") {
    model <- fit_conjugate_model_binary(model_data, num_samples)
  }
  if (data_type == "count") {
    model <- fit_conjugate_model_count(model_data, num_samples)
  }
  return(model)
}

fit_conjugate_model_binary <- function(model_data, num_samples) {
  # Draw samples from beta posterior
  chain <- matrix(0, num_samples, model_data$num_rows)
  for (i in 1:model_data$num_rows) {
    chain[, model_data$dyad_ids[i]] <- rbeta(
      num_samples,
      model_data$prior_edge_alpha + model_data$event[i],
      model_data$prior_edge_beta + model_data$divisor[i] - model_data$event[i]
    )
  }

  # Draw predictions from binomial distribution
  event_preds <- t(sapply(1:num_samples, function(i) {
    rbinom(model_data$num_edges, model_data$divisor, chain[i, ])
  }))

  chain <- qlogis(chain)
  return(list(chain=chain, event_preds=event_preds))
}

fit_conjugate_model_count <- function(model_data, num_samples) {
  # Draw samples from gamma posterior
  chain <- matrix(0, num_samples, model_data$num_rows)
  for (i in 1:model_data$num_rows) {
    if (model_data$divisor[i] > 0) {
      chain[, model_data$dyad_ids[i]] <- rgamma(
        num_samples,
        model_data$prior_edge_alpha + model_data$event[i],
        model_data$prior_edge_beta/model_data$divisor[i] + 1
      )/model_data$divisor[i]
    } else {
      chain[, model_data$dyad_ids[i]] <- rgamma(
        num_samples,
        model_data$prior_edge_alpha,
        model_data$prior_edge_beta
      )
    }
  }

  # Draw predictions from Poisson distribution
  event_preds <- t(sapply(1:num_samples, function(i) {
    rpois(model_data$num_edges, chain[i, ] * model_data$divisor)
  }))

  chain <- log(chain)
  return(list(chain=chain, event_preds=event_preds))
}
