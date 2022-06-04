stan_model_binary_code <- "
data {
  int<lower=0> N; // Number of data points
  int<lower=0> K_fixed; // Number of fixed effect parameters
  int y[N]; // Outcome for each data point (presence/absence)
  int divisor[N]; // Duration of each observation
  matrix[N, K_fixed] X; // Design matrix for fixed effects
}

parameters {
  vector[K_fixed] beta_fixed; // Parameters for fixed effects.
}

transformed parameters {
  vector[N] lprobs;
  lprobs = X * beta_fixed;
}

model {
  // Main model
  target += binomial_lpmf(y | divisor, inv_logit(lprobs));

  // Priors
  target += normal_lpdf(beta_fixed | 0, 1);
}

generated quantities {
  int y_pred[N];
  y_pred = binomial_rng(divisor, inv_logit(lprobs));
}
"

model_filepath <- cmdstanr::write_stan_file(stan_model_binary_code)
model <- cmdstanr::cmdstan_model(model_filepath)
model_data <- fit_edge$model_data
model_data$node_names <- NULL
fit <- model$sample (
  data = model_data,
  chains = 4,
  parallel_chains = 4
)
