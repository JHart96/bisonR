data {
  int<lower=0> N; // Number of data points
  int<lower=0> K_fixed; // Number of fixed effect parameters
  int y[N]; // Outcome for each data point (presence/absence)
  int duration000[N]; // Duration of each observation
  matrix[N, K_fixed] X; // Design matrix for fixed effects
}

parameters {
  vector[K_fixed] beta_fixed; // Parameters for fixed effects.
}

model {
  vector[N] lprobs;
  lprobs = X * beta_fixed;

  // Main model
  target += binomial_lpmf(y | duration000, inv_logit(lprobs));

  // Priors
  target += normal_lpdf(beta_fixed | 0, 1);
}
