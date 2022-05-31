data {
  int<lower=0> N; // Number of data points
  int<lower=0> K_fixed; // Number of fixed effect parameters
  int y[N]; // Outcome for each data point (presence/absence)
  vector[N] divisor; // Duration of each observation
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
  target += poisson_lpmf(y | exp(lprobs) .* divisor);

  // Priors
  target += normal_lpdf(beta_fixed | 0, 2.5);
}

generated quantities {
  int y_pred[N];
  y_pred = poisson_rng(exp(lprobs) .* divisor);
}
