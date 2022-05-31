data {
  int<lower=0> N; // Number of data points
  int<lower=0> K_fixed; // Number of fixed effect parameters
  int<lower=0> K_random; // Number of random effect parameters
  int<lower=0> R; // Number of random effect groups
  int y[N]; // Outcome for each data point (presence/absence)
  int divisor[N]; // Duration of each observation
  matrix[N, K_fixed] X; // Design matrix for fixed effects
  matrix[N, K_random] Z; // Design matrix for random effects.
  int<lower=0, upper=R> G[K_random]; // Index for groupings for random effects
}

parameters {
  vector[K_fixed] beta_fixed; // Parameters for fixed effects.
  vector[K_random] beta_random; // Parameters for random effects.
  vector[R] H_mu; // Hyperpriors for random effects (mean).
  vector<lower=0>[R] H_sigma; // Hyperpriors for random effects (std. dev.).
}

transformed parameters {
  vector[N] lprobs;
  lprobs = X * beta_fixed + Z * beta_random;
}

model {
  // Main model
  target += binomial_lpmf(y | divisor, inv_logit(lprobs));

  // Priors
  target += normal_lpdf(beta_fixed | 0, 1.0);
  target += normal_lpdf(beta_random | H_mu[G], H_sigma[G]);

  // Hyperpriors
  target += normal_lpdf(H_mu | 0, 1);
  target += exponential_lpdf(H_sigma | 1);
}

generated quantities {
  int y_pred[N];
  y_pred = binomial_rng(divisor, inv_logit(lprobs));
}
