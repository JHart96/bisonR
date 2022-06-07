data {
  int<lower=0> N; // Number of observations
  int<lower=0> M; // Number of dyads
  array[N] int<lower=0> dyad_ids; // Indices for dyads.
  int<lower=0> K_fixed; // Number of fixed effect parameters
  array[N] real<lower=0> y; // Duration for each observation (> 0). Needs to be bounded above to length of observation, 1 for now but needs changing.
  array[M] int<lower=0> k; // Number of events for each dyad. Needs to be changed to have M dyads or something.
  vector<lower=0>[M] divisor; // Number of observation periods for each dyad.
  matrix[N, K_fixed] X; // Design matrix for fixed effects
  int<lower=0> K_random; // Number of random effect parameters
  int<lower=0> R; // Number of random effect groups
  matrix[N, K_random] Z; // Design matrix for random effects.
  array[K_random] int<lower=0, upper=R> G; // Index for groupings for random effects
  real prior_fixed_mu;
  real<lower=0> prior_fixed_sigma;
  real prior_random_mean_mu;
  real<lower=0> prior_random_mean_sigma;
  real<lower=0> prior_random_std_sigma;
}

parameters {
  vector<lower=0>[M] lambda; // Parameters for dyad-level mean event rates.
  vector[K_fixed] beta_fixed; // Parameters for fixed effects.
  vector[K_random] beta_random; // Parameters for random effects.
  vector[R] H_mu; // Hyperpriors for random effects (mean).
  vector<lower=0>[R] H_sigma; // Hyperpriors for random effects (std. dev.).
}

transformed parameters {
  vector[N] log_pn = X * beta_fixed; // Calculate observation-level social preferences.
  if (K_random > 0) {
    log_pn += Z * beta_random;
  }
}

model {
  // Main model
  target += exponential_lpdf(y | lambda[dyad_ids] ./ exp(log_pn));
  target += poisson_lpmf(k | lambda .* divisor);

  // Priors
  target += exponential_lpdf(lambda | 1);
  target += normal_lpdf(beta_fixed | prior_fixed_mu, prior_fixed_sigma);

  if (K_random > 0) {
    target += normal_lpdf(beta_random | H_mu[G], H_sigma[G]);
    // Hyperpriors
    target += normal_lpdf(H_mu | prior_random_mean_mu, prior_random_mean_sigma);
    target += normal_lpdf(H_sigma | 0, prior_random_std_sigma);
  }
}

generated quantities {
  vector[N] mu = exp(log_pn) ./ lambda[dyad_ids]; // Calculate observation-level mean event times just in case.
}

