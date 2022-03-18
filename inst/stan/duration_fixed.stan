data {
  int<lower=0> N; // Number of observations
  int<lower=0> M; // Number of dyads
  int<lower=0> dyad_ids[N]; // Indices for dyads.
  int<lower=0> K_fixed; // Number of fixed effect parameters
  real<lower=0, upper=1> y[N]; // Duration for each observation (> 0). Needs to be bounded above to length of observation, 1 for now but needs changing.
  int<lower=0> k[M]; // Number of events for each dyad. Needs to be changed to have M dyads or something.
  vector<lower=0>[M] d; // Number of observation periods for each dyad.
  matrix[N, K_fixed] X; // Design matrix for fixed effects
}

parameters {
  vector<lower=0>[M] lambda; // Parameters for dyad-level mean event rates.
  vector[K_fixed] beta_fixed; // Parameters for fixed effects.
}

transformed parameters {
  vector[N] log_pn = X * beta_fixed; // Calculate observation-level social preferences.
}

model {
  // Main model
  target += exponential_lpdf(y | lambda[dyad_ids] ./ exp(log_pn));
  target += poisson_lpmf(k | lambda .* d);

  // Priors
  target += exponential_lpdf(lambda | 1);
  target += normal_lpdf(beta_fixed | 0, 2.5);
}
/*
generated quantities {
  vector[N] mu = exp(log_pn) ./ lambda[dyad_ids]; // Calculate observation-level mean event times just in case.
}
*/
