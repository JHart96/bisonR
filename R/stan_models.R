require(cmdstanr)

build_stan_model <- function(model_code) {
  model_filepath <- cmdstanr::write_stan_file(model_code)
  model <- cmdstanr::cmdstan_model(model_filepath)
}

stan_model_binary_fixed_code <- "
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

stan_model_binary_mixed_code <- "
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
"

stan_model_count_fixed_code <- "
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
"

stan_model_count_mixed_code <- "
data {
  int<lower=0> N; // Number of data points
  int<lower=0> K_fixed; // Number of fixed effect parameters
  int<lower=0> K_random; // Number of random effect parameters
  int<lower=0> R; // Number of random effect groups
  int y[N]; // Outcome for each data point (presence/absence)
  vector[N] divisor; // Duration of each observation
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
  target += poisson_lpmf(y | exp(lprobs) .* divisor);

  // Priors
  target += normal_lpdf(beta_fixed | 0, 2.5);
  target += normal_lpdf(beta_random | H_mu[G], H_sigma[G]);

  // Hyperpriors
  target += normal_lpdf(H_mu | 0, 2.5);
  target += exponential_lpdf(H_sigma | 1);
}

generated quantities {
  int y_pred[N];
  y_pred = poisson_rng(exp(lprobs) .* divisor);
}
"

stan_model_duration_fixed_code <- "
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
"

stan_model_dyadic_regression_code <- "
data {
  int<lower=0> N; // Number of data points
  int<lower=0> num_nodes; // Number of nodes for multimembership effects
  int<lower=0> K_fixed; // Number of fixed effect parameters
  vector[N] edge_mu; // Mean weight for each edge
  matrix[N, N] edge_cov; // Covariance matrix for Gaussian approximation of edge weights
  matrix[N, K_fixed] X; // Design matrix for fixed effects
  int node_ids_1[N]; // Node 1 IDs for multimembership terms
  int node_ids_2[N]; // Node 2 IDs for multimembership terms
  int<lower=0, upper=1> include_multimembership; // Boolean variable for including multimembership effects
}

parameters {
  vector[K_fixed] beta_fixed; // Parameters for fixed effects.
  real<lower=0> sigma;
  vector[include_multimembership ? num_nodes : 0] mm_nodes;
  real<lower=0> sigma_mm[include_multimembership];
}

transformed parameters {
  vector[N] predictor;
  predictor =  X * beta_fixed;
  if (include_multimembership == 1) {
    predictor = predictor + mm_nodes[node_ids_1] + mm_nodes[node_ids_2];
  }
}

model {
  // Main model
  edge_mu ~ multi_normal(predictor, edge_cov + diag_matrix(rep_vector(sigma, N)));

  // Priors
  beta_fixed ~ normal(0, 1);
  sigma ~ normal(0, 1);

  if (include_multimembership == 1) {
    mm_nodes ~ normal(0, sigma_mm);
    sigma_mm ~ normal(0, 1);
  }
}

generated quantities {
  vector[N] edge_pred;
  edge_pred = multi_normal_rng(predictor, edge_cov + diag_matrix(rep_vector(sigma, N)));
}
"

stan_model_nodal_regression_code <- "
data {
  int num_nodes; // Number of nodes
  int num_fixed; // Number of fixed effects
  vector[num_nodes] centrality_mu; // Means of Gaussian approximation of logit/log edge weights
  matrix[num_nodes, num_nodes] centrality_cov; // Covariance matrix of Gaussian approximation of logit/log edge weights
  matrix[num_nodes, num_fixed] design_fixed; // Design matrix for fixed effects
}

parameters {
  vector[num_fixed] beta_fixed; // Parameters for fixed effects.
  real<lower=0> sigma;
}

transformed parameters {
  vector[num_nodes] predictor;
  predictor = design_fixed * beta_fixed;
}

model {
  centrality_mu ~ multi_normal(predictor, centrality_cov + diag_matrix(rep_vector(sigma, num_nodes)));
  beta_fixed ~ normal(0, 1);
  sigma ~ normal(0, 1);
}

generated quantities {
  vector[num_nodes] centrality_pred;
  centrality_pred = multi_normal_rng(predictor, centrality_cov + diag_matrix(rep_vector(sigma, num_nodes)));
}
"
