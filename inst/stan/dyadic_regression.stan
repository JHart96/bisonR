data {
  int<lower=0> N; // Number of data points
  int<lower=0> num_nodes; // Number of nodes for multimembership effects
  int<lower=0> K_fixed; // Number of fixed effect parameters
  int<lower=0> K_random; // Number of random effect parameters
  vector[N] edge_mu; // Mean weight for each edge
  matrix[N, N] edge_cov; // Covariance matrix for Gaussian approximation of edge weights
  matrix[N, K_fixed] X; // Design matrix for fixed effects
  array[N] int node_ids_1; // Node 1 IDs for multimembership terms
  array[N] int node_ids_2; // Node 2 IDs for multimembership terms
  int<lower=0, upper=1> include_multimembership; // Boolean variable for including multimembership effects
  int<lower=0> R; // Number of random effect groups
  matrix[N, K_random] Z; // Design matrix for random effects.
  array[K_random] int<lower=0, upper=R> G; // Index for groupings for random effects

}
parameters {
  vector[K_fixed] beta_fixed; // Parameters for fixed effects.
  real<lower=0> sigma;
  vector[include_multimembership ? num_nodes : 0] mm_nodes;
  array[include_multimembership] real<lower=0> sigma_mm;
  vector[K_random] beta_random; // Parameters for random effects.
  vector[R] H_mu; // Hyperpriors for random effects (mean).
  vector<lower=0>[R] H_sigma; // Hyperpriors for random effects (std. dev.).
}
transformed parameters {
  vector[N] predictor;
  predictor =  X * beta_fixed;
  if (K_random > 0) {
    predictor += Z * beta_random;
  }
  if (include_multimembership == 1) {
    predictor += mm_nodes[node_ids_1] + mm_nodes[node_ids_2];
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
  if (K_random > 0) {
    beta_random ~ normal(H_mu[G], H_sigma[G]);
    H_mu ~ normal(0, 2.5);
    H_sigma ~ normal(0, 1);
  }
}
generated quantities {
  vector[N] edge_pred;
  edge_pred = multi_normal_rng(predictor, edge_cov + diag_matrix(rep_vector(sigma, N)));
}
