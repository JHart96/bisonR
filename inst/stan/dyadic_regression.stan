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
