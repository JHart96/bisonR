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
