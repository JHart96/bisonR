data {
  int<lower=0> num_rows; // Number of data points
  int<lower=0> num_nodes; // Number of nodes for multimembership effects
  int<lower=0> num_fixed; // Number of fixed effect parameters
  int<lower=0> num_random; // Number of random effect parameters
  int<lower=0> num_random_groups; // Number of random effect groups

  vector[num_rows] edge_mu; // Mean weight for each edge
  matrix[num_rows, num_rows] edge_cov; // Covariance matrix for random_group_indexaussian approximation of edge weights
  matrix[num_rows, num_fixed] design_fixed; // Design matrix for fixed effects
  matrix[num_rows, num_random] design_random; // Design matrix for random effects.
  array[num_rows] int node_ids_1; // Node 1 IDs for multimembership terms
  array[num_rows] int node_ids_2; // Node 2 IDs for multimembership terms
  int<lower=0, upper=1> include_multimembership; // Boolean variable for including multimembership effects
  array[num_random] int<lower=0, upper=num_random_groups> random_group_index; // Index for groupings for random effects

  real prior_fixed_mu;
  real<lower=0> prior_fixed_sigma;
  real prior_random_mean_mu;
  real<lower=0> prior_random_mean_sigma;
  real<lower=0> prior_random_std_sigma;
  real<lower=0> prior_error_sigma;
  real<lower=0> prior_multimembership_sigma;
}
parameters {
  vector[num_fixed] beta_fixed; // Parameters for fixed effects.
  real<lower=0> sigma;
  vector[include_multimembership ? num_nodes : 0] mm_nodes;
  array[include_multimembership] real<lower=0> sigma_mm;
  vector[num_random] beta_random; // Parameters for random effects.
  vector[num_random_groups] random_group_mu; // Hyperpriors for random effects (mean).
  vector<lower=0>[num_random_groups] random_group_sigma; // Hyperpriors for random effects (std. dev.).
}
transformed parameters {
  vector[num_rows] predictor;
  predictor =  design_fixed * beta_fixed;
  if (num_random > 0) {
    predictor += design_random * beta_random;
  }
  if (include_multimembership == 1) {
    predictor += mm_nodes[node_ids_1] + mm_nodes[node_ids_2];
  }
}
model {
  edge_mu ~ multi_normal(predictor, edge_cov + diag_matrix(rep_vector(sigma, num_rows)));
  beta_fixed ~ normal(prior_fixed_mu, prior_fixed_sigma);
  sigma ~ normal(0, prior_error_sigma);
  if (include_multimembership == 1) {
    mm_nodes ~ normal(0, sigma_mm);
    sigma_mm ~ normal(0, prior_multimembership_sigma);
  }
  if (num_random > 0) {
    beta_random ~ normal(random_group_mu[random_group_index], random_group_sigma[random_group_index]);
    random_group_mu ~ normal(prior_random_mean_mu, prior_random_mean_sigma);
    random_group_sigma ~ normal(0, prior_random_std_sigma);
  }
}
generated quantities {
  vector[num_rows] edge_pred;
  edge_pred = multi_normal_rng(predictor, edge_cov + diag_matrix(rep_vector(sigma, num_rows)));
}
