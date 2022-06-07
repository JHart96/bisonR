data {
  int num_nodes; // Number of nodes
  int num_fixed; // Number of fixed effects
  int num_random; // Number of random effects
  int num_random_groups; // Number of random effect groups
  vector[num_nodes] centrality_mu; // Means of Gaussian approximation of logit/log edge weights
  matrix[num_nodes, num_nodes] centrality_cov; // Covariance matrix of Gaussian approximation of logit/log edge weights
  matrix[num_nodes, num_fixed] design_fixed; // Design matrix for fixed effects
  matrix[num_nodes, num_random] design_random; // Design matrix for random effects
  array[num_random] int<lower=0, upper=num_random_groups> random_group_index; // Index for groupings for random effects
  real prior_fixed_mu;
  real<lower=0> prior_fixed_sigma;
  real prior_random_mean_mu;
  real<lower=0> prior_random_mean_sigma;
  real<lower=0> prior_random_std_sigma;
  real<lower=0> prior_error_sigma;
}
parameters {
  vector[num_fixed] beta_fixed; // Parameters for fixed effects.
  real<lower=0> sigma;
  vector[num_random] beta_random; // Parameters for random effects.
  vector[num_random_groups] random_group_mu; // Hyperpriors for random effects (mean).
  vector<lower=0>[num_random_groups] random_group_sigma; // Hyperpriors for random effects (std. dev.).
}
transformed parameters {
  vector[num_nodes] predictor;
  predictor = rep_vector(0, num_nodes);
  if (num_fixed > 0) predictor += design_fixed * beta_fixed;
  if (num_random > 0) predictor += design_random * beta_random;
}
model {
  centrality_mu ~ multi_normal(predictor, centrality_cov + diag_matrix(rep_vector(sigma, num_nodes)));
  beta_fixed ~ normal(prior_fixed_mu, prior_fixed_sigma);
  sigma ~ normal(0, prior_error_sigma);
  if (num_random > 0) {
    beta_random ~ normal(random_group_mu[random_group_index], random_group_sigma[random_group_index]);
    random_group_mu ~ normal(prior_random_mean_mu, prior_random_mean_sigma);
    random_group_sigma ~ normal(0, prior_random_std_sigma);
  }
}
generated quantities {
  vector[num_nodes] centrality_pred;
  centrality_pred = multi_normal_rng(predictor, centrality_cov + diag_matrix(rep_vector(sigma, num_nodes)));
}
