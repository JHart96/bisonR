data {
  int<lower=0> num_rows; // Number of data points
  int<lower=0> num_edges; // Number of edge weights to estimate
  int<lower=0> num_fixed; // Number of fixed effect parameters
  int<lower=0> num_random; // Number of random effect parameters
  int<lower=0> num_random_groups; // Number of random effect groups

  array[num_rows] int event; // Outcome for each data point (presence/absence)
  array[num_rows] int divisor; // Duration of each observation
  array[num_rows] int dyad_ids; // Dyad IDs of each observation for indexing edge weights
  matrix[num_rows, num_fixed] design_fixed; // Design matrix for fixed effects
  matrix[num_rows, num_random] design_random; // Design matrix for random effects.
  array[num_random] int<lower=0, upper=num_random_groups> random_group_index; // Index for groupings for random effects

  real prior_edge_mu; // Prior mean for fixed effects
  real<lower=0> prior_edge_sigma; // Prior standard deviation for fixed effects
  real prior_fixed_mu; // Prior mean for fixed effects
  real<lower=0> prior_fixed_sigma; // Prior standard deviation for fixed effects
  real prior_random_mean_mu; // Prior mean on centralisation of random effects
  real<lower=0> prior_random_mean_sigma; // Prior standard deviation on centralisation of random effects
  real<lower=0> prior_random_std_sigma; // Prior standard deviation on dispersion of random effects
}

parameters {
  vector[num_edges] edge_weight; // Parameters for edge weights.
  vector[num_fixed] beta_fixed; // Parameters for fixed effects.
  vector[num_random] beta_random; // Parameters for random effects.
  vector[num_random_groups] random_group_mu; // Hyperpriors for random effects (mean).
  vector<lower=0>[num_random_groups] random_group_sigma; // Hyperpriors for random effects (std. dev.).
}

transformed parameters {
  vector[num_rows] predictor;
  predictor = rep_vector(0, num_rows);
  if (num_edges > 0) {
    predictor += edge_weight[dyad_ids];
  }
  if (num_fixed > 0) {
    predictor += design_fixed * beta_fixed;
  }
  if (num_random > 0) {
    predictor += design_random * beta_random;
  }
}

model {
  // Main model
  event ~ binomial(divisor, inv_logit(predictor));

  // Priors
  if (num_edges > 0) {
    edge_weight ~ normal(prior_edge_mu, prior_edge_sigma);
  }

  if (num_fixed > 0) {
    beta_fixed ~ normal(prior_fixed_mu, prior_fixed_sigma);
  }

  if (num_random > 0) {
    beta_random ~ normal(random_group_mu[random_group_index], random_group_sigma[random_group_index]);
    // Hyperpriors
    random_group_mu ~ normal(prior_random_mean_mu, prior_random_mean_sigma);
    random_group_sigma ~ normal(0, prior_random_std_sigma);
  }
}

generated quantities {
  array[num_rows] int event_pred;
  event_pred = binomial_rng(divisor, inv_logit(predictor));
}
