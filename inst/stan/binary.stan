// BISoN binary edge model.
// Generated from the same template as the count model; keep the two in sync.
data {
  int<lower=0> num_rows; // Number of data points
  int<lower=0> num_edges; // Number of edge weights to estimate
  int<lower=0> num_fixed; // Number of fixed effect parameters
  int<lower=0> num_random; // Number of random effect parameters
  int<lower=0> num_random_groups; // Number of random effect groups

  array[num_rows] int<lower=0> event; // Social events for each data point
  array[num_rows] int<lower=0> divisor; // Observation effort for each data point
  array[num_rows] int dyad_ids; // Dyad IDs of each observation for indexing edge weights
  matrix[num_rows, num_fixed] design_fixed; // Design matrix for fixed effects
  matrix[num_rows, num_random] design_random; // Design matrix for random effects.
  array[num_random] int<lower=0, upper=num_random_groups> random_group_index; // Index for groupings for random effects

  real prior_edge_mu; // Prior mean for edge weights
  real<lower=0> prior_edge_sigma; // Prior standard deviation for edge weights (or for their spread, with partial pooling)
  real prior_fixed_mu; // Prior mean for fixed effects
  real<lower=0> prior_fixed_sigma; // Prior standard deviation for fixed effects
  real prior_random_mean_mu; // Prior mean on centralisation of random effects
  real<lower=0> prior_random_mean_sigma; // Prior standard deviation on centralisation of random effects
  real<lower=0> prior_random_std_sigma; // Prior standard deviation on dispersion of random effects
  real<lower=0> prior_zero_prob_alpha; // Prior alpha on zero inflation
  real<lower=0> prior_zero_prob_beta; // Prior beta on zero inflation

  int<lower=0, upper=1> priors_only; // Whether to sample from only the priors
  int<lower=0, upper=1> partial_pooling; // Whether to pool edge weight estimates
  int<lower=0, upper=1> zero_inflated; // Whether to use zero-inflated edge model
}

parameters {
  vector[num_edges] edge_raw; // Standardised edge weights (non-centred parameterisation)
  vector[num_fixed] beta_fixed; // Parameters for fixed effects.
  vector[num_random] random_raw; // Standardised random effects (non-centred parameterisation)
  vector[num_random_groups] random_group_mu; // Hyperpriors for random effects (mean).
  vector<lower=0>[num_random_groups] random_group_sigma; // Hyperpriors for random effects (std. dev.).
  vector<lower=0>[partial_pooling] edge_sigma; // Spread of edge weights when pooling.
  vector<lower=0, upper=1>[zero_inflated] zero_prob; // Probability of a structural zero.
}

transformed parameters {
  vector[num_edges] edge_weight; // Edge weights on the link scale
  vector[num_random] beta_random; // Random effects
  vector[num_rows] predictor;

  if (partial_pooling == 1) {
    edge_weight = prior_edge_mu + edge_sigma[1] * edge_raw;
  } else {
    edge_weight = prior_edge_mu + prior_edge_sigma * edge_raw;
  }
  if (num_random > 0) {
    beta_random = random_group_mu[random_group_index] + random_group_sigma[random_group_index] .* random_raw;
  }

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
  if (!priors_only) {
    if (zero_inflated == 0) {
      event ~ binomial_logit(divisor, predictor);
    } else {
      for (i in 1:num_rows) {
        if (event[i] == 0) {
          target += log_mix(zero_prob[1], 0, binomial_logit_lpmf(0 | divisor[i], predictor[i]));
        } else {
          target += log1m(zero_prob[1]) + binomial_logit_lpmf(event[i] | divisor[i], predictor[i]);
        }
      }
    }
  }

  // Priors
  edge_raw ~ std_normal();
  if (partial_pooling == 1) {
    edge_sigma ~ normal(0, prior_edge_sigma);
  }
  if (zero_inflated == 1) {
    zero_prob ~ beta(prior_zero_prob_alpha, prior_zero_prob_beta);
  }
  if (num_fixed > 0) {
    beta_fixed ~ normal(prior_fixed_mu, prior_fixed_sigma);
  }
  if (num_random > 0) {
    random_raw ~ std_normal();
    random_group_mu ~ normal(prior_random_mean_mu, prior_random_mean_sigma);
    random_group_sigma ~ normal(0, prior_random_std_sigma);
  }
}

generated quantities {
  array[num_rows] int event_pred;
  vector[num_rows] log_lik;
  for (i in 1:num_rows) {
    if (zero_inflated == 1 && bernoulli_rng(zero_prob[1])) {
      event_pred[i] = 0;
    } else {
      event_pred[i] = binomial_rng(divisor[i], inv_logit(predictor[i]));
    }
    if (zero_inflated == 0) {
      log_lik[i] = binomial_logit_lpmf(event[i] | divisor[i], predictor[i]);
    } else if (event[i] == 0) {
      log_lik[i] = log_mix(zero_prob[1], 0, binomial_logit_lpmf(0 | divisor[i], predictor[i]));
    } else {
      log_lik[i] = log1m(zero_prob[1]) + binomial_logit_lpmf(event[i] | divisor[i], predictor[i]);
    }
  }
}
