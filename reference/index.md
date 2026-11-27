# Package index

## Fitting edge weight models

Estimate the network, with uncertainty, from observation data.

- [`bison_model()`](https://bisonr.org/reference/bison_model.md) : Fit
  an edge weight model
- [`get_default_priors()`](https://bisonr.org/reference/get_default_priors.md)
  : Get the default priors for an edge model
- [`prior_check()`](https://bisonr.org/reference/prior_check.md) : Plot
  prior distributions
- [`prior_predictive_check()`](https://bisonr.org/reference/prior_predictive_check.md)
  : Prior predictive check

## Checking models

MCMC diagnostics, predictive checks and model comparison.

- [`plot_trace()`](https://bisonr.org/reference/plot_trace.md) : Trace
  plots of MCMC chains
- [`plot_predictions()`](https://bisonr.org/reference/plot_predictions.md)
  : Posterior predictive checks
- [`model_comparison()`](https://bisonr.org/reference/model_comparison.md)
  : Compare models with leave-one-out cross-validation

## Exploring the network

- [`summary(`*`<bison_model>`*`)`](https://bisonr.org/reference/summary.bison_model.md)
  : Summarise an edge weight model
- [`get_edgelist()`](https://bisonr.org/reference/get_edgelist.md) : Get
  the edge list with uncertainty
- [`get_effects()`](https://bisonr.org/reference/get_effects.md) : Get
  the estimates of other effects in an edge model
- [`draw_edgelist_samples()`](https://bisonr.org/reference/draw_edgelist_samples.md)
  : Draw samples of the edge list from the posterior
- [`plot_network()`](https://bisonr.org/reference/plot_network.md) :
  Plot the network with uncertainty
- [`bison_to_igraph()`](https://bisonr.org/reference/bison_to_igraph.md)
  : Convert a fitted edge model into a list of igraph networks

## Network metrics

- [`extract_metric()`](https://bisonr.org/reference/extract_metric.md) :
  Extract posterior samples of network metrics
- [`plot_metric()`](https://bisonr.org/reference/plot_metric.md) : Plot
  the posterior distribution of a network metric

## Regression

Use network quantities in brms models while keeping their uncertainty.

- [`bison_brm()`](https://bisonr.org/reference/bison_brm.md) : Fit a
  regression with a network metric
- [`bison_brm_get_prior()`](https://bisonr.org/reference/bison_brm_get_prior.md)
  : List the priors of a bison_brm() regression
- [`bison_mice()`](https://bisonr.org/reference/bison_mice.md) : Prepare
  network data for multiple imputation with the mice package

## Edge mixture models

Look for different types of social relationship.

- [`bison_mixture()`](https://bisonr.org/reference/bison_mixture.md) :
  Fit an edge mixture model
- [`summary(`*`<bison_mixture>`*`)`](https://bisonr.org/reference/summary.bison_mixture.md)
  : Summarise an edge mixture model
- [`get_network_component_probabilities()`](https://bisonr.org/reference/get_network_component_probabilities.md)
  : Get the probability of each number of components in the network
- [`get_edge_component_probabilities()`](https://bisonr.org/reference/get_edge_component_probabilities.md)
  : Get the probability of each edge belonging to each component
- [`get_component_means()`](https://bisonr.org/reference/get_component_means.md)
  : Get the mean edge weight of each mixture component
- [`plot_mixture()`](https://bisonr.org/reference/plot_mixture.md) :
  Plot the components of an edge mixture model

## Preparing data

- [`convert_gbi_to_bison()`](https://bisonr.org/reference/convert_gbi_to_bison.md)
  : Convert a group-by-individual matrix to bisonR format
- [`convert_duration_to_binary()`](https://bisonr.org/reference/convert_duration_to_binary.md)
  : Convert duration data to binary
- [`simulate_bison_model()`](https://bisonr.org/reference/simulate_bison_model.md)
  : Simulate observation data for an edge model
- [`simulate_bison_model_mixture()`](https://bisonr.org/reference/simulate_bison_model_mixture.md)
  : Simulate observation data with several types of edge

## Saving models

- [`bison_save()`](https://bisonr.org/reference/bison_save.md) : Save a
  fitted bisonR model to disk

- [`bison_load()`](https://bisonr.org/reference/bison_load.md) :

  Load a bisonR model saved with
  [`bison_save()`](https://bisonr.org/reference/bison_save.md)
