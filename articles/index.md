# Articles

### Getting started

Start here for a first analysis and the data format used by bisonR.

- [Getting started](https://bisonr.org/articles/getting_started.md):

  A first analysis: fitting an edge weight model, checking it, plotting
  the network, and using it in regressions and network metrics.

- [Preparing your data](https://bisonr.org/articles/preparing_data.md):

  The data format used by bisonR: events and observation effort,
  aggregation, directed data, unobserved dyads and group-by-individual
  matrices.

### Fitting models

Priors, diagnostics and additional effects in edge weight models.

- [Choosing and checking priors](https://bisonr.org/articles/priors.md):

  What each prior controls, the link scales, and checking priors with
  prior_check() and prior predictive checks.

- [Model diagnostics](https://bisonr.org/articles/diagnostics.md):

  MCMC diagnostics and posterior predictive checks, with a worked
  example of excess zeros.

- [Edge models with additional
  effects](https://bisonr.org/articles/edge_model_effects.md):

  Fixed effects, random effects and partial pooling in edge weight
  models, and how they change the interpretation of edge weights.

### Types of data

Worked analyses for directed networks, group-based observations and
count data.

- [Directed networks](https://bisonr.org/articles/directed_networks.md):

  Fitting, plotting and analysing directed networks, including in- and
  out-strength and dyadic regression.

- [Group-based observations](https://bisonr.org/articles/group_data.md):

  Gambit-of-the-group data: converting group-by-individual matrices,
  group-level effects and pairs that could not have been seen together.

- [Count data and continuous observation
  effort](https://bisonr.org/articles/count_data.md):

  Counts of interactions with continuous observation time, such as focal
  follows, and the interpretation of edge weights as rates.

### Analyses

Network metrics, regression, comparing networks, mixture models and
working with posterior networks.

- [Network metrics](https://bisonr.org/articles/network_metrics.md):

  Definitions of the node, edge and global metrics available in bisonR,
  and how to extract their posterior distributions.

- [Regression with brms](https://bisonr.org/articles/bison_brms.md):

  Dyadic, nodal and global regressions with bison_brm(), including
  formula syntax, scaling and priors.

- [Nodal regression](https://bisonr.org/articles/nodal_regression.md):

  A worked nodal regression relating an individual trait to centrality,
  with the metric as response or predictor.

- [Comparing
  networks](https://bisonr.org/articles/comparing_networks.md):

  Comparing two or more networks using global metrics, global
  regressions and edges present in both networks.

- [Edge mixture models](https://bisonr.org/articles/mixture_models.md):

  Looking for different types of social relationship with edge mixture
  models.

- [Working with posterior
  networks](https://bisonr.org/articles/posterior_networks.md):

  Using posterior draws directly: igraph networks, metrics not provided
  by bisonR, community detection, edge contrasts and saving models.

### Help

- [Installing cmdstanr on
  Windows](https://bisonr.org/articles/windows_installation.md):

  Solutions to common problems when installing cmdstanr and CmdStan on
  Windows.
