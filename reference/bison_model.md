# Fit an edge model to data

Fit an edge model to data

## Usage

``` r
bison_model(
  formula,
  data,
  model_type = c("binary", "count"),
  directed = FALSE,
  partial_pooling = FALSE,
  zero_inflated = FALSE,
  priors = NULL,
  refresh = 0,
  mc_cores = 4,
  iter_sampling = 1000,
  iter_warmup = 1000,
  priors_only = FALSE
)
```

## Arguments

- formula:

  Formula specifying social events and sampling effort on the LHS and
  edge weights, fixed, and random effects on the RHS.

- data:

  Aggregated or disaggregated dataframe of dyadic observations.

- model_type:

  "binary" or "count", specifying the type of edge weight model to use.

- directed:

  `TRUE` or `FALSE` specifying whether the network is directed or not.

- partial_pooling:

  Whether to pool edge weights so that information is shared between
  edges.

- zero_inflated:

  Whether to use a zero-inflated model to model excess zeroes.

- priors:

  List of priors in the format supplied by
  [`get_default_priors()`](https://bisonr.org/reference/get_default_priors.md).

- refresh:

  Frequency of messages printed while running the sampler.

- mc_cores:

  Number of cores to use when running the sampler.

- iter_sampling:

  Number of iterations to use for posterior samples.

- iter_warmup:

  Number of iterations to use for warmup (will not be used for samples).

- priors_only:

  Whether to use priors as posteriors or to allow the posteriors to be
  updated by data.

## Value

An S3 edge model object containing edge samples and processed data.

## Details

Fits a BISoN edge weight model to a user-provided dataframe. The
function supports either aggregated (at the dyad-level) or disaggregated
(at the observation-level) dataframes. Node names or IDs need to be
formatted as factors with the same levels.

The type of edge model and the interpretation of edge weights used
depends on `model_type`, and will change the interpretation of the edge
weights.
