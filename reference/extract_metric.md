# Extract samples from posterior metrics for a fitted edge weight model

Extract samples from posterior metrics for a fitted edge weight model

## Usage

``` r
extract_metric(obj, metric_string, num_draws = 1000, standardise = FALSE)
```

## Arguments

- obj:

  Fitted S3 edge weight model.

- metric_string:

  Character name of metric function (e.g. node_strength).

- num_draws:

  Number of posterior draws to use.

- standardise:

  TRUE/FALSE indicating whether to mean-center the metrics within
  sample.

## Value

A matrix of metric samples where each column corresponds to a node, and
each row corresponds to a posterior draw.
