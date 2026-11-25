# Retrieves an edgelist with uncertainty for a fitted edge weight model object

Retrieves an edgelist with uncertainty for a fitted edge weight model
object

## Usage

``` r
get_edgelist(obj, ci = 0.9, transform = TRUE)
```

## Arguments

- obj:

  An S3 edge model object.

- ci:

  Credible interval to use in summary, based on quantiles.

- transform:

  `TRUE` or `FALSE` specifying whether to transform the edge weights
  from the internal link function scale.

## Value

A `data.frame` object with columns of node IDs, median, lower, and upper
bounds.
