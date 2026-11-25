# Generates a summary object for an edge model object

Generates a summary object for an edge model object

## Usage

``` r
# S3 method for class 'bison_model'
summary(object, ci = 0.9, transform = TRUE, ...)
```

## Arguments

- object:

  An S3 edge model object to be summarised.

- ci:

  Credible interval to use in summary, based on quantiles.

- transform:

  `TRUE` or `FALSE` specifying whether to transform the edge weights
  from the internal link function scale.

- ...:

  Additional arguments
