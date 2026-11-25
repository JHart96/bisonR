# Draw samples from edgelist posterior

Draw samples from edgelist posterior

## Usage

``` r
draw_edgelist_samples(obj, num_draws)
```

## Arguments

- obj:

  An S3 edge model object to be summarised.

- num_draws:

  Number of sample distributions to draw from the posterior.

## Value

A `data.frame` of sample draws from the posteriors, where each column
corresponds to a posterior draw.
