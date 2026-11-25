# Sociogram plot with uncertainty of a fitted edge weight model object

Sociogram plot with uncertainty of a fitted edge weight model object

## Usage

``` r
plot_network(obj, ci = 0.9, lwd = 2, threshold = NULL)
```

## Arguments

- obj:

  An S3 edge model object to be summarised.

- ci:

  Credible interval to use in plot, based on quantiles.

- lwd:

  Line width scaling for edge weights

- threshold:

  Edge weights below this will be removed from the visualisation (for
  plotting purposes only). NULL if no threshold should be applied.
