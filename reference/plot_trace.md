# Traceplot of MCMC chains

Traceplot of MCMC chains

## Usage

``` r
plot_trace(obj, par_ids = 1:5, ...)
```

## Arguments

- obj:

  A fitted S3 object for edge weight, dyadic regression, or nodal
  regression models.

- par_ids:

  A vector of parameter IDs to be plotted.

- ...:

  Additional arguments to be passed to base R plots.

## Details

Plots MCMC chains for a fitted model. If the model has too many
parameters to plot, only the first 12 will be plotted. To change the
parameters being plotted, pass the `par_ids` argument
