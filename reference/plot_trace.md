# Trace plots of MCMC chains

Plots the MCMC chains for some of the edge weights of a fitted model, to
check that the sampler worked. MCMC (Markov chain Monte Carlo) draws
samples from the posterior in several independent runs called chains.
The chains should overlap and show no trends; if they don't, the
posterior may not be reliable, and running more warmup iterations may
help.

## Usage

``` r
plot_trace(obj, par_ids = 1:5, ...)
```

## Arguments

- obj:

  A fitted `bison_model` object (fitted with Stan, not a conjugate
  model).

- par_ids:

  Indices of the edge weights to plot, in the order of
  [`get_edgelist()`](https://bisonr.org/reference/get_edgelist.md).
  Indices beyond the number of edges are ignored.

- ...:

  Additional arguments passed to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

`NULL`, invisibly. The function is called for its plots.

## See also

[`plot_predictions()`](https://bisonr.org/reference/plot_predictions.md)
for checking model fit.

## Examples

``` r
sim_data <- simulate_bison_model("binary", aggregated = TRUE)
fit_edge <- bison_model(
  (event | effort) ~ dyad(node_1_id, node_2_id),
  data = sim_data$df_sim, model_type = "binary",
  priors = get_default_priors("binary")
)
#> Running MCMC with 4 parallel chains...
#> 
#> Chain 1 finished in 0.2 seconds.
#> Chain 2 finished in 0.3 seconds.
#> Chain 3 finished in 0.2 seconds.
#> Chain 4 finished in 0.2 seconds.
#> 
#> All 4 chains finished successfully.
#> Mean chain execution time: 0.2 seconds.
#> Total execution time: 0.4 seconds.
#> 
plot_trace(fit_edge, par_ids = 1:4)
```
