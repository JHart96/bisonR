# Save a fitted bisonR model to disk

Fitted models contain a CmdStan fit object whose draws live in temporary
CSV files, so [`saveRDS()`](https://rdrr.io/r/base/readRDS.html) on its
own does not produce a file that can be reloaded in a new R session.
`bison_save()` reads all draws into memory first and then saves the
complete object as a single `.rds` file. Use it to keep a fitted model
so it doesn't need to be refitted later.

## Usage

``` r
bison_save(bisonobj, file)
```

## Arguments

- bisonobj:

  A fitted model, such as the output of
  [`bison_model()`](https://bisonr.org/reference/bison_model.md).

- file:

  Path of the `.rds` file to write. `.rds` is appended if missing.

## Value

The path to the saved file, invisibly.

## See also

[`bison_load()`](https://bisonr.org/reference/bison_load.md)

## Examples

``` r
sim_data <- simulate_bison_model("binary", aggregated = TRUE)
fit_edge <- bison_model(
  (event | effort) ~ dyad(node_1_id, node_2_id),
  data = sim_data$df_sim,
  model_type = "binary_conjugate",
  priors = get_default_priors("binary_conjugate")
)
path <- tempfile(fileext = ".rds")
bison_save(fit_edge, path)
fit_reloaded <- bison_load(path)
```
