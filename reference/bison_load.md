# Load a bisonR model saved with `bison_save()`

Reads a fitted model saved with
[`bison_save()`](https://bisonr.org/reference/bison_save.md), so it can
be used in a new R session in the same way as a newly fitted model.
Models saved with earlier versions of bisonR (0.4.x), which kept the
Stan fit in a separate file, are also read.

## Usage

``` r
bison_load(file)
```

## Arguments

- file:

  Path of an `.rds` file written by
  [`bison_save()`](https://bisonr.org/reference/bison_save.md). `.rds`
  is appended if the file isn't found without it.

## Value

The fitted model object.

## See also

[`bison_save()`](https://bisonr.org/reference/bison_save.md)

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
