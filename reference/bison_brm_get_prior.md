# List the priors of a bison_brm() regression

Shows the parameters of a
[`bison_brm()`](https://bisonr.org/reference/bison_brm.md) regression
that can be given priors, and their default priors, using
[`brms::get_prior()`](https://paulbuerkner.com/brms/reference/default_prior.html).
Use it before fitting to choose priors, which can then be set with the
`prior` argument of
[`bison_brm()`](https://bisonr.org/reference/bison_brm.md). A prior is
the distribution that describes what is plausible for a parameter before
seeing the data.

## Usage

``` r
bison_brm_get_prior(
  formula,
  edgemodel_list,
  data_list,
  scale = c("none", "centre", "z"),
  subset = NULL,
  z_score = FALSE,
  ...
)
```

## Arguments

- formula:

  A brms formula containing one `bison()` term, as for
  [`bison_brm()`](https://bisonr.org/reference/bison_brm.md).

- edgemodel_list:

  A fitted `bison_model` object, or a list of them.

- data_list:

  A data frame of regression variables, or a list of them, as for
  [`bison_brm()`](https://bisonr.org/reference/bison_brm.md).

- scale, subset, z_score:

  As for [`bison_brm()`](https://bisonr.org/reference/bison_brm.md). The
  scale of the network variable affects the default priors, so use the
  same values as in
  [`bison_brm()`](https://bisonr.org/reference/bison_brm.md).

- ...:

  Additional arguments passed to
  [`brms::get_prior()`](https://paulbuerkner.com/brms/reference/default_prior.html),
  such as `family`.

## Value

A data frame of priors (a `brmsprior` object) from
[`brms::get_prior()`](https://paulbuerkner.com/brms/reference/default_prior.html),
with one row per parameter or class of parameters.

## See also

[`bison_brm()`](https://bisonr.org/reference/bison_brm.md),
[`brms::set_prior()`](https://paulbuerkner.com/brms/reference/set_prior.html)

## Examples

``` r
sim_data <- simulate_bison_model("binary", aggregated = TRUE)
df <- sim_data$df_sim
fit_edge <- bison_model(
  (event | effort) ~ dyad(node_1_id, node_2_id),
  data = df,
  model_type = "binary_conjugate",
  priors = get_default_priors("binary_conjugate")
)
bison_brm_get_prior(
  bison(edge_weight(node_1_id, node_2_id)) ~ age_diff + (1 | mm(node_1_id, node_2_id)),
  fit_edge,
  df
)
#>                    prior     class      coef                group resp dpar
#>  student_t(3, -0.5, 2.5) Intercept                                         
#>                   (flat)         b                                         
#>                   (flat)         b  age_diff                               
#>     student_t(3, 0, 2.5)        sd                                         
#>     student_t(3, 0, 2.5)        sd           mmnode_1_idnode_2_id          
#>     student_t(3, 0, 2.5)        sd Intercept mmnode_1_idnode_2_id          
#>     student_t(3, 0, 2.5)     sigma                                         
#>  nlpar lb ub tag       source
#>                       default
#>                       default
#>                  (vectorized)
#>         0             default
#>         0        (vectorized)
#>         0        (vectorized)
#>         0             default
```
