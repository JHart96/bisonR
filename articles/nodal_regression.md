# Nodal regression

A nodal regression asks whether a trait of individuals, such as age, sex
or rank, is related to their position in the network, such as their
centrality. This article shows how to fit nodal regressions with
[`bison_brm()`](https://bisonr.org/reference/bison_brm.md), with the
node metric as the response or as a predictor. The metric is calculated
for many posterior draws of the edge model (each draw is one possible
network that is consistent with the data), so the uncertainty in each
node’s metric is included in the regression.

See
[`vignette("bison_brms")`](https://bisonr.org/articles/bison_brms.md)
for an overview of
[`bison_brm()`](https://bisonr.org/reference/bison_brm.md) and
[`vignette("network_metrics")`](https://bisonr.org/articles/network_metrics.md)
for the definitions of the metrics.

``` r

library(bisonR)
```

## Simulated data

We simulate 20 individuals whose ages range from 2 to 15 years. Each
individual has a *sociality*, which increases with age and also varies
between individuals of the same age. On the logit (log-odds) scale, the
edge weight of each dyad is the sum of the two individuals’ sociality,
plus an intercept and dyad-level variation. Older individuals therefore
tend to have stronger edges.

``` r

num_nodes <- 20
nodes <- paste0("ID", sprintf("%02d", 1:num_nodes))
df_nodal <- data.frame(
  node = factor(nodes, levels = nodes),
  age = round(runif(num_nodes, 2, 15))
)

# Sociality increases by 0.15 per year of age (logit scale)
sociality <- 0.15 * (df_nodal$age - mean(df_nodal$age)) + rnorm(num_nodes, 0, 0.3)

dyads <- t(combn(num_nodes, 2))
logit_weight <- -2 + sociality[dyads[, 1]] + sociality[dyads[, 2]] +
  rnorm(nrow(dyads), 0, 0.5)
effort <- sample(5:30, nrow(dyads), replace = TRUE)

df_obs <- data.frame(
  node_1_id = factor(nodes[dyads[, 1]], levels = nodes),
  node_2_id = factor(nodes[dyads[, 2]], levels = nodes),
  event = rbinom(nrow(dyads), effort, plogis(logit_weight)),
  effort = effort
)
head(df_obs)
#>   node_1_id node_2_id event effort
#> 1      ID01      ID02     1      6
#> 2      ID01      ID03     0      8
#> 3      ID01      ID04     3     28
#> 4      ID01      ID05     1     29
#> 5      ID01      ID06     0     10
#> 6      ID01      ID07     2     21
```

`df_obs` contains the observations, with one row per dyad. `event` is
the number of sampling periods in which the pair was seen together, and
`effort` is the number of sampling periods in which they were observed.
`df_nodal` has one row per individual.

## Fitting the edge model

The edge model is fitted in the same way as for any other analysis.
There are no observation-level effects, so the conjugate binary model
can be used. Conjugate models calculate the posterior directly, which is
fast, but they cannot include extra effects:

``` r

fit_edge <- bison_model(
  (event | effort) ~ dyad(node_1_id, node_2_id),
  data = df_obs,
  model_type = "binary_conjugate",
  priors = get_default_priors("binary_conjugate")
)
```

## The node data frame

A nodal regression needs a data frame with one row per node and a column
containing the node names. The names must match those used in the edge
model, which are the levels of the node columns in the observation data
(here `ID01` to `ID20`).
[`bison_brm()`](https://bisonr.org/reference/bison_brm.md) gives an
error if the data frame contains a node that is not in the edge model.
The data frame can contain only some of the nodes, for example if a
trait is missing for some individuals.

``` r

head(df_nodal)
#>   node age
#> 1 ID01   5
#> 2 ID02   7
#> 3 ID03   9
#> 4 ID04  14
#> 5 ID05   5
#> 6 ID06  14
all(df_nodal$node %in% levels(df_obs$node_1_id))
#> [1] TRUE
```

## Node metric as the response

To test whether age predicts node strength, write the metric inside
`bison()` on the left-hand side of the formula. The column inside
`node_strength()` is the column of the data frame that contains the node
names:

``` r

fit_strength <- bison_brm(
  bison(node_strength(node)) ~ age,
  fit_edge,
  df_nodal,
  num_draws = 5, # Small for demonstration; use 100 or more in practice
  refresh = 0,
  silent = 2
)
summary(fit_strength)
#>  Family: gaussian 
#>   Links: mu = identity 
#> Formula: bison_node_strength ~ age 
#>    Data: structure(list(data = structure(list(node = struct (Number of observations: 20) 
#>   Draws: 10 chains, each with iter = 2000; warmup = 1000; thin = 1;
#>          total post-warmup draws = 10000
#> 
#> Regression Coefficients:
#>           Estimate Est.Error l-95% CI u-95% CI
#> Intercept     1.50      0.61     0.32     2.71
#> age           0.21      0.06     0.09     0.34
#> 
#> Further Distributional Parameters:
#>       Estimate Est.Error l-95% CI u-95% CI
#> sigma     0.92      0.19     0.61     1.36
#> 
#> Draws were sampled using sample(hmc). Overall Rhat and ESS estimates
#> are not informative for brm_multiple models and are hence not displayed.
#> Please see ?brm_multiple for how to assess convergence of such models.
```

Node strength is the sum of a node’s edge weights on the natural scale
(here, probabilities), so it is the expected number of associates in a
sampling period. The `age` coefficient is the expected change in
strength per year of age. Its 95% credible interval (the range that
contains the coefficient with 95% posterior probability) is above zero,
which is consistent with the simulation. Convergence of the individual
models can be checked as described in
[`vignette("bison_brms")`](https://bisonr.org/articles/bison_brms.md).

### Standardising the metric

With `scale = "z"`, the metric is converted to z-scores within each
posterior draw: the mean across nodes is subtracted and the result is
divided by the standard deviation across nodes:

``` r

fit_strength_z <- bison_brm(
  bison(node_strength(node)) ~ age,
  fit_edge,
  df_nodal,
  num_draws = 5, # Small for demonstration; use 100 or more in practice
  scale = "z",
  refresh = 0,
  silent = 2
)
summary(fit_strength_z)
#>  Family: gaussian 
#>   Links: mu = identity 
#> Formula: bison_node_strength ~ age 
#>    Data: structure(list(data = structure(list(node = struct (Number of observations: 20) 
#>   Draws: 10 chains, each with iter = 2000; warmup = 1000; thin = 1;
#>          total post-warmup draws = 10000
#> 
#> Regression Coefficients:
#>           Estimate Est.Error l-95% CI u-95% CI
#> Intercept    -1.80      0.49    -2.71    -0.81
#> age           0.20      0.05     0.09     0.29
#> 
#> Further Distributional Parameters:
#>       Estimate Est.Error l-95% CI u-95% CI
#> sigma     0.75      0.16     0.49     1.11
#> 
#> Draws were sampled using sample(hmc). Overall Rhat and ESS estimates
#> are not informative for brm_multiple models and are hence not displayed.
#> Please see ?brm_multiple for how to assess convergence of such models.
```

The `age` coefficient is now the change in strength, in standard
deviations, per year of age. Standardising makes coefficients easier to
compare between metrics. It also makes priors easier to set, because the
scale of the response is known in advance. Standardising also removes
differences in the overall level of the metric between posterior draws,
so the model describes the relative position of nodes within each
network rather than their absolute values. Which of these is appropriate
depends on the question.

## Node metric as a predictor

A node metric can also be used as a predictor, for example to test
whether more social individuals have higher values of some outcome. As
an illustration, age is used as the response here:

``` r

fit_predictor <- bison_brm(
  age ~ bison(node_strength(node)),
  fit_edge,
  df_nodal,
  num_draws = 5, # Small for demonstration; use 100 or more in practice
  refresh = 0,
  silent = 2
)
summary(fit_predictor)
#>  Family: gaussian 
#>   Links: mu = identity 
#> Formula: age ~ bison_node_strength 
#>    Data: structure(list(data = structure(list(node = struct (Number of observations: 20) 
#>   Draws: 10 chains, each with iter = 2000; warmup = 1000; thin = 1;
#>          total post-warmup draws = 10000
#> 
#> Regression Coefficients:
#>                     Estimate Est.Error l-95% CI u-95% CI
#> Intercept               1.74      2.24    -2.60     6.23
#> bison_node_strength     2.01      0.57     0.84     3.12
#> 
#> Further Distributional Parameters:
#>       Estimate Est.Error l-95% CI u-95% CI
#> sigma     2.98      0.54     2.14     4.24
#> 
#> Draws were sampled using sample(hmc). Overall Rhat and ESS estimates
#> are not informative for brm_multiple models and are hence not displayed.
#> Please see ?brm_multiple for how to assess convergence of such models.
```

The `bison_node_strength` coefficient is the expected change in the
response per unit of strength. Each posterior draw of the metric is
treated as a possible value of the predictor, a separate model is fitted
to each, and the results are combined. Uncertainty in the metric
therefore widens the credible intervals.

## Choosing a metric

The metric should follow from the biological question. Some common
choices are:

- `node_strength`: how much an individual associates in total. Strength
  depends on the number of other individuals in the network, so it is
  not directly comparable between networks of different sizes.
- `node_eigen`: eigenvector centrality, which is high when an individual
  has strong connections to individuals that are themselves well
  connected.
- `node_betweenness` and `node_closeness`: an individual’s position on
  paths through the network, which may be relevant to questions about
  transmission.
- Thresholded metrics such as `node_degree[0.3]`: the number of edges
  with a weight at or above the threshold (here 0.3), calculated in each
  posterior draw.

Thresholded metrics are written in the same way, with the threshold in
square brackets:

``` r

fit_degree <- bison_brm(
  bison(node_degree[0.3](node)) ~ age,
  fit_edge,
  df_nodal,
  num_draws = 5, # Small for demonstration; use 100 or more in practice
  family = poisson(),
  refresh = 0,
  silent = 2
)
summary(fit_degree)
#>  Family: poisson 
#>   Links: mu = log 
#> Formula: bison_node_degree_0_3 ~ age 
#>    Data: structure(list(data = structure(list(node = struct (Number of observations: 20) 
#>   Draws: 10 chains, each with iter = 2000; warmup = 1000; thin = 1;
#>          total post-warmup draws = 10000
#> 
#> Regression Coefficients:
#>           Estimate Est.Error l-95% CI u-95% CI
#> Intercept    -0.21      0.47    -1.18     0.65
#> age           0.14      0.04     0.07     0.23
#> 
#> Draws were sampled using sample(hmc). Overall Rhat and ESS estimates
#> are not informative for brm_multiple models and are hence not displayed.
#> Please see ?brm_multiple for how to assess convergence of such models.
```

Degree is a count, so a Poisson family is used here, and the `age`
coefficient is on the log scale. The choice of threshold is arbitrary
and results can change with it, so check whether they are sensitive to
the threshold (see
[`vignette("network_metrics")`](https://bisonr.org/articles/network_metrics.md)).

## Priors

Priors describe which parameter values are plausible before the data are
seen. By default, brms uses flat priors on regression coefficients,
which treat all values as equally plausible.
[`bison_brm_get_prior()`](https://bisonr.org/reference/bison_brm_get_prior.md)
takes the same formula, edge model and data as
[`bison_brm()`](https://bisonr.org/reference/bison_brm.md) and returns
the default priors:

``` r

prior <- bison_brm_get_prior(
  bison(node_strength(node)) ~ age,
  fit_edge,
  df_nodal
)
prior
#>                   prior     class coef group resp dpar nlpar lb ub tag
#>                  (flat)         b                                     
#>                  (flat)         b  age                                
#>  student_t(3, 3.3, 2.5) Intercept                                     
#>    student_t(3, 0, 2.5)     sigma                             0       
#>        source
#>       default
#>  (vectorized)
#>       default
#>       default
```

The result is a `brmsprior` object that lists the classes and
coefficient names that priors can be set for. Priors are then set with
[`brms::set_prior()`](https://paulbuerkner.com/brms/reference/set_prior.html)
(or by modifying this object) and passed to
[`bison_brm()`](https://bisonr.org/reference/bison_brm.md) with the
`prior` argument. Here a normal prior with standard deviation 1 is
placed on the `age` coefficient:

``` r

prior <- brms::set_prior("normal(0, 1)", class = "b", coef = "age")

fit_strength_prior <- bison_brm(
  bison(node_strength(node)) ~ age,
  fit_edge,
  df_nodal,
  num_draws = 5, # Small for demonstration; use 100 or more in practice
  prior = prior,
  refresh = 0,
  silent = 2
)
brms::prior_summary(fit_strength_prior)
#>                   prior     class coef group resp dpar nlpar lb ub tag  source
#>                  (flat)         b                                      default
#>            normal(0, 1)         b  age                                    user
#>  student_t(3, 3.5, 2.5) Intercept                                      default
#>    student_t(3, 0, 2.5)     sigma                             0        default
```

Priors should be chosen on the scale of the response. For example, a
sensible prior on the effect of age on strength depends on the units of
age and on the range of strength values. This is one reason that
`scale = "z"` can make priors easier to set.
