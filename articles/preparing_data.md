# Preparing your data

``` r

library(bisonR)
```

This article describes the data format used by bisonR and how to convert
common types of data to it.

## The basic format

[`bison_model()`](https://bisonr.org/reference/bison_model.md) takes a
data frame in which each row contains observations of one pair of
individuals (a dyad). This is similar to an edge list, but contains the
observation data from which edge weights are estimated rather than the
edge weights themselves:

| Column | Meaning |
|----|----|
| Node 1, node 2 | The two individuals in the dyad. Names or IDs. |
| Event | How many social events were seen: 0/1 for binary data, or a count. |
| Effort | How much the dyad was observed: number of sampling periods, hours of observation, etc. |
| Anything else | Optional predictors, such as location, observer or individual traits. |

``` r

df <- data.frame(
  node_1 = c("Ada", "Ada", "Ada", "Bea", "Bea", "Cat"),
  node_2 = c("Bea", "Cat", "Dot", "Cat", "Dot", "Dot"),
  together = c(8, 1, 0, 5, 2, 9),
  periods = c(10, 10, 4, 10, 4, 12)
)
df
#>   node_1 node_2 together periods
#> 1    Ada    Bea        8      10
#> 2    Ada    Cat        1      10
#> 3    Ada    Dot        0       4
#> 4    Bea    Cat        5      10
#> 5    Bea    Dot        2       4
#> 6    Cat    Dot        9      12
```

Here `together` is the number of sampling periods in which each pair
were seen together, and `periods` is the number of sampling periods in
which they could have been seen together. These columns are given on the
left-hand side of the formula:

``` r

fit <- bison_model(
  (together | periods) ~ dyad(node_1, node_2),
  data = df,
  model_type = "binary_conjugate",
  priors = get_default_priors("binary_conjugate")
)
summary(fit)
#> === Fitted BISoN edge model ===
#> Data type: binary
#> Formula: (together | periods) ~ dyad(node_1, node_2)
#> Number of nodes: 4
#> Number of dyads: 6
#> Directed: FALSE
#> === Edge list summary ===
#>             median    5%   95%
#> Ada <-> Bea  0.765 0.528 0.922
#> Ada <-> Cat  0.146 0.032 0.369
#> Ada <-> Dot  0.132 0.011 0.454
#> Bea <-> Cat  0.501 0.275 0.725
#> Bea <-> Dot  0.503 0.187 0.808
#> Cat <-> Dot  0.719 0.506 0.886
```

Ada and Dot were observed in only 4 periods, so their edge weight is
less certain than that of Ada and Cat (10 periods), although both pairs
have low point estimates.

### What “effort” means

Effort is how much each dyad was observed, for example the number of
scans that included both individuals, or the number of hours they were
both under observation. It tells the model how much information there is
about each edge weight: the more a dyad was observed, the more certain
its edge weight.

If every row has the same effort, it can be written as a number instead
of a column: `(event | 1) ~ dyad(node_1, node_2)`.

### Node names

Node names can be characters, numbers or factors. If both node columns
are factors with the same levels, every level is treated as a node,
including any that do not appear in the data. Otherwise, the nodes are
the names that appear in the data, plus the levels of a node column that
is a factor.

## Binary or count data?

- Use a **binary** model (`model_type = "binary"`) when each sampling
  period records whether a social event happened or not, such as being
  seen in the same group during a scan. The edge weight is the
  probability of an event in a sampling period. Events can’t exceed
  effort.
- Use a **count** model (`model_type = "count"`) when you count
  interactions, such as the number of grooming bouts during focal
  follows. The edge weight is the rate of events per unit effort (for
  example, per hour). Effort doesn’t need to be a whole number.

``` r

df_counts <- data.frame(
  node_1 = c("Ada", "Ada", "Bea"),
  node_2 = c("Bea", "Cat", "Cat"),
  grooming_bouts = c(12, 0, 3),
  hours_observed = c(5.5, 2, 6.25)
)
fit_counts <- bison_model(
  (grooming_bouts | hours_observed) ~ dyad(node_1, node_2),
  data = df_counts,
  model_type = "count_conjugate",
  priors = get_default_priors("count_conjugate")
)
get_edgelist(fit_counts)
#>   node_1 node_2 median    5%   95%
#> 1    Ada    Bea  1.923 1.169 2.966
#> 2    Ada    Cat  0.235 0.017 1.033
#> 3    Bea    Cat  0.498 0.180 1.065
```

bisonR doesn’t yet have a model for the length of interactions (for
example, minutes spent grooming). For now,
[`convert_duration_to_binary()`](https://bisonr.org/reference/convert_duration_to_binary.md)
can turn this kind of data into binary data by splitting time into
fixed-length periods. This is an approximation and should be used with
caution.

## Aggregated or disaggregated?

Data can have one row per observation (disaggregated) or one row per
dyad, with events and effort summed (aggregated). Both give the same
edge weights, but differ in what else can be included in the model:

- Observation-level predictors, such as the location or time of day of
  each sampling period, need disaggregated data, since they vary between
  rows of the same dyad.
- Aggregated data are much faster to fit with Stan. (Conjugate models
  give the same results for aggregated and disaggregated data.)

Data can be aggregated as far as possible while keeping the predictors
that are needed. For example, if location is the only observation-level
predictor, the data can be aggregated to one row per dyad per location:

``` r

sim <- simulate_bison_model("binary", aggregated = FALSE)
df_obs <- sim$df_sim
nrow(df_obs)
#> [1] 135

df_by_location <- aggregate(cbind(event, effort) ~ node_1_id + node_2_id + location, data = df_obs, FUN = sum)
nrow(df_by_location)
#> [1] 96
head(df_by_location)
#>   node_1_id node_2_id location event effort
#> 1         2         3        1     0      2
#> 2         2         4        1     2      3
#> 3         3         4        1     0      2
#> 4         3         5        1     1      2
#> 5         2         6        1     0      1
#> 6         3         6        1     0      1
```

## Undirected and directed networks

By default, networks are **undirected**: the edge between A and B is the
same as the edge between B and A, and rows can list the pair in either
order.

``` r

df_undirected <- data.frame(
  node_1 = c("Ada", "Bea"), node_2 = c("Bea", "Ada"), event = c(3, 2), effort = c(5, 5)
)
fit_undirected <- bison_model(
  (event | effort) ~ dyad(node_1, node_2), data = df_undirected,
  model_type = "binary_conjugate", priors = get_default_priors("binary_conjugate")
)
#> 1 dyad(s) have rows in both orders (e.g. A, B and B, A). Because directed = FALSE, both orders are combined into one edge. If node_1 and node_2 are an actor and a receiver, use directed = TRUE.
get_edgelist(fit_undirected) # One edge, using all 5 events in 10 periods
#>   node_1 node_2 median   5%   95%
#> 1    Ada    Bea  0.501 0.27 0.733
```

For behaviours with an actor and a receiver, such as grooming or
aggression, use `directed = TRUE`. Then `node_1` is the actor and
`node_2` the receiver, and A → B and B → A are separate edges with their
own weights:

``` r

fit_directed <- bison_model(
  (event | effort) ~ dyad(node_1, node_2), data = df_undirected, directed = TRUE,
  model_type = "binary_conjugate", priors = get_default_priors("binary_conjugate")
)
get_edgelist(fit_directed)
#>   node_1 node_2 median    5%   95%
#> 1    Ada    Bea  0.582 0.272 0.847
#> 2    Bea    Ada  0.420 0.152 0.733
```

The effort for each direction should be how much that behaviour could
have been observed in that direction, for example in hours or sampling
periods. Edges are only estimated for the directions that appear in the
data, so rows with 0 events should be included for directions in which
the behaviour was never seen.

## Dyads that were never observed

Dyads with no rows in the data are not included in the model. If both
individuals were present but never seen together, include a row with 0
events and the corresponding effort. Dyads that could never have been
observed together should be left out.

## Group-by-individual data

Gambit-of-the-group data (where individuals seen in the same group are
assumed to be associating) are often stored as a group-by-individual
(GBI) matrix, where each row is a group or sampling event and each
column an individual, with 1 if the individual was present:

``` r

gbi <- matrix(rbinom(30 * 6, 1, 0.4), nrow = 30, ncol = 6,
              dimnames = list(NULL, c("Ada", "Bea", "Cat", "Dot", "Eve", "Fay")))
head(gbi)
#>      Ada Bea Cat Dot Eve Fay
#> [1,]   1   0   0   0   1   1
#> [2,]   1   1   0   1   0   1
#> [3,]   0   0   1   1   0   1
#> [4,]   1   1   0   0   1   1
#> [5,]   0   0   0   0   1   1
#> [6,]   1   0   0   0   0   1
```

[`convert_gbi_to_bison()`](https://bisonr.org/reference/convert_gbi_to_bison.md)
converts this to the bisonR format. For each pair of individuals and
each group in which at least one of the pair was present, it adds a row
with `event = 1` if both were present and `event = 0` if only one was.
Groups in which neither was present are not included. Each row is one
sampling period, so the effort is 1:

``` r

df_gbi <- convert_gbi_to_bison(gbi)
head(df_gbi)
#>   node_1 node_2 event group_id
#> 1    Ada    Bea     0        1
#> 2    Ada    Bea     1        2
#> 3    Ada    Bea     1        4
#> 4    Ada    Bea     0        6
#> 5    Ada    Bea     0        7
#> 6    Ada    Bea     0        8

fit_gbi <- bison_model(
  (event | 1) ~ dyad(node_1, node_2),
  data = df_gbi,
  model_type = "binary_conjugate",
  priors = get_default_priors("binary_conjugate")
)
```

Group-level information (such as location) can be added with
`group_properties`, and individual traits with `individual_properties`.
If some pairs could not have been seen together, for example because
they were in the study area at different times, use
`individual_constraints` to give each individual a value. Only pairs
with the same value are then included.

Rows from the same group are not independent, so a random effect of
group can be added to a (non-conjugate) model with `+ (1 | group_id)`.

## Common errors

[`bison_model()`](https://bisonr.org/reference/bison_model.md) checks
the data before fitting. The most common errors are:

- **Columns not found**: the names in the formula must exactly match
  columns in the data.
- **Missing values**: remove rows with `NA` in any column used by the
  model, or impute (fill in) the missing values.
- **More events than sampling periods**: in binary models, the number of
  events cannot exceed the effort. Use a count model for counts of
  interactions.
- **Self-loops**: rows where both nodes are the same individual aren’t
  allowed.
