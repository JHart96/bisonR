# Convert a group-by-individual matrix to bisonR format

Group-by-individual (GBI) data record which individuals were seen
together in each group or sampling event, as a matrix with one row per
group and one column per individual. Use this function to prepare such
data for [`bison_model()`](https://bisonr.org/reference/bison_model.md).
It produces one row for each group and each pair of individuals where
*at least one* of the pair was present. The `event` column is 1 if both
were present and 0 if only one was. Groups where neither individual was
seen carry no information about the pair and are left out.

## Usage

``` r
convert_gbi_to_bison(
  gbi,
  group_properties = NULL,
  individual_properties = NULL,
  individual_constraints = NULL
)
```

## Arguments

- gbi:

  A matrix with R rows (groups) and C columns (individuals) of 0s and
  1s. If it has column names they are used as node names; otherwise
  nodes are numbered `1:C`.

- group_properties:

  Optional vector of length R with a property of each group (e.g.
  location), added as column `group_property`.

- individual_properties:

  Optional vector of length C with a property of each individual, added
  as columns `node_1_property` and `node_2_property`.

- individual_constraints:

  Optional vector of length C. Only pairs of individuals with the same
  value are included (for example, individuals that could have been seen
  together because they were in the study area at the same time).

## Value

A data frame with columns `node_1`, `node_2`, `event` and `group_id`,
plus any property columns, ready for use with a binary
[`bison_model()`](https://bisonr.org/reference/bison_model.md) with an
effort of 1, e.g. `(event | 1) ~ dyad(node_1, node_2)`.

## See also

[`bison_model()`](https://bisonr.org/reference/bison_model.md)

## Examples

``` r
gbi <- matrix(rbinom(20 * 5, 1, 0.4), 20, 5, dimnames = list(NULL, LETTERS[1:5]))
df <- convert_gbi_to_bison(gbi)
head(df)
#>   node_1 node_2 event group_id
#> 1      A      B     0        4
#> 2      A      B     1        5
#> 3      A      B     0        6
#> 4      A      B     0        8
#> 5      A      B     0        9
#> 6      A      B     0       11
```
