# Convert group-by-individual matrix to bisonR format

Convert group-by-individual matrix to bisonR format

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

  Group-by-individual (R x C) matrix corresponding to R groups of C
  individuals.

- group_properties:

  A vector of length R describing properties of the groups to be added
  to the dataframe.

- individual_properties:

  A vector of length C describing properties of the individuals to be
  added to the dataframe.

- individual_constraints:

  A vector of length C describing constraints of the individuals.

## Value

A dataframe for use in bisonR, where each row corresponds to a possible
social event between two individuals.
