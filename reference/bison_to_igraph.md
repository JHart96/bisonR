# Convert bisonR edgemodel posterior to igraph network objects

Convert bisonR edgemodel posterior to igraph network objects

## Usage

``` r
bison_to_igraph(obj, num_draws = 1000)
```

## Arguments

- obj:

  bisonR edge model

- num_draws:

  Number of draws to sample from network posterior

## Value

List of igraph objects representing posterior distribution of networks
