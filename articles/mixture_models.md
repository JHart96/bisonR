# Edge mixture models

## Introduction

Edge mixture models can be used to detect different types of social
connections by looking for clustering in edge weights. bisonR provides
support for edge mixture models with the `bison_mixture` function. There
are two key outputs from a bisonR edge mixture model. First it generates
a posterior probability distribution over the number of components
(clusters, interpreted associal connection types), from 1 to a maximum
$`K`$, over the entire network. Secondly, for each edge it calculates
the probability that each dyad belongs to each possible component. To
show how this all works, first we’ll load in bisonR and simulate some
data:

``` r

library(bisonR)
#> Loading required package: cmdstanr
#> This is cmdstanr version 0.9.0.9002
#> - CmdStanR documentation and vignettes: mc-stan.org/cmdstanr
#> - CmdStan path: /home/runner/.cmdstan/cmdstan-2.40.0
#> - CmdStan version: 2.40.0
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
library(igraph)
#> 
#> Attaching package: 'igraph'
#> The following objects are masked from 'package:dplyr':
#> 
#>     as_data_frame, groups, union
#> The following objects are masked from 'package:stats':
#> 
#>     decompose, spectrum
#> The following object is masked from 'package:base':
#> 
#>     union
library(brms)
#> Loading required package: Rcpp
#> Loading 'brms' package (version 2.23.0). Useful instructions
#> can be found by typing help('brms'). A more detailed introduction
#> to the package is available through vignette('brms_overview').
#> 
#> Attaching package: 'brms'
#> The following object is masked from 'package:stats':
#> 
#>     ar
library(mclust)
#> Package 'mclust' version 6.1.4
#> Type 'citation("mclust")' for citing this R package in publications.
#> 
#> Attaching package: 'mclust'
#> The following object is masked from 'package:brms':
#> 
#>     me
#> The following object is masked from 'package:dplyr':
#> 
#>     count
```

In this example we simulate data from a binary edge model with 2
underlying clusters (also called components), with dyads belong to
either cluster with equal probability:

``` r

sim_data <- simulate_bison_model_mixture("binary", num_components = 2, component_weights = c(0.5, 0.5))
df <- sim_data$df_sim
head(df)
#> # A tibble: 6 × 5
#> # Groups:   node_1_id [1]
#>   node_1_id node_2_id event duration age_diff
#>   <fct>     <fct>     <dbl>    <dbl>    <dbl>
#> 1 1         2            10       11        1
#> 2 1         3             0        1        1
#> 3 1         4            13       13        1
#> 4 1         5            10       10        1
#> 5 1         6             3       12        1
#> 6 1         7            14       16        1
```

Since the data are already aggregated (e.g. a single row for each dyad),
we can use the binary conjugate model for speed:

``` r

fit_edge <- bison_model(
  (event | duration) ~ dyad(node_1_id, node_2_id), 
  data=df, 
  model_type="binary_conjugate"
)
#> No priors set by user, using default priors instead. We recommend setting and checking priors explicitly for reliable inference.
```

We won’t conduct the usual checks for the sake of this tutorial, but at
this point we’d recommend conducting diagnostic checks to ensure model
fit. Once the model has fit well, we can use the
[`bison_mixture()`](https://bisonr.org/reference/bison_mixture.md)
function to fit the edge mixture model:

``` r

# verbose=FALSE for tutorial purposes
fit_mixture <- bison_mixture(fit_edge, num_components=5, verbose=FALSE)
summary(fit_mixture)
#> === Fitted dyadic mixture model ===
#> Maximum number of components: 
#> Best model: 1 components
#> Probability of best model: 73.5%
#> === Component probabilities ===
#>            1     2     3     4 5
#> P(K=k) 0.735 0.231 0.019 0.015 0
#> === Component means for best model (K = 1) ===
#>              1
#> mean 0.5471098
#> === Edge component probabilities for best model (K = 1) ===
#>         1
#> 1 <-> 2 1
#> 1 <-> 3 1
#> 1 <-> 4 1
#> 1 <-> 5 1
#> 1 <-> 6 1
#> 1 <-> 7 1
#> ...
```

There’s quite a bit going on in the summary, so let’s break it down. The
top section shows the model that fits best and the probability assigned
to that model. The second section shows the probabilities of the models
corresponding to each number of components ($`K = 1, 2, 3, ...`$).
You’ll notice the model with the highest probability is the same as the
best fitting model from the section above. The final section shows the
edge-level probabilities of component membership for the corresponding
dyad. This can be interpreted as the probability a dyad belongs to a
given component/cluster/social connection type.

The information shown in the summary above can also be accessed using
the functions
[`get_network_component_probabilities()`](https://bisonr.org/reference/get_network_component_probabilities.md)
and
[`get_edge_component_probabilities()`](https://bisonr.org/reference/get_edge_component_probabilities.md).
These two functions make it possible to access the mixture model output
programmatically for downstream analysis or plotting.

To get network-level probabilities of the number of components, we can
use
[`get_network_component_probabilities()`](https://bisonr.org/reference/get_network_component_probabilities.md)
function:

``` r

get_network_component_probabilities(fit_mixture)
#>   num_components  probability
#> 1              1 0.7345254805
#> 2              2 0.2312023503
#> 3              3 0.0188747949
#> 4              4 0.0152140086
#> 5              5 0.0001833657
```

To use probabilities over component membership for a given number of
components, we can use the
[`get_edge_component_probabilities()`](https://bisonr.org/reference/get_edge_component_probabilities.md)
function, where `3` is the number of components we assume to exist in
the network:

``` r

get_edge_component_probabilities(fit_mixture, 3)
#>     node_1 node_2 P(K = 1) P(K = 2) P(K = 3)
#> 1        1      2    0.015    0.425    0.560
#> 2        1      3    0.605    0.260    0.135
#> 3        1      4    0.000    0.310    0.690
#> 4        1      5    0.000    0.340    0.660
#> 5        1      6    0.790    0.155    0.055
#> 6        1      7    0.030    0.460    0.510
#> 7        1      8    0.935    0.060    0.005
#> 8        1      9    0.925    0.075    0.000
#> 9        1     10    0.925    0.070    0.005
#> 10       1     11    0.850    0.125    0.025
#> 11       1     12    0.895    0.095    0.010
#> 12       1     13    0.960    0.040    0.000
#> 13       1     14    0.650    0.250    0.100
#> 14       1     15    0.005    0.270    0.725
#> 15       2      3    0.860    0.115    0.025
#> 16       2      4    0.015    0.455    0.530
#> 17       2      5    0.085    0.430    0.485
#> 18       2      6    0.035    0.355    0.610
#> 19       2      7    0.000    0.385    0.615
#> 20       2      8    0.935    0.055    0.010
#> 21       2      9    0.795    0.190    0.015
#> 22       2     10    0.930    0.070    0.000
#> 23       2     11    0.925    0.075    0.000
#> 24       2     12    0.010    0.450    0.540
#> 25       2     13    0.920    0.075    0.005
#> 26       2     14    0.090    0.410    0.500
#> 27       2     15    0.935    0.055    0.010
#> 28       3      4    0.170    0.450    0.380
#> 29       3      5    0.915    0.080    0.005
#> 30       3      6    0.005    0.490    0.505
#> 31       3      7    0.000    0.420    0.580
#> 32       3      8    0.000    0.375    0.625
#> 33       3      9    0.230    0.430    0.340
#> 34       3     10    0.000    0.395    0.605
#> 35       3     11    0.835    0.135    0.030
#> 36       3     12    0.930    0.060    0.010
#> 37       3     13    0.940    0.055    0.005
#> 38       3     14    0.910    0.085    0.005
#> 39       3     15    0.940    0.060    0.000
#> 40       4      5    0.625    0.220    0.155
#> 41       4      6    0.035    0.405    0.560
#> 42       4      7    0.015    0.420    0.565
#> 43       4      8    0.935    0.060    0.005
#> 44       4      9    0.000    0.310    0.690
#> 45       4     10    0.005    0.405    0.590
#> 46       4     11    0.935    0.065    0.000
#> 47       4     12    0.245    0.455    0.300
#> 48       4     13    0.030    0.425    0.545
#> 49       4     14    0.940    0.055    0.005
#> 50       4     15    0.945    0.050    0.005
#> 51       5      6    0.085    0.415    0.500
#> 52       5      7    0.910    0.080    0.010
#> 53       5      8    0.660    0.255    0.085
#> 54       5      9    0.000    0.300    0.700
#> 55       5     10    0.925    0.070    0.005
#> 56       5     11    0.020    0.450    0.530
#> 57       5     12    0.005    0.380    0.615
#> 58       5     13    0.000    0.370    0.630
#> 59       5     14    0.005    0.405    0.590
#> 60       5     15    0.015    0.510    0.475
#> 61       6      7    0.000    0.350    0.650
#> 62       6      8    0.010    0.370    0.620
#> 63       6      9    0.670    0.245    0.085
#> 64       6     10    0.065    0.485    0.450
#> 65       6     11    0.090    0.440    0.470
#> 66       6     12    0.620    0.295    0.085
#> 67       6     13    0.000    0.300    0.700
#> 68       6     14    0.080    0.410    0.510
#> 69       6     15    0.945    0.045    0.010
#> 70       7      8    0.940    0.055    0.005
#> 71       7      9    0.000    0.435    0.565
#> 72       7     10    0.940    0.050    0.010
#> 73       7     11    0.025    0.475    0.500
#> 74       7     12    0.000    0.380    0.620
#> 75       7     13    0.000    0.360    0.640
#> 76       7     14    0.005    0.400    0.595
#> 77       7     15    0.890    0.090    0.020
#> 78       8      9    0.000    0.280    0.720
#> 79       8     10    0.945    0.055    0.000
#> 80       8     11    0.000    0.310    0.690
#> 81       8     12    0.925    0.065    0.010
#> 82       8     13    0.015    0.460    0.525
#> 83       8     14    0.925    0.075    0.000
#> 84       8     15    0.005    0.445    0.550
#> 85       9     10    0.000    0.405    0.595
#> 86       9     11    0.915    0.065    0.020
#> 87       9     12    0.000    0.320    0.680
#> 88       9     13    0.810    0.120    0.070
#> 89       9     14    0.000    0.435    0.565
#> 90       9     15    0.005    0.380    0.615
#> 91      10     11    0.005    0.455    0.540
#> 92      10     12    0.935    0.065    0.000
#> 93      10     13    0.020    0.380    0.600
#> 94      10     14    0.000    0.255    0.745
#> 95      10     15    0.935    0.065    0.000
#> 96      11     12    0.000    0.425    0.575
#> 97      11     13    0.000    0.415    0.585
#> 98      11     14    0.665    0.220    0.115
#> 99      11     15    0.960    0.030    0.010
#> 100     12     13    0.960    0.035    0.005
#> 101     12     14    0.025    0.475    0.500
#> 102     12     15    0.000    0.425    0.575
#> 103     13     14    0.910    0.090    0.000
#> 104     13     15    0.000    0.285    0.715
#> 105     14     15    0.225    0.460    0.315
```

It’s often useful to know where the mixture components are location, so
we can get their posterior means using the
[`get_component_means()`](https://bisonr.org/reference/get_component_means.md)
function for a given number of components. Here we’ll calculate the
posterior means of the components for the mixture model with 3
components:

``` r

get_component_means(fit_mixture, 3)
#>              50%         5%       95%
#> K = 1 0.09508668 0.03845895 0.1298062
#> K = 2 0.84673683 0.22354973 0.9533372
#> K = 3 0.90170414 0.78050683 0.9976114
```

## Using components in downstream analysis

Now say that we’re interested in the number of strong vs weak partners.
Specifically we want to know if the strength of connections only to
strong partners predicts a node-level trait. This is a non-standard
analysis so will require some manual data wrangling to fit the right
kind of model.

For demonstration purposes, we’ll create a dataframe that might be used
for this kind of analysis:

``` r

df_nodal <- data.frame(node_id=as.factor(levels(df$node_1_id)), trait=rnorm(15))
df_nodal
#>    node_id        trait
#> 1        1 -0.163956869
#> 2        2  0.417166376
#> 3        3 -1.026529143
#> 4        4  0.154483632
#> 5        5  0.049465437
#> 6        6 -0.418917019
#> 7        7  1.124576040
#> 8        8 -0.724778621
#> 9        9  0.106686091
#> 10      10  0.783754076
#> 11      11 -0.721932177
#> 12      12  0.465342170
#> 13      13  0.057968075
#> 14      14  0.811572887
#> 15      15 -0.003290682
```

This dataframe will provide the information about traits that we need
for the analysis. Next, we need to calculate strength using only
individuals categorised as strong according to the mixture model. Since
both component membership (strong or weak) *and* edge weights are
probabilistic, we have to sample from both posteriors as we build a new
posterior of mixture-based strength:

``` r

# Set number of nodes and number of posterior samples
num_nodes <- fit_edge$num_nodes
num_draws <- 5 # Keep short for demo purposes.

# Create a list of igraph networks from edgemodel to represent network posterior
nets <- bison_to_igraph(fit_edge, num_draws)

# Create an empty matrix to hold strengths of top mixture component
mix_strengths <- matrix(0, num_draws, num_nodes)

# Create an empty list for imputed versions of the dataframe
imputed_data <- list()

# Loop through each posterior sample and calculate strength of top mixture.
for (i in 1:num_draws) {
  # Calculate edge components (1 if strong, 0 if weak)
  edge_components <- 1 * (fit_mixture$edge_component_samples[i, 2, ] == 2)
  
  # If the edge is strong, don't change edge weight, but if it's weak then set to zero.
  E(nets[[i]])$weight <- E(nets[[i]])$weight * edge_components
  
  # Change the metric values of the imputed data to be mixture-based strength
  imputed_data[[i]] <- df_nodal
  imputed_data[[i]]$mix_strength <- strength(nets[[i]])
}
imputed_data[[2]]
#>    node_id        trait mix_strength
#> 1        1 -0.163956869     5.022821
#> 2        2  0.417166376     6.020381
#> 3        3 -1.026529143     5.139991
#> 4        4  0.154483632     6.716043
#> 5        5  0.049465437     7.994020
#> 6        6 -0.418917019     8.381191
#> 7        7  1.124576040     8.500934
#> 8        8 -0.724778621     5.236292
#> 9        9  0.106686091     8.102017
#> 10      10  0.783754076     6.080982
#> 11      11 -0.721932177     5.932365
#> 12      12  0.465342170     6.218257
#> 13      13  0.057968075     6.985885
#> 14      14  0.811572887     5.875518
#> 15      15 -0.003290682     5.162626
```

Now we have a set of imputed dataframes containing the posteriors of
mixture-based node strength, we can fit the model using brm_multiple.

``` r

fit_brm <- brm_multiple(trait ~ mix_strength, imputed_data, silent=2, refresh=0)
summary(fit_brm)
#>  Family: gaussian 
#>   Links: mu = identity 
#> Formula: trait ~ mix_strength 
#>    Data: imputed_data (Number of observations: 15) 
#>   Draws: 20 chains, each with iter = 2000; warmup = 1000; thin = 1;
#>          total post-warmup draws = 20000
#> 
#> Regression Coefficients:
#>              Estimate Est.Error l-95% CI u-95% CI
#> Intercept       -1.17      0.96    -3.07     0.75
#> mix_strength     0.19      0.14    -0.10     0.48
#> 
#> Further Distributional Parameters:
#>       Estimate Est.Error l-95% CI u-95% CI
#> sigma     0.64      0.15     0.42     0.99
#> 
#> Draws were sampled using sampling(NUTS). Overall Rhat and ESS estimates
#> are not informative for brm_multiple models and are hence not displayed.
#> Please see ?brm_multiple for how to assess convergence of such models.
```

This was a short tutorial on how to use
[`bison_mixture()`](https://bisonr.org/reference/bison_mixture.md) to
fit mixture models in bisonR. If you have any questions please get in
touch.
