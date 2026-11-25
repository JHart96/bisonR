# Network metrics

Network metrics can be extracted using the function
`extract_metric(bison_model, metric_name)`. The `metric_name` variable
is made up of two parts: the network feature and the network metric. For
example, `node_eigen` will calculate the eigenvector centrality of each
node.

The following feature/metric pairs are currently supported in bisonR:

- `edge` - Edge-level properties
  - `weight` - Edge weight
- `node` - Node-level properties
  - `degree[t]` - Node degree (where t is the binary threshold for an
    edge)
  - `strength` - Node strength
  - `eigen` - Node eigenvector centrality
  - `betweenness` - Node betweenness
  - `closeness` - Node closeness
  - `clustering` - Node-level clustering coefficient
- `global` - Network-level properties
  - `density` - Weighted density
  - `cv` - Coefficient of variation of edge weights - also known as
    social differentiation
  - `std` - Standard deviation of edge weights
  - `diameter` - Diameter of the network
  - `clustering` - Clustering coefficient of the network - also known as
    transitivity

## Examples

See below for examples of how to generate the network metrics listed
above.

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
```

``` r

sim_data <- simulate_bison_model("binary", aggregated = TRUE)
df <- sim_data$df_sim

fit_edge <- bison_model(
  (event | duration) ~ dyad(node_1_id, node_2_id),
  data=df,
  model_type="binary_conjugate",
  priors=get_default_priors("binary_conjugate")
)
```

### Edge

#### Edge weight

``` r

metrics <- extract_metric(fit_edge, "edge_weight", num_draws=10)
metrics[1:6, 1:5] # Preview first few samples and edges
#>            [,1]      [,2]      [,3]      [,4]      [,5]
#> [1,] -1.7763303 1.2071965 1.3689566 3.2731301 3.4882468
#> [2,]  0.9357479 0.9649296 1.0682155 1.1780742 0.8553972
#> [3,] -3.6258711 1.8834824 2.0402320 0.2606149 2.6947782
#> [4,] -2.7593545 3.0038231 1.1405372 1.8913551 5.6115489
#> [5,] -1.8913159 0.9642772 0.8170626 4.0066255 1.8581266
#> [6,]  0.1826523 1.2889817 1.4426996 0.8583163 2.2183483
```

### Node

#### Node degree

``` r

metrics <- extract_metric(fit_edge, "node_degree[0.2]")
metrics[1:6, 1:5] # Preview first few samples and nodes
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    8    7    5    6    8
#> [2,]    8    7    8    8    4
#> [3,]    8    5    7    8    7
#> [4,]    9    6    5    8    8
#> [5,]    7    6    8    7    6
#> [6,]    8    5    8    7    6
```

#### Node strength

``` r

metrics <- extract_metric(fit_edge, "node_strength", num_draws=10)
metrics[1:6, 1:5] # Preview first few samples and nodes
#>          [,1]     [,2]     [,3]     [,4]     [,5]
#> [1,] 5.250600 3.605330 4.044063 4.328342 4.147407
#> [2,] 6.390623 3.869768 4.278853 5.145356 4.451589
#> [3,] 6.153401 3.094590 4.164048 4.940242 4.975292
#> [4,] 6.359865 3.531231 5.004913 5.160148 4.461221
#> [5,] 5.756916 3.451609 4.778303 4.609269 4.311387
#> [6,] 5.174642 3.617823 5.412316 4.863402 3.786787
```

#### Node eigenvector centrality

``` r

metrics <- extract_metric(fit_edge, "node_eigen")
metrics[1:6, 1:5] # Preview first few samples and nodes
#>      [,1]      [,2]      [,3]      [,4]      [,5]
#> [1,]    1 0.6998075 0.8274551 0.7619760 0.7679569
#> [2,]    1 0.5417456 0.8373102 0.8503530 0.8265666
#> [3,]    1 0.6760519 0.8898821 0.9472859 0.9077447
#> [4,]    1 0.4732591 0.7720060 0.8812934 0.7572672
#> [5,]    1 0.4138459 0.8436132 0.9317520 0.8012959
#> [6,]    1 0.5517734 0.8733358 0.9190946 0.6936717
```

#### Node betweenness

``` r

metrics <- extract_metric(fit_edge, "node_betweenness")
metrics[1:6, 1:5] # Preview first few samples and nodes
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    4    0    0    4    1
#> [2,]    9    0    0    1    6
#> [3,]    7    0    0    6    0
#> [4,]   14    0    0    3    6
#> [5,]   12    0    3    1    2
#> [6,]    8    2    5    0    3
```

#### Node closeness

``` r

metrics <- extract_metric(fit_edge, "node_closeness")
metrics[1:6, 1:5] # Preview first few samples and nodes
#>            [,1]       [,2]       [,3]       [,4]       [,5]
#> [1,] 0.06313667 0.03415140 0.05568646 0.05583207 0.04863124
#> [2,] 0.05638048 0.03776727 0.05059571 0.05298703 0.05238034
#> [3,] 0.06271996 0.03884989 0.05684554 0.05761289 0.05130179
#> [4,] 0.05566669 0.02823852 0.05793152 0.04768573 0.04804053
#> [5,] 0.05791287 0.05182419 0.05418053 0.06348239 0.05845334
#> [6,] 0.06342691 0.03999048 0.05167693 0.05367846 0.04884259
```

#### Node-level clustering coefficient

``` r

metrics <- extract_metric(fit_edge, "node_clustering[0.2]")
metrics[1:6, 1:5] # Preview first few samples and nodes
#>           [,1]      [,2]      [,3]      [,4]      [,5]
#> [1,] 0.6111111 0.8333333 0.6666667 0.6666667 0.8000000
#> [2,] 0.7857143 0.7142857 0.7619048 0.9047619 0.9333333
#> [3,] 0.7142857 0.6071429 0.7333333 0.7333333 0.6666667
#> [4,] 0.7142857 0.6000000 0.6190476 0.6666667 0.8000000
#> [5,] 0.6428571 0.7000000 0.6666667 0.6785714 0.6666667
#> [6,] 0.5833333 0.8000000 0.7619048 0.5714286 0.9000000
```

### Global

#### Global coefficient of variation

``` r

metrics <- extract_metric(fit_edge, "global_cv")
metrics[1:6] # Preview first few samples
#> [1] 0.6620834 0.7752250 0.7524315 0.6514228 0.6969299 0.8203320
```

#### Global density

``` r

metrics <- extract_metric(fit_edge, "global_density")
metrics[1:6] # Preview first few samples
#> [1] 0.3930515 0.4329158 0.4216413 0.4753094 0.4301271 0.4341689
```

#### Global standard deviation of edge weights

``` r

metrics <- extract_metric(fit_edge, "global_std")
metrics[1:6] # Preview first few samples
#> [1] 0.3093545 0.3094215 0.3027798 0.3243845 0.3032741 0.3233578
```

#### Global diameter

``` r

metrics <- extract_metric(fit_edge, "global_diameter")
metrics[1:6] # Preview first few samples
#> [1] 0.3318503 0.2553142 0.4628383 0.2985160 0.4401158 0.5826086
```

#### Global clustering coefficient

``` r

metrics <- extract_metric(fit_edge, "global_clustering[0.2]")
metrics[1:6] # Preview first few samples
#> [1] 0.5827338 0.7288136 0.6580645 0.7853403 0.8000000 0.8111588
```
