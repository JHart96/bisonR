nodes <- c("A", "B", "C", "D")
dyads <- rbind(c("A", "B"), c("B", "C"), c("A", "C"), c("C", "D"))
weights <- c(0.5, 0.25, 0.1, 0.8)

test_that("metrics are correct on a network with known edge weights", {
  fit <- fake_bison_model(nodes, dyads, weights)
  strength <- extract_metric(fit, "node_strength", num_draws = 3)
  expect_equal(colnames(strength), nodes)
  expect_equal(unname(strength[1, ]), c(0.6, 0.75, 1.15, 0.8))

  degree <- extract_metric(fit, "node_degree[0.2]", num_draws = 2)
  expect_equal(unname(degree[1, ]), c(1, 2, 2, 1))

  expect_equal(extract_metric(fit, "global_density", num_draws = 2)[1, 1], mean(weights))
  expect_equal(extract_metric(fit, "global_cv", num_draws = 2)[1, 1], sd(weights) / mean(weights))
  expect_equal(extract_metric(fit, "global_std", num_draws = 2)[1, 1], sd(weights))

  # Triangle A-B-C exists only if the 0.1 edge counts
  expect_equal(extract_metric(fit, "global_clustering[0.05]", num_draws = 1)[1, 1],
               igraph::transitivity(igraph::make_graph(~ A - B - C - A, C - D), "global"))
  expect_equal(extract_metric(fit, "node_clustering[0.2]", num_draws = 1)[1, ], c(A = 0, B = 0, C = 0, D = 0))

  edge_weight <- extract_metric(fit, "edge_weight", num_draws = 2)
  expect_equal(unname(edge_weight[1, ]), weights)
  expect_equal(unname(extract_metric(fit, "edge_weight", num_draws = 2, transform = FALSE)[1, ]), qlogis(weights))
})

test_that("nodes with no edges are kept, in node order", {
  fit <- fake_bison_model(c(nodes, "E"), dyads, weights)
  strength <- extract_metric(fit, "node_strength", num_draws = 2)
  expect_equal(colnames(strength), c(nodes, "E"))
  expect_equal(unname(strength[1, 5]), 0)
  net <- bison_to_igraph(fit, num_draws = 1)[[1]]
  expect_equal(igraph::V(net)$name, c(nodes, "E"))
})

test_that("node metrics stay aligned with node names when rows are shuffled (issue #63)", {
  set.seed(10)
  sim <- simulate_bison_model("binary", aggregated = TRUE)
  df <- sim$df_sim[sample(nrow(sim$df_sim)), ]
  levels(df$node_1_id) <- levels(df$node_2_id) <- LETTERS[1:10]
  fit <- bison_model((event | duration) ~ dyad(node_1_id, node_2_id), df, "binary_conjugate",
                     priors = conjugate_priors("binary"))
  strength <- extract_metric(fit, "node_strength", num_draws = 200)
  edgelist <- get_edgelist(fit)
  by_hand <- sapply(LETTERS[1:10], function(n) {
    sum(edgelist$median[edgelist$node_1 == n | edgelist$node_2 == n])
  })
  expect_equal(colnames(strength), LETTERS[1:10])
  expect_equal(unname(apply(strength, 2, median)), unname(by_hand), tolerance = 0.15)
})

test_that("standardise centres metrics within each draw", {
  fit <- fake_bison_model(nodes, dyads, weights)
  x <- extract_metric(fit, "node_strength", num_draws = 4, standardise = TRUE)
  expect_equal(unname(rowMeans(x)), rep(0, 4))
})

test_that("all metric names work and unknown metrics fail", {
  fit <- fake_bison_model(nodes, dyads, weights)
  for (metric in c("node_eigen", "node_betweenness", "node_closeness", "global_diameter")) {
    expect_true(all(is.finite(extract_metric(fit, metric, num_draws = 2))), info = metric)
  }
  expect_error(extract_metric(fit, "not_a_metric"), "does not exist")
})

test_that("plot_metric handles continuous and discrete metrics", {
  set.seed(11)
  sim <- simulate_bison_model("binary", aggregated = TRUE)
  fit <- bison_model((event | duration) ~ dyad(node_1_id, node_2_id), sim$df_sim, "binary_conjugate",
                     priors = conjugate_priors("binary"))
  pdf(NULL)
  on.exit(dev.off())
  expect_silent(plot_metric(extract_metric(fit, "global_cv", num_draws = 200)))
  expect_silent(plot_metric(extract_metric(fit, "node_degree[0.5]", num_draws = 200)[, 1:2]))
})
