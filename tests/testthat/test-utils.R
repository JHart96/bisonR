test_that("convert_gbi_to_bison records co-occurrences correctly", {
  set.seed(1)
  gbi <- matrix(rbinom(200 * 10, 1, 0.25), 200, 10)
  df <- convert_gbi_to_bison(gbi)
  expect_equal(df$event, gbi[cbind(df$group_id, df$node_1)] * gbi[cbind(df$group_id, df$node_2)])
  # Every row has at least one of the pair present
  expect_true(all(gbi[cbind(df$group_id, df$node_1)] + gbi[cbind(df$group_id, df$node_2)] >= 1))
  # Every informative group is included for each pair
  expect_equal(sum(df$node_1 == 1 & df$node_2 == 2), sum(gbi[, 1] | gbi[, 2]))
})

test_that("convert_gbi_to_bison uses column names and respects constraints", {
  gbi <- matrix(c(1, 1, 0, 1, 0, 1, 1, 1, 1), 3, 3, dimnames = list(NULL, c("x", "y", "z")))
  df <- convert_gbi_to_bison(gbi, group_properties = c("a", "b", "c"),
                             individual_properties = c(1, 2, 3),
                             individual_constraints = c(1, 1, 2))
  expect_equal(unique(paste(df$node_1, df$node_2)), "x y")
  expect_equal(colnames(df), c("node_1", "node_2", "event", "group_id", "group_property",
                               "node_1_property", "node_2_property"))
})

test_that("convert_duration_to_binary rounds to whole sampling periods", {
  df <- data.frame(event = c(1.1, 2.9, 0), duration = c(3.2, 4.8, 1))
  converted <- convert_duration_to_binary(df, "event", "duration", 1)
  expect_equal(converted$event, c(1, 3, 0))
  expect_equal(converted$duration, c(3, 5, 1))
})

test_that("bison_to_igraph returns networks aligned with the edge list", {
  sim <- simulate_bison_model("binary", aggregated = TRUE)
  df <- sim$df_sim
  levels(df$node_1_id) <- levels(df$node_2_id) <- LETTERS[1:10]
  fit <- bison_model((event | duration) ~ dyad(node_1_id, node_2_id), df, "binary_conjugate",
                     priors = conjugate_priors("binary"))
  nets <- bison_to_igraph(fit, num_draws = 10)
  expect_length(nets, 10)
  expect_s3_class(nets[[1]], "igraph")
  expect_equal(igraph::as_edgelist(nets[[1]])[, 1], get_edgelist(fit)$node_1)
  expect_equal(igraph::as_edgelist(nets[[1]])[, 2], get_edgelist(fit)$node_2)
  expect_true(all(igraph::E(nets[[1]])$weight > 0 & igraph::E(nets[[1]])$weight < 1))
  # A single draw works too (issue #73)
  expect_length(bison_to_igraph(fit, num_draws = 1), 1)
})

test_that("draw_edgelist_samples returns the requested draws", {
  sim <- simulate_bison_model("binary", aggregated = TRUE)
  fit <- bison_model((event | duration) ~ dyad(node_1_id, node_2_id), sim$df_sim, "binary_conjugate",
                     priors = conjugate_priors("binary"))
  draws <- draw_edgelist_samples(fit, num_draws = 3)
  expect_equal(colnames(draws), c("node_1", "node_2", "draw.1", "draw.2", "draw.3"))
  expect_true(all(draws$draw.1 < 0 | draws$draw.1 >= 0))
  natural <- draw_edgelist_samples(fit, num_draws = 2, transform = TRUE)
  expect_true(all(natural$draw.1 > 0 & natural$draw.1 < 1))
})

test_that("plot_network draws every node and supports thresholds", {
  sim <- simulate_bison_model("count", aggregated = TRUE)
  fit <- bison_model((event | duration) ~ dyad(node_1_id, node_2_id), sim$df_sim, "count_conjugate",
                     priors = conjugate_priors("count"))
  pdf(NULL)
  on.exit(dev.off())
  layout <- plot_network(fit, threshold = 1e6)
  expect_equal(nrow(layout), 10)
  expect_silent(plot_network(fit, threshold = 1, threshold_on = "lower", layout = layout))
})

test_that("save and load round-trips a conjugate model", {
  sim <- simulate_bison_model("binary", aggregated = TRUE)
  fit <- bison_model((event | duration) ~ dyad(node_1_id, node_2_id), sim$df_sim, "binary_conjugate",
                     priors = conjugate_priors("binary"))
  path <- file.path(withr::local_tempdir(), "model")
  saved <- bison_save(fit, path)
  expect_true(file.exists(saved))
  expect_equal(bison_load(path)$edge_samples, fit$edge_samples)
})
