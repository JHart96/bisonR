test_that("convert_gbi_to_bison works", {
  gbi <- matrix(rbinom(20 * 10, 1, 0.25), 20, 10)
  df <- convert_gbi_to_bison(gbi)

  # Correct length
  expect_true(nrow(df) == 20 * 10 * 9 * 0.5)

  # Correct number of nodes
  expect_true(length(unique(c(df$node_1, df$node_2))) == ncol(gbi))

  # Entries are correct
  results <- rep(0, nrow(df))
  for (i in 1:nrow(df)) {
    results[i] <- gbi[ df[i, ]$group_id, df[i, ]$node_1] * gbi[ df[i, ]$group_id, df[i, ]$node_2] == df[i, ]$event
  }
  expect_true(all(results == TRUE))
})

test_that("convert_gbi_to_bison with constraints/properties works", {
  gbi <- matrix(rbinom(20 * 10, 1, 0.25), 20, 10)
  group_properties <- rnorm(20)
  individual_properties <- rnorm(10)
  individual_constraints <- sample(1:2, 10, replace=TRUE)
  df <- expect_warning(
    convert_gbi_to_bison(gbi, group_properties, individual_properties, individual_constraints),
    regexp=NA
  )
})

test_that("convert duration to binary work", {
  sim_data <- simulate_bison_model("binary", aggregated = TRUE)
  df <- sim_data$df_sim
  df_true <- df
  df$event <- abs(rnorm(nrow(df), df$event, 0.1))
  df$duration <- abs(rnorm(nrow(df), df$duration, 0.1))
  df_converted <- convert_duration_to_binary(df, "event", "duration", 1)
  expect_true(all(df_true$event == df_converted$event))
  expect_true(all(df_true$duration == df_converted$duration))
})

test_that("bison to igraph conversion works", {
  sim_data <- bisonR::simulate_bison_model("binary", aggregated = TRUE)
  df_sim <- sim_data$df_sim
  levels(df_sim$node_1_id) <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")
  levels(df_sim$node_2_id) <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")

  fit_edge <- bison_model(
    (event | duration) ~ dyad(node_1_id, node_2_id),
    data=df_sim,
    model_type="binary_conjugate",
    priors=get_default_priors("binary_conjugate")
  )

  nets <- expect_warning(
    bison_to_igraph(fit_edge, num_draws=10),
    regexp=NA
  )

  expect_true(length(nets) == 10)
  expect_true(class(nets[[1]]) == "igraph")

  expect_true(all(igraph::as_edgelist(nets[[1]])[, 1] == get_edgelist(fit_edge)$node_1))
  expect_true(all(igraph::as_edgelist(nets[[1]])[, 2] == get_edgelist(fit_edge)$node_2))

})
