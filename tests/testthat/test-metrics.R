test_that("network metrics work", {
  sim_data <- simulate_bison_model("binary", aggregated = TRUE)
  df <- sim_data$df_sim

  fit_edge <- bison_model(
    (event | duration) ~ dyad(node_1_id, node_2_id),
    data=df,
    model_type="binary_conjugate",
    priors=get_default_priors("binary_conjugate")
  )

  x <- expect_warning(
    extract_metric(fit_edge, "edge_weight", num_draws=10),
    regexp=NA
  )
  expect_equal(dim(x)[1], 10)

  x <- expect_warning(
    extract_metric(fit_edge, "node_strength", num_draws=50),
    regexp=NA
  )
  expect_equal(dim(x)[1], 50)

  net <- bison_to_igraph(fit_edge, 2)[[1]]
  expect_true(cor(igraph::strength(net), x[1, ]) > 0)

  x <- expect_error(
    extract_metric(fit_edge, "not_a_metric", standardise=TRUE)
  )

  x <- expect_warning(
    extract_metric(fit_edge, "node_strength", standardise=TRUE),
    regexp=NA
  )

  x <- expect_warning(
    extract_metric(fit_edge, "node_betweenness"),
    regexp=NA
  )

  x <- expect_warning(
    extract_metric(fit_edge, "node_eigen"),
    regexp=NA
  )

  x <- expect_warning(
    extract_metric(fit_edge, "node_degree[0.2]"),
    regexp=NA
  )

  x <- expect_warning(
    extract_metric(fit_edge, "node_clustering[0.2]"),
    regexp=NA
  )

  x <- expect_warning(
    extract_metric(fit_edge, "global_cv"),
    regexp=NA
  )

  x <- expect_warning(
    extract_metric(fit_edge, "global_density"),
    regexp=NA
  )

  x <- expect_warning(
    extract_metric(fit_edge, "global_std"),
    regexp=NA
  )

  x <- expect_warning(
    extract_metric(fit_edge, "global_diameter"),
    regexp=NA
  )

  x <- expect_warning(
    extract_metric(fit_edge, "global_clustering[0.2]"),
    regexp=NA
  )

  expect_warning (
    plot_metric(x),
    regexp=NA
  )


})
