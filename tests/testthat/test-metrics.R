test_that("network metrics work", {
  sim_data <- simulate_edge_model("binary", aggregated = TRUE)
  df <- sim_data$df_sim

  fit_edge <- edge_model(
    (event | duration) ~ dyad(node_1_id, node_2_id),
    data=df,
    data_type="binary_conjugate",
    priors=get_default_priors("binary_conjugate")
  )

  x <- expect_warning(
    draw_node_metric_samples(fit_edge, "strength"),
    regexp=NA
  )

  fit_edge <- edge_model(
    (event | duration) ~ dyad(node_1_id, node_2_id),
    data=df,
    data_type="count_conjugate",
    priors=get_default_priors("binary_conjugate")
  )

  x <- expect_error(
    draw_node_metric_samples(fit_edge, "not_a_metric", standardise=TRUE),
    regexp="Network metric does not exist"
  )

  x <- expect_warning(
    draw_node_metric_samples(fit_edge, "strength", standardise=TRUE),
    regexp=NA
  )

  x <- expect_warning(
    draw_node_metric_samples(fit_edge, "betweenness"),
    regexp=NA
  )

  x <- expect_warning(
    draw_node_metric_samples(fit_edge, "eigenvector"),
    regexp=NA
  )

  x <- expect_warning(
    draw_node_metric_samples(fit_edge, "degree[0.5]"),
    regexp=NA
  )

  x <- expect_warning(
    draw_network_metric_samples(fit_edge, "social_differentiation"),
    regexp=NA
  )

  x <- expect_warning(
    draw_network_metric_samples(fit_edge, "weighted_density"),
    regexp=NA
  )

  x <- expect_warning(
    draw_network_metric_samples(fit_edge, "standard_deviation"),
    regexp=NA
  )

  expect_output(
    print(x)
  )

  expect_warning(
    plot(x),
    regexp=NA
  )
})