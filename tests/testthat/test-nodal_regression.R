test_that("nodal regression works", {
  sim_data <- simulate_edge_model("binary", aggregated = TRUE)
  df <- sim_data$df_sim
  head(df)

  fit_edge <- edge_model(
    (event | duration) ~ dyad(node_1_id, node_2_id),
    data=df,
    data_type="binary_conjugate",
    priors=get_default_priors("binary_conjugate")
  )

  df_nodal <- data.frame(node_id=levels(df$node_1_id))
  node_metric <- apply(draw_node_metric_samples(fit_edge, "strength"), 2, mean)
  df_nodal$trait <- rnorm(nrow(df_nodal), 0, 0.1)

  expect_warning(
    fit_nodal <- nodal_regression(
      strength(node_id) ~ trait,
      fit_edge,
      df_nodal
    ),
    regexp=NA
  )

  expect_warning(
    summary(fit_nodal),
    regexp=NA
  )

  expect_warning(
    plot_predictions(fit_nodal, num_draws=20, type="density"),
    regexp=NA
  )

  expect_warning(
    plot_predictions(fit_nodal, num_draws=20, type="marginal"),
    regexp=NA
  )

  expect_warning(
    fit_nodal <- nodal_regression(
      trait ~ strength(node_id),
      fit_edge,
      df_nodal
    ),
    regexp=NA
  )

  expect_warning(
    summary(fit_nodal),
    regexp=NA
  )

  expect_warning(
    plot_predictions(fit_nodal, num_draws=20, type="density"),
    regexp=NA
  )

  expect_warning(
    plot_predictions(fit_nodal, num_draws=20, type="marginal"),
    regexp=NA
  )

})
