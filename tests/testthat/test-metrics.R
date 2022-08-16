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
    draw_network_metric_samples(fit_edge, "social_differentiation"),
    regexp=NA
  )
  plot(x)

  x <- expect_warning(
    draw_network_metric_samples(fit_edge, "weighted_density"),
    regexp=NA
  )
  plot(x)

  x <- expect_warning(
    draw_network_metric_samples(fit_edge, "standard_deviation"),
    regexp=NA
  )
  plot(x)
})
