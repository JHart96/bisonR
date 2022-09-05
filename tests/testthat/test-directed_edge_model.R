test_that("directed edge models work", {
  library(igraph)

  sim_data <- simulate_bison_model("binary", aggregated = FALSE)
  df <- sim_data$df_sim
  head(df)

  fit_edge <- expect_warning(
    bison_model (
      (event | duration) ~ dyad(node_1_id, node_2_id),
      data=df,
      data_type="binary",
      directed=T,
      priors=get_default_priors("binary")
    ),
    regexp=NA
  )

  expect_warning(
    get_edgelist(fit_edge),
    regexp=NA
  )

  expect_warning(
    summary(fit_edge),
    regexp=NA
  )

  expect_output(
    print(summary(fit_edge))
  )

  expect_warning(
    plot_network(fit_edge, lwd=1),
    regexp=NA
  )

  fit_mixture <- expect_warning(
    edge_mixture(fit_edge, num_components=3, verbose=TRUE),
    regexp=NA
  )

  expect_warning(
    extract_metric(fit_edge, "global_cv"),
    regexp=NA
  )

})
