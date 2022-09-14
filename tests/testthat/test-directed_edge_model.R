test_that("directed edge models work", {
  library(igraph)

  sim_data <- simulate_bison_model("binary", aggregated = FALSE)
  df <- sim_data$df_sim
  head(df)

  fit_edge <- expect_warning(
    bison_model (
      (event | duration) ~ dyad(node_1_id, node_2_id),
      data=df,
      model_type="binary",
      directed=TRUE,
      priors=get_default_priors("binary")
    ),
    regexp=NA
  )

  get_dyad_ids(df$node_1_id, df$node_2_id, fit_edge$dyad_to_idx, directed=TRUE)
  get_dyad_ids(df$node_2_id, df$node_1_id, fit_edge$dyad_to_idx, directed=TRUE)


  expect_warning(
    get_edgelist(fit_edge),
    regexp=NA
  )

  expect_warning(
    plot_predictions(fit_edge),
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
    bison_mixture(fit_edge, num_components=3, verbose=TRUE),
    regexp=NA
  )

  expect_warning(
    extract_metric(fit_edge, "global_cv"),
    regexp=NA
  )

})
