test_that("binary conjugate works", {
  sim_data <- simulate_bison_model("binary", aggregated = TRUE)
  df <- sim_data$df_sim

  priors = expect_warning(
    get_default_priors("binary_conjugate"),
    regexp=NA
  )

  expect_warning(
    prior_check(priors, "binary_conjugate"),
    regexp=NA
  )

  expect_warning(
    prior_predictive_check(
      (event | duration) ~ dyad(node_1_id, node_2_id),
      data=df,
      model_type="binary_conjugate",
      priors=priors
    ),
    regexp=NA
  )

  expect_warning(
    fit_edge <- bison_model(
      (event | duration) ~ dyad(node_1_id, node_2_id),
      data=df,
      data_type="binary_conjugate",
      priors=priors
    ),
    regexp=NA
  )

  expect_warning(
    plot_predictions(fit_edge, num_draws=20, type="density"),
    regexp=NA
  )

})


test_that("count conjugate works", {
  sim_data <- simulate_bison_model("count", aggregated = TRUE)
  df <- sim_data$df_sim

  priors = expect_warning(
    get_default_priors("count_conjugate"),
    regexp=NA
  )

  expect_warning(
    prior_check(priors, "count_conjugate"),
    regexp=NA
  )

  expect_warning(
    prior_predictive_check(
      (event | duration) ~ dyad(node_1_id, node_2_id),
      data=df,
      model_type="count_conjugate",
      priors=priors
    ),
    regexp=NA
  )

  expect_warning(
    fit_edge <- bison_model(
      (event | duration) ~ dyad(node_1_id, node_2_id),
      data=df,
      data_type="count_conjugate",
      priors=priors
    ),
    regexp=NA
  )

  expect_warning(
    plot_predictions(fit_edge, num_draws=20, type="density"),
    regexp=NA
  )

})
