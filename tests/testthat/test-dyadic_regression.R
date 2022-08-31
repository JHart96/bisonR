test_that("dyadic regression works", {
  sim_data <- simulate_edge_model("binary", aggregated = TRUE)
  df <- sim_data$df_sim
  head(df)

  fit_edge <- edge_model(
    (event | duration) ~ dyad(node_1_id, node_2_id),
    data=df,
    data_type="binary_conjugate",
    priors=get_default_priors("binary_conjugate")
  )

  coef <- 1
  df$trait <- rnorm(nrow(df), coef * fit_edge$chain[sample(1:100, size=1), ], 0.1)

  expect_warning(
    fit_dyad <- dyadic_regression(
      dyad(node_1_id, node_2_id) ~ trait,
      fit_edge,
      df,
      mm=FALSE
    ),
    regexp=NA
  )

  expect_warning(
    summary(fit_dyad),
    regexp=NA
  )

  expect_warning(
    plot_predictions(fit_dyad, num_draws=20, type="density"),
    regexp=NA
  )

  expect_warning(
    plot_predictions(fit_dyad, num_draws=20, type="marginal"),
    regexp=NA
  )

  coef <- 1
  df$trait <- rnorm(nrow(df), coef * fit_edge$chain[sample(1:1000, size=1), ], 0.1)

  expect_warning(
    fit_dyad <- dyadic_regression(
      trait ~ dyad(node_1_id, node_2_id),
      fit_edge,
      df,
      mm=FALSE
    ),
    regexp=NA
  )

  expect_warning(
    plot_predictions(fit_dyad, num_draws=20, type="density"),
    regexp=NA
  )

  expect_warning(
    plot_predictions(fit_dyad, num_draws=20, type="marginal"),
    regexp=NA
  )

  expect_warning(
    summary(fit_dyad),
    regexp=NA
  )

  ests <- fit_dyad$fit$summary()[4, 6:7]
  expect_equal(
    ests[1] < coef && ests[2] > coef,
    TRUE
  )
})
