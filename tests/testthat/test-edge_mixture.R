test_that("binary edge mixtures work", {
  sim_data <- simulate_bison_model_mixture("binary", num_components = 2, component_weights = c(0.5, 0.5))
  df <- sim_data$df_sim

  fit_edge <- bison_model(
    (event | duration) ~ dyad(node_1_id, node_2_id),
    data=df,
    model_type="binary_conjugate",
    priors=get_default_priors("binary_conjugate")
  )

  fit_mixture <- expect_warning(
    bison_mixture(fit_edge, num_components=3, verbose=TRUE),
    regexp=NA
  )

  fit_mixture_summary <- expect_warning(summary(fit_mixture), regexp=NA)
  expect_output(print(fit_mixture))
  expect_output(print(fit_mixture_summary))

  expect_warning(
    get_edge_component_probabilities(fit_mixture, 2),
    regexp=NA
  )

  expect_warning(
    get_network_component_probabilities(fit_mixture),
    regexp=NA
  )

  expect_warning(
    get_component_means(fit_mixture, 3),
    regexp=NA
  )
})

test_that("count edge mixtures work", {
  sim_data <- simulate_bison_model_mixture("count", num_components = 2, component_weights = c(0.5, 0.5))
  df <- sim_data$df_sim

  fit_edge <- bison_model(
    (event | duration) ~ dyad(node_1_id, node_2_id),
    data=df,
    model_type="count_conjugate",
    priors=get_default_priors("count_conjugate")
  )

  fit_mixture <- expect_warning(
    bison_mixture(fit_edge, num_components=3, verbose=TRUE),
    regexp=NA
  )

  fit_mixture_summary <- expect_warning(summary(fit_mixture), regexp=NA)
  expect_output(print(fit_mixture))
  expect_output(print(fit_mixture_summary))

  expect_warning(
    get_edge_component_probabilities(fit_mixture, 2),
    regexp=NA
  )

  expect_warning(
    get_network_component_probabilities(fit_mixture),
    regexp=NA
  )

  expect_warning(
    get_component_means(fit_mixture, 3),
    regexp=NA
  )
})
