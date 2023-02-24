test_that("Binary edge model parameter estimation", {
  library(dplyr)
  library(igraph)

  set.seed(123)

  # Load data in with minimal effects
  sim_data <- simulate_bison_model("binary", aggregated = TRUE, location_effect = FALSE, age_diff_effect = FALSE)
  df <- sim_data$df_sim
  df$group_id <- sample(1:4, nrow(df), replace=TRUE)

  priors = get_default_priors("binary")

  expect_warning(
    prior_check(priors, "binary"),
    regexp=NA
  )

  expect_warning(
    prior_predictive_check(
      (event | duration) ~ dyad(node_1_id, node_2_id),
      data=df,
      model_type="binary",
      priors=priors
    ),
    regexp=NA
  )

  # Fit model to the data
  fit_edge <- expect_warning(
    bison_model(
      (event | duration) ~ dyad(node_1_id, node_2_id) + age_diff + (1 | group_id),
      data=df,
      model_type="binary"
    ),
    regexp=NA
  )

  # Fit model to the data
  fit_pool <- expect_warning(
    bison_model(
      (event | duration) ~ dyad(node_1_id, node_2_id),
      data=df,
      model_type="binary",
      partial_pooling=TRUE,
      zero_inflated=TRUE
    ),
    regexp=NA
  )

  fit_null <- expect_warning (
    bison_model(
      (event | duration) ~ 1,
      data=df,
      model_type="binary",
      priors=get_default_priors("binary")
    ),
    regexp=NA
  )

  fit_compare <- suppressWarnings(model_comparison(list(non_random_model=fit_edge, random_model=fit_null)))

  expect_output(print(fit_compare))

  # Extract estimates and true values.
  true <- sim_data$df_true %>%
    mutate(edge_weight=edge_weight) %>%
    select(node_1=node_1_id, node_2=node_2_id, true=edge_weight)
  est <- get_edgelist(fit_edge, transform=FALSE) %>%
    select(node_1, node_2, est=median, est_lb="5%", est_ub="95%")
  comparison <- left_join(true, est, by=c("node_1", "node_2"))

  # Calculate correlation between estimates and true values.
  proportion_within_ci <- mean(comparison$true < comparison$est_ub & comparison$est > comparison$est_lb)
  expect_gt(proportion_within_ci, 0.9)

  # Check that plots don't produce warnings
  expect_warning(plot_predictions(fit_edge), regexp=NA)
  expect_warning(plot_network(fit_edge, lwd=10, threshold=0.1), regexp=NA)
  expect_warning(plot_trace(fit_edge, par_ids=1), regexp=NA)
  expect_output(print(summary(fit_edge)))
})

test_that("Count edge model parameter estimation", {
  library(dplyr)

  set.seed(123)

  # Load data in with minimal effects
  sim_data <- simulate_bison_model("count", aggregated = TRUE, location_effect = FALSE, age_diff_effect = FALSE)
  df <- sim_data$df_sim

  # Fit model to the data
  fit_edge <- bison_model(
    (event | duration) ~ dyad(node_1_id, node_2_id),
    data=df,
    model_type="count"
  )

  # Extract estimates and true values.
  true <- sim_data$df_true %>%
    mutate(edge_weight=edge_weight) %>%
    select(node_1=node_1_id, node_2=node_2_id, true=edge_weight)
  est <- get_edgelist(fit_edge) %>%
    select(node_1, node_2, est=median, est_lb="5%", est_ub="95%")
  comparison <- left_join(true, est, by=c("node_1", "node_2"))

  # Calculate correlation between estimates and true values.
  proportion_within_ci <- mean(comparison$true < comparison$est_ub & comparison$est > comparison$est_lb)
  expect_gt(proportion_within_ci, 0.9)
})
