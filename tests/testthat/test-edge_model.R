test_that("Binary edge model parameter estimation", {
  library(dplyr)
  library(igraph)

  set.seed(123)

  # Load data in with minimal effects
  sim_data <- simulate_edge_model("binary", aggregated = TRUE, location_effect = FALSE, age_diff_effect = FALSE)
  df <- sim_data$df_sim

  # Fit model to the data
  fit_edge <- edge_model(
    (event | duration) ~ dyad(node_1_id, node_2_id),
    data=df,
    data_type="binary"
  )

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
  expect_warning(plot_trace(fit_edge, par_ids=1), regexp=NA)

  # Modify comparison dataframe to test dyadic regression.
  df_dyadic <- comparison
  dyadic_trait <- rnorm(nrow(df_dyadic), df_dyadic$true)
  df_dyadic$dyadic_trait <- dyadic_trait
  dyadic_coef_true <- lm(true ~ dyadic_trait, df_dyadic)$coefficient[[2]]

  fit_dyadic <- dyadic_regression(dyad(node_1, node_2) ~ dyadic_trait, fit_edge, df_dyadic, mm=FALSE)

  summary_dyadic <- summary(fit_dyadic)
  expect_equal(summary_dyadic[1, 2] < 0 & summary_dyadic[1, 3] > 0, TRUE)
  expect_equal(summary_dyadic[2, 2] < dyadic_coef_true & summary_dyadic[2, 3] > dyadic_coef_true, TRUE)

  # Check that plots don't produce warnings
  expect_warning(plot_predictions(fit_dyadic), regexp=NA)
  expect_warning(plot_trace(fit_dyadic, par_ids=1), regexp=NA)

  net <- graph_from_edgelist(as.matrix(comparison[, 1:2]), directed=FALSE)
  E(net)$weight <- plogis(comparison[, 3])
  nodal_metric <- strength(net)
  nodal_metric <- nodal_metric - mean(nodal_metric)
  df_nodal <- data.frame(node=factor(V(net), levels=1:length(V(net))), metric_true=nodal_metric)
  df_nodal$nodal_trait <- rnorm(nrow(df_nodal), 2 * df_nodal$metric)
  nodal_coef_true <- lm(metric_true ~ nodal_trait, df_nodal)$coefficient[[2]]

  fit_nodal <- nodal_regression(strength(node) ~ nodal_trait, fit_edge, df_nodal)
  summary_nodal <- summary(fit_nodal)
  expect_equal(summary_nodal[1, 2] <= 0 & summary_nodal[1, 3] >= 0, TRUE)
  expect_equal(summary_nodal[2, 2] <= nodal_coef_true & summary_nodal[2, 3] >= nodal_coef_true, TRUE)

  # Check that plots don't produce warnings
  expect_warning(plot_predictions(fit_nodal), regexp=NA)
  expect_warning(plot_trace(fit_nodal, par_ids=1), regexp=NA)
})

test_that("Count edge model parameter estimation", {
  library(dplyr)

  set.seed(123)

  # Load data in with minimal effects
  sim_data <- simulate_edge_model("count", aggregated = TRUE, location_effect = FALSE, age_diff_effect = FALSE)
  df <- sim_data$df_sim

  # Fit model to the data
  fit_edge <- edge_model(
    (event | duration) ~ dyad(node_1_id, node_2_id),
    data=df,
    data_type="count"
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

  # Modify comparison dataframe to test dyadic regression.
  df_dyadic <- comparison
  dyadic_trait <- rnorm(nrow(df_dyadic), df_dyadic$true)
  df_dyadic$dyadic_trait <- dyadic_trait
  dyadic_coef_true <- lm(true ~ dyadic_trait, df_dyadic)$coefficient[[2]]

  fit_dyadic <- dyadic_regression(dyad(node_1, node_2) ~ dyadic_trait, fit_edge, df_dyadic, mm=FALSE)

  summary_dyadic <- summary(fit_dyadic)
  expect_equal(summary_dyadic[1, 2] <= 0 & summary_dyadic[1, 3] >= 0, TRUE)
  expect_equal(summary_dyadic[2, 2] <= dyadic_coef_true & summary_dyadic[2, 3] >= dyadic_coef_true, TRUE)
})
