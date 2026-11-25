test_that("simulations return data in the expected format", {
  for (model_type in c("binary", "count")) {
    for (aggregated in c(TRUE, FALSE)) {
      sim <- simulate_bison_model(model_type, aggregated = aggregated)
      expect_named(sim, c("df_sim", "df_true"))
      expect_true(all(c("event", "node_1_id", "node_2_id", "duration") %in% colnames(sim$df_sim)))
      expect_equal(nrow(sim$df_true), 45)
      expect_true(all(sim$df_sim$event >= 0))
    }
  }
})

test_that("binary simulations never have more events than sampling periods", {
  sim <- simulate_bison_model("binary", aggregated = TRUE)
  expect_true(all(sim$df_sim$event <= sim$df_sim$duration))
})

test_that("mixture simulations work", {
  for (model_type in c("binary", "count")) {
    sim <- simulate_bison_model_mixture(model_type, num_components = 2, component_weights = c(0.5, 0.5))
    expect_equal(nrow(sim$df_true), 105)
  }
})

test_that("simulated data have an effort column, with duration kept for older code", {
  sim <- simulate_bison_model("count", aggregated = TRUE)
  expect_equal(sim$df_sim$effort, sim$df_sim$duration)
  expect_equal(colnames(sim$df_sim)[1:4], c("node_1_id", "node_2_id", "event", "effort"))
})
