test_that("duration simulation works", {
  expect_warning(
    simulate_edge_model("duration", aggregated = FALSE, location_effect = TRUE, age_diff_effect = TRUE),
    regexp=NA
  )
})
