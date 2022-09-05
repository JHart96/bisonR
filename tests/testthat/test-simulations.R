test_that("duration simulation works", {
  expect_warning(
    simulate_bison_model("binary", aggregated = FALSE, location_effect = TRUE, age_diff_effect = TRUE),
    regexp=NA
  )

  expect_warning(
    simulate_bison_model("binary", aggregated = TRUE, location_effect = TRUE, age_diff_effect = TRUE),
    regexp=NA
  )

  expect_warning(
    simulate_bison_model("count", aggregated = FALSE, location_effect = TRUE, age_diff_effect = TRUE),
    regexp=NA
  )

  expect_warning(
    simulate_bison_model("duration", aggregated = FALSE, location_effect = TRUE, age_diff_effect = TRUE),
    regexp=NA
  )

  expect_warning(
    simulate_bison_model_mixture("binary", num_components=2, component_weights=c(0.5, 0.5)),
    regexp=NA
  )

  expect_warning(
    simulate_bison_model_mixture("count", num_components=2, component_weights=c(0.5, 0.5)),
    regexp=NA
  )
})
