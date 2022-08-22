test_that("prior checks work", {
  for (model_type in c("binary", "count", "binary_conjugate", "count_conjugate")) {
    expect_warning(
      prior_check(get_default_priors(model_type), model_type),
      regexp=NA
    )
  }
})
