test_that("default priors can be plotted for every model type", {
  pdf(NULL)
  on.exit(dev.off())
  for (model_type in c("binary", "count", "binary_conjugate", "count_conjugate")) {
    expect_silent(prior_check(get_default_priors(model_type), model_type))
  }
})

test_that("priors are converted to Stan data", {
  params <- extract_prior_parameters(list(edge = "normal(-1, 2)", random_std = "half-normal(0.5)", zero_prob = "beta(1, 3)"))
  expect_equal(params$prior_edge_mu, -1)
  expect_equal(params$prior_edge_sigma, 2)
  expect_equal(params$prior_random_std_sigma, 0.5)
  expect_equal(params$prior_zero_prob_beta, 3)
})

test_that("prior densities for edge weights integrate to one on the natural scale", {
  binary <- get_density_fn("edge", "binary", "normal")
  expect_equal(integrate(function(x) binary(x, 0, 2.5), 0, 1)$value, 1, tolerance = 1e-4)
  # normal(0, 2.5) on the logit scale puts more density near 0 and 1 than at 0.5
  expect_gt(binary(0.02, 0, 2.5), binary(0.5, 0, 2.5))
  count <- get_density_fn("edge", "count", "normal")
  expect_equal(integrate(function(x) count(x, 0, 1), 0, Inf)$value, 1, tolerance = 1e-4)
})
