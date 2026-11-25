test_that("bison_brm fits dyadic, nodal and global regressions", {
  skip_if_not_installed("brms")
  skip_if_brms_disabled()
  skip_if_not_installed("mice")
  skip_if_no_cmdstan()
  skip_on_cran()

  sim <- simulate_bison_model("binary", aggregated = TRUE)
  df_sim <- sim$df_sim
  levels(df_sim$node_1_id) <- levels(df_sim$node_2_id) <- LETTERS[1:10]
  fit_edge <- bison_model((event | duration) ~ dyad(node_1_id, node_2_id), df_sim,
                          "binary_conjugate", priors = conjugate_priors("binary"))

  expect_s3_class(
    bison_brm_get_prior(age_diff ~ bison(edge_weight(node_1_id, node_2_id)), fit_edge, df_sim),
    "brmsprior"
  )

  # A dyadic regression without a multi-membership term triggers a warning
  warnings <- character(0)
  withCallingHandlers(
    quiet(bison_brm(bison(edge_weight(node_1_id, node_2_id)) ~ age_diff, fit_edge, df_sim,
                    num_draws = 2, silent = 2, refresh = 0)),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("multi-membership", warnings)))
  fit_dyadic <- suppressWarnings(quiet(bison_brm(
    bison(edge_weight(node_1_id, node_2_id)) ~ age_diff + (1 | mm(node_1_id, node_2_id)),
    fit_edge, df_sim, num_draws = 2, silent = 2, refresh = 0
  )))
  expect_s3_class(fit_dyadic, "brmsfit")

  df_nodal <- data.frame(node = factor(LETTERS[1:10]), trait = rnorm(10))
  fit_nodal <- suppressWarnings(quiet(bison_brm(trait ~ bison(node_eigen(node)), fit_edge, df_nodal,
                                                 num_draws = 2, silent = 2, refresh = 0, z_score = TRUE)))
  expect_s3_class(fit_nodal, "brmsfit")

  fit_degree <- suppressWarnings(quiet(bison_brm(trait ~ bison(node_degree[0.5](node)), fit_edge, df_nodal,
                                                  num_draws = 2, silent = 2, refresh = 0)))
  expect_s3_class(fit_degree, "brmsfit")

  df_global <- data.frame(bison_network = factor(c("1", "2")), condition = c("before", "after"))
  fit_global <- suppressWarnings(quiet(bison_brm(bison(global_cv(bison_network)) ~ condition,
                                                  list(fit_edge, fit_edge), df_global,
                                                  num_draws = 2, silent = 2, refresh = 0)))
  expect_s3_class(fit_global, "brmsfit")
})

test_that("bison_brm gives clear errors for mismatched data", {
  skip_if_not_installed("brms")
  skip_if_brms_disabled()
  skip_if_not_installed("mice")
  sim <- simulate_bison_model("binary", aggregated = TRUE)
  fit_edge <- bison_model((event | duration) ~ dyad(node_1_id, node_2_id), sim$df_sim,
                          "binary_conjugate", priors = conjugate_priors("binary"))
  df_nodal <- data.frame(node = c("1", "2", "not_a_node"), trait = 1:3)
  expect_error(
    bison_mice(fit_edge, df_nodal, "node", "node", "strength", num_draws = 2),
    "not_a_node"
  )
  expect_error(parse_bison_brms_formula(y ~ x), "bison\\(")
  expect_error(parse_bison_brms_formula(y ~ bison(strength(node))), "edge_, node_ or global_")
})

test_that("subset selects rows of the regression data", {
  df <- data.frame(node = c("A", "B", "C"), sex = c("F", "M", "F"))
  expect_equal(subset_bison_data(df, quote(sex == "F"), environment(), "node")$node, c("A", "C"))
  expect_equal(nrow(subset_bison_data(list(df, df), quote(sex == "M"), environment(), "node")[[2]]), 1)
  expect_null(subset_bison_data(NULL, NULL, environment(), "node"))
  expect_error(subset_bison_data(df, quote(sex == "X"), environment(), "node"), "doesn't select any rows")
  expect_error(subset_bison_data(df, quote(1), environment(), "node"), "logical expression")
  expect_error(subset_bison_data(df, quote(sex == "F"), environment(), "global"), "global metrics")
})

test_that("network variables can be centred or standardised", {
  x <- c(1, 2, 3, 6)
  expect_equal(scale_draw(x, "none"), x)
  expect_equal(mean(scale_draw(x, "centre")), 0)
  expect_equal(sd(scale_draw(x, "z")), 1)
  expect_error(scale_draw(c(1, 1), "z"), "two distinct values")
})

test_that("bison_brm_get_prior accepts subset and scale", {
  skip_if_not_installed("brms")
  skip_if_brms_disabled()
  skip_if_not_installed("mice")
  sim <- simulate_bison_model("binary", aggregated = TRUE)
  fit_edge <- bison_model((event | effort) ~ dyad(node_1_id, node_2_id), sim$df_sim,
                          "binary_conjugate", priors = conjugate_priors("binary"))
  df_nodal <- data.frame(node = as.character(1:10), group = rep(c("a", "b"), 5), trait = rnorm(10))
  expect_s3_class(
    bison_brm_get_prior(trait ~ bison(node_strength(node)), fit_edge, df_nodal,
                        scale = "centre", subset = group == "a"),
    "brmsprior"
  )
})
