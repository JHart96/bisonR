test_that("bison_brm fits dyadic, nodal and global regressions", {
  skip_if_not_installed("brms")
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
