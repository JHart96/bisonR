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

test_that("bison_mice matches regression rows to nodes and dyads by name, in any order", {
  skip_if_not_installed("mice")
  fit <- fake_bison_model(c("A", "B", "C", "D"), rbind(c("A", "B"), c("B", "C"), c("A", "C"), c("C", "D")),
                          c(0.5, 0.25, 0.1, 0.8))
  df_nodal <- data.frame(node = c("D", "A", "C", "B"))
  imputed <- mice::complete(bison_mice(fit, df_nodal, "node", "node", "strength", num_draws = 2), 1)
  expect_equal(imputed$bison_node_strength, c(0.8, 0.6, 1.15, 0.75))

  # Undirected dyads can be given in either order
  df_dyadic <- data.frame(n1 = c("D", "B"), n2 = c("C", "A"))
  imputed <- mice::complete(bison_mice(fit, df_dyadic, c("n1", "n2"), "edge", "weight", num_draws = 2), 1)
  expect_equal(imputed$bison_edge_weight, qlogis(c(0.8, 0.5)))
})

test_that("global metrics are scaled across networks and must match the edge models", {
  skip_if_not_installed("mice")
  fit_1 <- fake_bison_model(c("A", "B", "C"), rbind(c("A", "B"), c("B", "C")), c(0.2, 0.4))
  fit_2 <- fake_bison_model(c("A", "B", "C"), rbind(c("A", "B"), c("B", "C")), c(0.6, 0.8))
  df_global <- data.frame(bison_network = factor(1:2), condition = c("a", "b"))
  imputed <- mice::complete(
    bison_mice(list(fit_1, fit_2), df_global, "bison_network", "global", "density", num_draws = 2, scale = "z"), 1
  )
  expect_equal(imputed$bison_global_density, c(-1, 1) / sqrt(2))

  expect_error(
    bison_mice(list(fit_1, fit_2), df_global[1, ], "bison_network", "global", "density", num_draws = 2),
    "one row of data per edge model"
  )
  expect_error(
    bison_mice(list(fit_1, fit_2), df_global[2:1, ], "bison_network", "global", "density", num_draws = 2),
    "doesn't match the position"
  )
})

test_that("non-finite network variables and conflicting scale arguments are caught", {
  skip_if_not_installed("mice")
  fit <- fake_bison_model(c("A", "B", "C", "D"), rbind(c("A", "B"), c("B", "C")), c(0.5, 0.6))
  df_nodal <- data.frame(node = c("A", "B", "C", "D"))
  expect_error(
    suppressWarnings(bison_mice(fit, df_nodal, "node", "node", "closeness", num_draws = 2)),
    "isn't finite"
  )
  warnings <- character(0)
  withCallingHandlers(
    bison_mice(fit, df_nodal[1:3, , drop = FALSE], "node", "node", "strength",
               num_draws = 2, z_score = TRUE, scale = "centre"),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("overrides", warnings)))
  expect_error(parse_bison_brms_formula(y ~ (bison(node_strength(node)))), "without surrounding brackets")
})

test_that("scaling edge or node metrics within networks warns when bison_network is a predictor", {
  formula <- "bison_node_strength ~ bison_network"
  expect_warning(check_network_scaling("z", 2, "node", formula), "removes differences in its average")
  expect_silent(check_network_scaling("none", 2, "node", formula))
  expect_silent(check_network_scaling("z", 1, "node", formula))
  expect_silent(check_network_scaling("z", 2, "global", "bison_global_cv ~ bison_network"))
  expect_silent(check_network_scaling("z", 2, "node", "bison_node_strength ~ age"))
})
