test_that("binary conjugate posterior matches the analytic Beta posterior", {
  set.seed(1)
  df <- data.frame(n1 = c("A", "A", "A", "A"), n2 = c("B", "B", "B", "C"), event = c(1, 1, 1, 0), duration = 1)
  fit <- bison_model((event | duration) ~ dyad(n1, n2), df, "binary_conjugate",
                     priors = list(edge = "beta(1, 1)"), iter_sampling = 5000)
  expect_equal(dim(fit$chain), c(20000, 2))
  # A-B: 3 events in 3 periods -> Beta(4, 1); A-C: 0 in 1 -> Beta(1, 2)
  expect_equal(mean(plogis(fit$chain[, 1])), 4 / 5, tolerance = 0.01)
  expect_equal(mean(plogis(fit$chain[, 2])), 1 / 3, tolerance = 0.01)
})

test_that("aggregated and disaggregated data give the same conjugate posterior", {
  set.seed(2)
  sim <- simulate_bison_model("binary", aggregated = FALSE)
  df <- sim$df_sim
  agg <- aggregate(cbind(event, duration) ~ node_1_id + node_2_id, df, sum)
  f <- (event | duration) ~ dyad(node_1_id, node_2_id)
  fit_disagg <- bison_model(f, df, "binary_conjugate", priors = conjugate_priors("binary"))
  fit_agg <- bison_model(f, agg, "binary_conjugate", priors = conjugate_priors("binary"))
  e1 <- get_edgelist(fit_disagg)
  e2 <- get_edgelist(fit_agg)
  merged <- merge(e1, e2, by = c("node_1", "node_2"))
  expect_equal(nrow(merged), fit_agg$num_dyads)
  expect_equal(merged$median.x, merged$median.y, tolerance = 0.03)
})

test_that("count conjugate posterior matches the analytic Gamma posterior", {
  set.seed(3)
  df <- data.frame(n1 = c("A", "A", "B"), n2 = c("B", "B", "C"), event = c(4, 6, 0), duration = c(2, 3, 5))
  fit <- bison_model((event | duration) ~ dyad(n1, n2), df, "count_conjugate",
                     priors = list(edge = "gamma(1, 1)"), iter_sampling = 5000)
  # A-B: Gamma(1 + 10, 1 + 5), B-C: Gamma(1, 6)
  expect_equal(mean(exp(fit$chain[, 1])), 11 / 6, tolerance = 0.02)
  expect_equal(mean(exp(fit$chain[, 2])), 1 / 6, tolerance = 0.02)
})

test_that("conjugate predictions and log-likelihoods are per data row", {
  sim <- simulate_bison_model("count", aggregated = FALSE)
  fit <- bison_model((event | duration) ~ dyad(node_1_id, node_2_id), sim$df_sim, "count_conjugate",
                     priors = conjugate_priors("count"), iter_sampling = 100)
  expect_equal(dim(fit$event_preds), c(400, nrow(sim$df_sim)))
  expect_equal(dim(fit$log_lik), c(400, nrow(sim$df_sim)))
  expect_true(all(is.finite(fit$log_lik)))
})

test_that("priors_only conjugate draws come from the prior", {
  set.seed(4)
  sim <- simulate_bison_model("binary", aggregated = TRUE)
  fit <- bison_model((event | duration) ~ dyad(node_1_id, node_2_id), sim$df_sim, "binary_conjugate",
                     priors = list(edge = "beta(2, 8)"), priors_only = TRUE)
  expect_equal(mean(plogis(fit$chain)), 0.2, tolerance = 0.01)
})

test_that("conjugate models reject terms they can't fit", {
  sim <- simulate_bison_model("binary", aggregated = TRUE)
  expect_error(
    bison_model((event | duration) ~ dyad(node_1_id, node_2_id) + age_diff, sim$df_sim,
                "binary_conjugate", priors = conjugate_priors("binary")),
    "only support dyad"
  )
  expect_error(
    bison_model((event | duration) ~ 1, sim$df_sim, "binary_conjugate", priors = conjugate_priors("binary")),
    "need a dyad"
  )
})

test_that("conjugate models work with summaries, plots and model comparison", {
  sim <- simulate_bison_model("binary", aggregated = TRUE)
  fit <- bison_model((event | duration) ~ dyad(node_1_id, node_2_id), sim$df_sim, "binary_conjugate",
                     priors = conjugate_priors("binary"))
  expect_output(print(summary(fit)), "Edge list summary")
  expect_output(print(fit), "binary_conjugate")
  pdf(NULL)
  on.exit(dev.off())
  expect_silent(plot_predictions(fit, type = c("density", "stats", "point")))
  expect_silent(prior_predictive_check((event | duration) ~ dyad(node_1_id, node_2_id), sim$df_sim,
                                       "binary_conjugate", priors = conjugate_priors("binary")))
  expect_error(plot_trace(fit), "trace plots aren't needed")
  comparison <- suppressWarnings(model_comparison(list(a = fit, b = fit)))
  expect_equal(names(comparison), c("a", "b"))
})

test_that("model_comparison works when models have different numbers of draws", {
  sim <- simulate_bison_model("binary", aggregated = TRUE)
  f <- (event | duration) ~ dyad(node_1_id, node_2_id)
  fit_a <- bison_model(f, sim$df_sim, "binary_conjugate", priors = conjugate_priors("binary"), iter_sampling = 100)
  fit_b <- bison_model(f, sim$df_sim, "binary_conjugate", priors = conjugate_priors("binary"), iter_sampling = 75)
  expect_message(
    comparison <- suppressWarnings(model_comparison(list(a = fit_a, b = fit_b))),
    "different numbers of posterior draws"
  )
  expect_equal(names(comparison), c("a", "b"))

  sim_other <- simulate_bison_model("binary", aggregated = TRUE, num_nodes = 6)
  fit_c <- bison_model(f, sim_other$df_sim, "binary_conjugate", priors = conjugate_priors("binary"))
  expect_error(model_comparison(list(a = fit_a, c = fit_c)), "same data")
  expect_error(model_comparison(list(a = fit_a)), "at least two")
})
