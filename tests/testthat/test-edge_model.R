# These tests fit models with Stan, so they are slower and need CmdStan.

coverage <- function(true, edgelist) {
  mean(true >= edgelist[, 4] & true <= edgelist[, 5])
}

test_that("binary edge model recovers true edge weights", {
  skip_if_no_cmdstan()
  skip_on_cran()
  set.seed(123)
  sim <- simulate_bison_model("binary", aggregated = TRUE, location_effect = FALSE, age_diff_effect = FALSE)
  fit <- quiet(bison_model(
    (event | duration) ~ dyad(node_1_id, node_2_id), data = sim$df_sim,
    model_type = "binary", priors = get_default_priors("binary"), seed = 1
  ))
  est <- get_edgelist(fit, transform = FALSE)
  comparison <- merge(sim$df_true, est, by.x = c("node_1_id", "node_2_id"), by.y = c("node_1", "node_2"))
  expect_equal(nrow(comparison), 45)
  expect_gt(coverage(comparison$edge_weight, comparison[, c("node_1_id", "node_2_id", "median", "5%", "95%")]), 0.8)
  expect_gt(cor(comparison$edge_weight, comparison$median), 0.5)
})

test_that("count edge model recovers true edge weights", {
  skip_if_no_cmdstan()
  skip_on_cran()
  set.seed(124)
  sim <- simulate_bison_model("count", aggregated = TRUE, location_effect = FALSE, age_diff_effect = FALSE)
  fit <- quiet(bison_model(
    (event | duration) ~ dyad(node_1_id, node_2_id), data = sim$df_sim,
    model_type = "count", priors = get_default_priors("count"), seed = 1
  ))
  est <- get_edgelist(fit, transform = FALSE)
  comparison <- merge(sim$df_true, est, by.x = c("node_1_id", "node_2_id"), by.y = c("node_1", "node_2"))
  expect_gt(coverage(comparison$edge_weight, comparison[, c("node_1_id", "node_2_id", "median", "5%", "95%")]), 0.8)
  expect_gt(cor(comparison$edge_weight, comparison$median), 0.8)
})

test_that("Stan and conjugate binary models agree", {
  skip_if_no_cmdstan()
  skip_on_cran()
  set.seed(5)
  sim <- simulate_bison_model("binary", aggregated = TRUE, location_effect = FALSE, age_diff_effect = FALSE)
  f <- (event | duration) ~ dyad(node_1_id, node_2_id)
  # beta(1, 1) on the probability scale is close to normal(0, 1.6) on the logit scale
  fit_stan <- quiet(bison_model(f, sim$df_sim, "binary", priors = list(edge = "normal(0, 1.6)"), seed = 1))
  fit_conj <- bison_model(f, sim$df_sim, "binary_conjugate", priors = list(edge = "beta(1, 1)"))
  expect_equal(get_edgelist(fit_stan)$median, get_edgelist(fit_conj)$median, tolerance = 0.1)
})

test_that("fixed effects ending in 1 are included in the model (issue #80)", {
  skip_if_no_cmdstan()
  skip_on_cran()
  sim <- simulate_bison_model("binary", aggregated = FALSE)
  df <- sim$df_sim
  df$age1 <- df$age_diff
  fit <- quiet(bison_model(
    (event | duration) ~ dyad(node_1_id, node_2_id) + age1 + (1 | location), data = df,
    model_type = "binary", priors = get_default_priors("binary"), iter_sampling = 200, iter_warmup = 200, seed = 1
  ))
  expect_equal(fit$model_data$num_fixed, 1)
  expect_equal(ncol(fit$fit$draws("beta_fixed", format = "matrix")), 1)
  expect_equal(ncol(fit$fit$draws("beta_random", format = "matrix")), 5)
})

test_that("zero-inflated and partially pooled models sample sensibly", {
  skip_if_no_cmdstan()
  skip_on_cran()
  set.seed(6)
  sim <- simulate_bison_model("binary", aggregated = TRUE)
  fit <- quiet(bison_model(
    (event | duration) ~ dyad(node_1_id, node_2_id), data = sim$df_sim,
    model_type = "binary", partial_pooling = TRUE, zero_inflated = TRUE,
    priors = get_default_priors("binary"), seed = 1
  ))
  zero_prob <- fit$fit$draws("zero_prob", format = "matrix")
  expect_true(all(zero_prob >= 0 & zero_prob <= 1))
  diagnostics <- fit$fit$diagnostic_summary(quiet = TRUE)
  expect_lt(sum(diagnostics$num_divergent), 40) # under 1% of draws
  expect_equal(dim(fit$event_preds), c(4000, nrow(sim$df_sim)))
})

test_that("models without dyads (null models) can be compared with edge models", {
  skip_if_no_cmdstan()
  skip_on_cran()
  set.seed(7)
  sim <- simulate_bison_model("binary", aggregated = FALSE)
  fit_edge <- quiet(bison_model((event | duration) ~ dyad(node_1_id, node_2_id), sim$df_sim, "binary",
                                priors = get_default_priors("binary"), seed = 1))
  fit_null <- quiet(bison_model((event | duration) ~ 1, sim$df_sim, "binary",
                                priors = get_default_priors("binary"), seed = 1))
  expect_error(get_edgelist(fit_null), "no edge weights")
  comparison <- suppressWarnings(model_comparison(list(edge = fit_edge, null = fit_null)))
  expect_equal(names(comparison), c("edge", "null"))
  expect_gt(comparison[["edge"]], 0.5)
  pdf(NULL)
  on.exit(dev.off())
  expect_silent(plot_trace(fit_edge, par_ids = 1:2))
  expect_silent(plot_predictions(fit_edge, type = c("density", "point")))
  expect_silent(plot_predictions(fit_null, type = "stats"))
})

test_that("directed models keep both directions and plot", {
  skip_if_no_cmdstan()
  skip_on_cran()
  df <- data.frame(
    n1 = c("A", "B", "A", "C", "B", "C"), n2 = c("B", "A", "C", "A", "C", "B"),
    event = c(5, 1, 3, 3, 0, 4), duration = 6
  )
  fit <- quiet(bison_model((event | duration) ~ dyad(n1, n2), df, "binary", directed = TRUE,
                           priors = get_default_priors("binary"), seed = 1))
  el <- get_edgelist(fit)
  expect_equal(nrow(el), 6)
  expect_gt(el$median[1], el$median[2]) # A -> B stronger than B -> A
  expect_output(print(summary(fit)), "A -> B")
  pdf(NULL)
  on.exit(dev.off())
  expect_silent(plot_network(fit))
})

test_that("models can be saved and loaded", {
  skip_if_no_cmdstan()
  skip_on_cran()
  sim <- simulate_bison_model("binary", aggregated = TRUE)
  fit <- quiet(bison_model((event | duration) ~ dyad(node_1_id, node_2_id), sim$df_sim, "binary",
                           priors = get_default_priors("binary"), iter_sampling = 100, iter_warmup = 100))
  path <- withr::local_tempfile(fileext = ".rds")
  bison_save(fit, path)
  reloaded <- bison_load(path)
  expect_equal(get_edgelist(reloaded), get_edgelist(fit))
  expect_equal(dim(reloaded$fit$draws("edge_weight")), dim(fit$fit$draws("edge_weight")))
})
