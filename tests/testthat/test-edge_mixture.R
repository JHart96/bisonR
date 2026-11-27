for (model_type in c("binary", "count")) {
  test_that(paste(model_type, "edge mixtures find well-separated components"), {
    set.seed(20)
    sim <- simulate_bison_model_mixture(model_type, num_components = 2, component_weights = c(0.5, 0.5))
    fit_edge <- bison_model((event | duration) ~ dyad(node_1_id, node_2_id), sim$df_sim,
                            paste0(model_type, "_conjugate"), priors = conjugate_priors(model_type))
    fit_mixture <- bison_mixture(fit_edge, num_components = 3, verbose = FALSE, num_draws = 50)

    probs <- get_network_component_probabilities(fit_mixture)
    expect_equal(sum(probs$probability), 1)
    if (model_type == "count") {
      # Binary edge posteriors have long tails on the logit scale, which makes the
      # number of components less identifiable, so only check recovery for counts.
      expect_equal(which.max(probs$probability), 2)
    }

    # Components are ordered by their means, in every draw
    means <- fit_mixture$component_mean_samples[[2]]
    expect_true(all(means[, 1] < means[, 2], na.rm = TRUE))
    component_means <- get_component_means(fit_mixture, 2)
    expect_lt(component_means[1, 1], component_means[2, 1])

    edge_probs <- get_edge_component_probabilities(fit_mixture, 2)
    expect_equal(colnames(edge_probs), c("node_1", "node_2", "P(component 1)", "P(component 2)"))
    expect_equal(unname(rowSums(edge_probs[, 3:4])), rep(1, nrow(edge_probs)))

    expect_output(print(summary(fit_mixture)), "Maximum number of components: 3")
    expect_error(get_edge_component_probabilities(fit_mixture, 4), "between 1 and 3")
    expect_error(get_edge_component_probabilities(fit_mixture, 0), "between 1 and 3")
    expect_error(get_component_means(fit_mixture, 4), "between 1 and 3")
  })
}

test_that("plot_mixture draws the components and returns the assignments", {
  set.seed(21)
  sim <- simulate_bison_model_mixture("count", num_components = 2, component_weights = c(0.5, 0.5))
  fit_edge <- bison_model((event | effort) ~ dyad(node_1_id, node_2_id), sim$df_sim,
                          "count_conjugate", priors = conjugate_priors("count"))
  fit_mixture <- bison_mixture(fit_edge, num_components = 3, verbose = FALSE, num_draws = 30)
  pdf(NULL)
  on.exit(dev.off())
  out <- plot_mixture(fit_mixture, num_components = 2)
  expect_equal(nrow(out), fit_edge$num_dyads)
  expect_true(all(out$component %in% 1:2))
  expect_silent(plot_mixture(fit_mixture))
  expect_error(plot_mixture(fit_mixture, num_components = 4), "between 1 and 3")
})

test_that("the mixture criterion is not case sensitive and draws are recorded", {
  set.seed(22)
  sim <- simulate_bison_model_mixture("count", num_components = 2, component_weights = c(0.5, 0.5))
  fit_edge <- bison_model((event | effort) ~ dyad(node_1_id, node_2_id), sim$df_sim,
                          "count_conjugate", priors = conjugate_priors("count"))
  fit_mixture <- bison_mixture(fit_edge, num_components = 2, criterion = "icl", verbose = FALSE, num_draws = 10)
  expect_equal(fit_mixture$criterion, "ICL")
  expect_length(fit_mixture$draw_ids, 10)
  expect_equal(rownames(get_component_means(fit_mixture, 2)), c("Component 1", "Component 2"))
})
