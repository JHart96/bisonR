test_that("brms works", {
  sim_data <- bisonR::simulate_bison_model("binary", aggregated = TRUE)
  df_sim <- sim_data$df_sim
  levels(df_sim$node_1_id) <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")
  levels(df_sim$node_2_id) <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")

  fit_edge <- bison_model(
    (event | duration) ~ dyad(node_1_id, node_2_id),
    data=df_sim,
    model_type="binary_conjugate",
    priors=get_default_priors("binary_conjugate")
  )

  fit_brm <- expect_error (
    suppressWarnings(bison_brm(
      age_diff ~ bison(edge_weight(node_1_id, node_2_id)),
      list(fit_edge, fit_edge),
      list(df_sim, df_sim),
      num_draws=5,
      silent=2,
      refresh=0
    )),
    regexp=NA
  )

  fit_brm <- expect_error (
    suppressWarnings(bison_brm(
      age_diff ~ bison(edge_weight(node_1_id, node_2_id)),
      list(fit_edge, fit_edge),
      list(df_sim, df_sim),
      num_draws=5,
      silent=2,
      refresh=0,
      sample_prior="only",
      prior=brms::prior(normal(0, 10), class="b")
    )),
    regexp=NA
  )

  priors <- expect_warning(bison_brm_get_prior(
      age_diff ~ bison(edge_weight(node_1_id, node_2_id)),
      list(fit_edge, fit_edge),
      list(df_sim, df_sim)
    ),
    regexp=NA
  )

  df_nodal <- data.frame(node=as.factor(c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")))
  df_nodal$trait <- rnorm(10)

  expect_error(
    suppressWarnings(bison_brm(
      trait ~ bison(node_eigen(node)),
      fit_edge,
      df_nodal,
      num_draws=5,
      silent=2,
      refresh=0,
      z_score=TRUE
    )),
    regexp=NA
  )

  fit_brm <- expect_error(
    suppressWarnings(bison_brm(
      trait ~ bison(node_eigen(node)) + bison_network,
      list(fit_edge, fit_edge),
      list(df_nodal, df_nodal),
      num_draws=5,
      silent=2,
      refresh=0
    )),
    regexp=NA
  )

  df_global <- data.frame(bison_network=as.factor(c("1", "2")), condition=c("before", "after"))

  fit_brm <- expect_error(
    suppressWarnings(bison_brm(
      bison(global_cv(bison_network)) ~ condition,
      list(fit_edge, fit_edge),
      df_global,
      num_draws=5,
      silent=2,
      refresh=0
    )),
    regexp=NA
  )
})
