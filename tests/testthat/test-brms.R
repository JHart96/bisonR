test_that("brms works", {
  sim_data <- simulate_bison_model("binary", aggregated = TRUE)
  df <- sim_data$df_sim
  levels(df$node_1_id) <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")
  levels(df$node_2_id) <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")

  fit_edge <- bison_model(
    (event | duration) ~ dyad(node_1_id, node_2_id),
    data=df,
    data_type="binary_conjugate",
    priors=get_default_priors("binary_conjugate")
  )

  expect_error (
    suppressWarnings(bison_brm(
      age_diff ~ bison(edge_weight(node_1_id, node_2_id)),
      fit_edge,
      df,
      num_draws=5,
      silent=2,
      refresh=0
    )),
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
      refresh=0
    )),
    regexp=NA
  )

  expect_error(
    suppressWarnings(bison_brm(
      bison(node_strength(node)) ~ trait,
      fit_edge,
      df_nodal,
      num_draws=5,
      silent=2,
      refresh=0
    )),
    regexp=NA
  )
})
