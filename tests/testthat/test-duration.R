test_that("multiplication works", {
  sim_data <- simulate_bison_model("duration", aggregated = TRUE, location_effect = FALSE, age_diff_effect = FALSE)
  df <- sim_data$df
  duration_data <- sim_data$df_agg

  priors <- get_default_priors("duration")
  # priors$rate <- "half-normal(0.05)"
  # priors$edge <- "normal(-5, 5)"

  # The model seems obsessed with getting estimates close to 1.0. What is going on?

  fit_bison <- bison_model(
    (event | duration) ~ dyad(node_1_id, node_2_id),
    data=df,
    model_type="duration",
    duration_data=duration_data,
    priors=priors,
    directed=TRUE
  )

  plot_predictions(fit_bison, type="point")

  # Why does true edge weight not correlate with SRI? It should.

  plot(
    plogis(sim_data$df_true$edge_weight),
    df$event/df$duration
  )
  abline(a=0, b=1)

  plot(
    fit_bison$fit$summary("edge_weight")[, "mean"][[1]],
    sim_data$df_true$edge_weight
  )
  abline(a=0, b=1)

  min(fit_bison$fit$summary("rate")[, "mean"][[1]])

})
