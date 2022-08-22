test_that("directed edge models work", {
  library(igraph)

  sim_data <- simulate_edge_model("binary", aggregated = FALSE)
  df <- sim_data$df_sim
  head(df)

  fit_edge <- expect_warning(
    edge_model (
      (event | duration) ~ dyad(node_1_id, node_2_id),
      data=df,
      data_type="binary",
      directed=T,
      priors=get_default_priors("binary")
    ),
    regexp=NA
  )

  # If a dyad ID is missing, don't fit an edge for it.
  # So we need to calculate dyad IDs based on which dyads are present.
  # Make dyad_to_idx an edgelist instead I think, and make it first.

  expect_warning(
    get_edgelist(fit_edge),
    regexp=NA
  )

  expect_warning(
    summary(fit_edge),
    regexp=NA
  )

  expect_output(
    print(summary(fit_edge))
  )

  expect_warning(
    plot_network(fit_edge, lwd=1),
    regexp=NA
  )

  df_dyadic <- get_edgelist(fit_edge, transform=FALSE)
  dyadic_trait <- rnorm(nrow(df_dyadic))
  df_dyadic$group_id <- as.factor(sample(1:4, size=nrow(df_dyadic), replace=TRUE))
  df_dyadic$dyadic_trait <- dyadic_trait
  expect_warning(
    fit_dyadic <- dyadic_regression(dyad(node_1, node_2) ~ dyadic_trait, fit_edge, df_dyadic, mm=FALSE),
    regexp=NA
  )

  edgelist <- get_edgelist(fit_edge, transform=FALSE)
  net <- igraph::graph_from_edgelist(as.matrix(edgelist[, 1:2]), directed=FALSE)
  igraph::E(net)$weight <- plogis(edgelist[, 3])
  nodal_metric <- igraph::strength(net)
  nodal_metric <- nodal_metric - mean(nodal_metric)
  df_nodal <- data.frame(node=factor(V(net), levels=1:length(V(net))), metric_true=nodal_metric)
  df_nodal$nodal_trait <- rnorm(nrow(df_nodal), 2 * df_nodal$metric)
  df_nodal$group_id <- sample(1:2, nrow(df_nodal), replace=TRUE)

  expect_warning(
    nodal_regression(strength(node) ~ nodal_trait, fit_edge, df_nodal),
    regexp=NA
  )

  fit_mixture <- expect_warning(
    edge_mixture(fit_edge, num_components=3, verbose=TRUE),
    regexp=NA
  )

  expect_warning(
    draw_network_metric_samples(fit_edge, "social_differentiation"),
    regexp=NA
  )

})
