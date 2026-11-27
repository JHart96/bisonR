utils::globalVariables(c("node_1_id", "node_2_id"))
utils::globalVariables(c("num_locations", "age_diff", "locations", "age_1", "age_2"))

#' Simulate observation data for an edge model
#'
#' Simulates observations of social events in an undirected network with known
#' edge weights. Use it to try out bisonR functions, to learn the format of the
#' data, or to check that a model recovers the true values before using it on
#' real data.
#'
#' @details Each dyad's true edge weight is drawn from a standard normal
#' distribution on the link scale (log-odds for binary data, log rates for
#' count data). Each node gets an age, and each observation takes place at one
#' of five locations, each with its own effect. Binary data have one sampling
#' period (effort 1) per observation; count data have an observation effort
#' drawn uniformly between 1 and `max_obs`.
#'
#' @param model_type Type of data to simulate, `"binary"` (whether an event
#'   occurred in each sampling period) or `"count"` (the number of events).
#' @param aggregated If `TRUE`, observations are summed to one row per dyad,
#'   with the location column removed. If `FALSE`, there is one row per
#'   observation.
#' @param location_effect Whether location affects the probability or rate of
#'   events.
#' @param age_diff_effect Whether the age difference between the two nodes
#'   affects the probability or rate of events (by 0.2 per unit of age
#'   difference on the link scale).
#' @param num_nodes Number of nodes in the simulated network.
#' @param max_obs Maximum number of observations per dyad. Each dyad gets between 1 and `max_obs`.
#'
#' @return A list with `df_sim`, the simulated observations, and `df_true`, the
#'   true edge weights on the link scale. `df_sim` has columns `node_1_id`,
#'   `node_2_id`, `event` and `effort` (observation effort), and predictors. It
#'   also has a `duration` column, identical to `effort`, for code written for
#'   earlier versions; new code should use `effort`. The true edge weights in
#'   `df_true` don't include the effects of age difference or location.
#' @seealso [simulate_bison_model_mixture()], [bison_model()]
#' @examples
#' sim_data <- simulate_bison_model("binary", aggregated = FALSE)
#' head(sim_data$df_sim)  # Observations
#' head(sim_data$df_true) # True edge weights (on the logit scale)
#' @export
simulate_bison_model <- function(model_type, aggregated, location_effect=TRUE, age_diff_effect=TRUE,
                                 num_nodes = 10, max_obs = 5) {
  if (!is.character(model_type) || length(model_type) != 1 || !model_type %in% c("binary", "count")) {
    stop('`model_type` must be "binary" or "count".', call. = FALSE)
  }
  num_locations <- 5
  edge_weights <- matrix(rnorm(num_nodes^2, 0, 1), num_nodes, num_nodes)
  edge_weights <- edge_weights * upper.tri(edge_weights)
  ages <- rnorm(num_nodes, 20, 3)
  locations <- rnorm(num_locations, 0, 1)

  if (model_type %in% c("binary", "count")) {
    df_sim <- data.frame(event=numeric(), node_1_id=numeric(), node_2_id=numeric(), age_diff=numeric(), age_1=numeric(), age_2=numeric(), location=numeric(), duration=numeric())
    df_true <- data.frame(node_1_id=numeric(), node_2_id=numeric(), edge_weight=numeric(), age_diff=numeric(), age_1=numeric(), age_2=numeric())
  }

  for (i in 1:num_nodes) {
    for (j in 1:num_nodes) {
      if (i < j) {
        # Generate simulated dataframe
        age_diff <- ages[i] - ages[j]
        for (k in 1:sample.int(max_obs, 1)) {
          location_id <- sample.int(num_locations, 1)
          predictor <- edge_weights[i, j]
          if (age_diff_effect) {
            predictor <- predictor + 0.2 * age_diff
          }
          if (location_effect) {
            predictor <- predictor + locations[location_id]
          }
          if (model_type == "binary") {
            event <- rbinom(1, 1, plogis(predictor))
            df_sim <- append_sim_row(df_sim, list(event=event, node_1_id=i, node_2_id=j, age_diff=age_diff, age_1=ages[i], age_2=ages[j], location=location_id, duration=1))
          }
          if (model_type == "count") {
            duration <- runif(1, 1, max_obs)
            event <- rpois(1, exp(predictor) * duration)
            df_sim <- append_sim_row(df_sim, list(event=event, node_1_id=i, node_2_id=j, age_diff=age_diff, age_1=ages[i], age_2=ages[j], location=location_id, duration=duration))
          }
        }

        # Set true dataframe
        df_true <- append_sim_row(df_true, list(node_1_id=i, node_2_id=j, edge_weight=edge_weights[i, j], age_1=ages[i], age_2=ages[j], age_diff=age_diff))
      }
    }
  }
  if (aggregated) {
    df_sim <- dplyr::summarise(dplyr::group_by(df_sim, node_1_id, node_2_id), event=sum(event), duration=sum(duration), age_diff=mean(age_diff), age_1=mean(age_1), age_2=mean(age_2), .groups = "drop")
  }
  df_sim$node_1_id <- factor(df_sim$node_1_id, levels=1:num_nodes)
  df_sim$node_2_id <- factor(df_sim$node_2_id, levels=1:num_nodes)
  df_true$node_1_id <- factor(df_true$node_1_id, levels=1:num_nodes)
  df_true$node_2_id <- factor(df_true$node_2_id, levels=1:num_nodes)
  if (model_type %in% c("binary", "count")) {

    return(list(df_sim=add_effort_column(df_sim), df_true=df_true))
  }
}

#' Simulate observation data with several types of edge
#'
#' Simulates observations of social events in an undirected network of 15
#' nodes, where each edge belongs to one of several components (types of
#' relationship) with different mean edge weights. Use it to try out
#' [bison_mixture()], or to check that it recovers the true number of
#' components.
#'
#' @details Component means are evenly spaced between -2.5 and 2.5 on the link
#' scale (log-odds for binary data, log rates for count data), and edge weights
#' within a component have a standard deviation of 0.5. Each dyad is observed
#' between 1 and 20 times, and the data are summed to one row per dyad.
#'
#' @param model_type Type of data to simulate, `"binary"` or `"count"`.
#' @param num_components Number of components in the simulation (at least 2).
#' @param component_weights Vector giving the proportion of edges in each
#'   component, e.g. `c(0.5, 0.5)`. Must sum to 1.
#' @return A list with `df_sim`, the simulated observations, and `df_true`, the
#'   true edge weights on the link scale. `df_sim` has columns `node_1_id`,
#'   `node_2_id`, `event` and `effort` (observation effort). It
#'   also has a `duration` column, identical to `effort`, for code written for
#'   earlier versions; new code should use `effort`.
#' @seealso [bison_mixture()], [simulate_bison_model()]
#' @examples
#' sim_data <- simulate_bison_model_mixture(
#'   "binary", num_components = 2, component_weights = c(0.5, 0.5)
#' )
#' head(sim_data$df_sim)
#' @export
simulate_bison_model_mixture <- function(model_type, num_components, component_weights) {
  num_nodes <- 15
  num_edges <- 0.5 * num_nodes * (num_nodes - 1)
  max_obs <- 20
  component_ids <- extraDistr::rcat(num_edges, component_weights)
  component_mu <- seq(-2.5, 2.5, 5/(num_components - 1))
  component_sigma <- rep(0.5, num_components)
  edge_weights <- matrix(0, num_nodes, num_nodes)
  edge_weights[upper.tri(edge_weights)] <- rnorm(num_edges, mean=component_mu[component_ids], sd=component_sigma[component_ids])

  df_sim <- data.frame(event=numeric(), node_1_id=numeric(), node_2_id=numeric(), duration=numeric())
  df_true <- data.frame(node_1_id=numeric(), node_2_id=numeric(), edge_weight=numeric())

  for (i in 1:num_nodes) {
    for (j in 1:num_nodes) {
      if (i < j) {
        for (k in 1:sample.int(max_obs, 1)) {
          predictor <- edge_weights[i, j]
          if (model_type == "binary") {
            event <- rbinom(1, 1, plogis(predictor))
            df_sim <- append_sim_row(df_sim, list(event=event, node_1_id=i, node_2_id=j, duration=1))
          }
          if (model_type == "count") {
            duration <- runif(1, 1, max_obs)
            event <- rpois(1, exp(predictor) * duration)
            df_sim <- append_sim_row(df_sim, list(event=event, node_1_id=i, node_2_id=j, duration=duration))
          }
        }

        # Set true dataframe
        df_true <- append_sim_row(df_true, list(node_1_id=i, node_2_id=j, edge_weight=edge_weights[i, j]))
      }
    }
  }
  df_sim <- dplyr::summarise(dplyr::group_by(df_sim, node_1_id, node_2_id), event=sum(event), duration=sum(duration), .groups = "drop")
  df_sim$node_1_id <- factor(df_sim$node_1_id, levels=1:num_nodes)
  df_sim$node_2_id <- factor(df_sim$node_2_id, levels=1:num_nodes)
  df_true$node_1_id <- factor(df_true$node_1_id, levels=1:num_nodes)
  df_true$node_2_id <- factor(df_true$node_2_id, levels=1:num_nodes)
  return(list(df_sim=add_effort_column(df_sim), df_true=df_true))
}


# Add an `effort` column (observation effort) and keep `duration` as a copy for
# code written for earlier versions.
add_effort_column <- function(df_sim) {
  df_sim <- as.data.frame(df_sim)
  df_sim$effort <- df_sim$duration
  other <- setdiff(colnames(df_sim), c("event", "effort", "duration"))
  df_sim[, c(other[other %in% c("node_1_id", "node_2_id")], "event", "effort",
             other[!other %in% c("node_1_id", "node_2_id")], "duration")]
}

# Add one row to a simulated data frame, matching values to columns by name
# (assigning a list by position would put values in the wrong columns if the
# order differed).
append_sim_row <- function(df, row) {
  df[nrow(df) + 1, names(row)] <- row
  df
}
