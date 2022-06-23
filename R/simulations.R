#' Simulates observations from an edge model
#'
#' @return
#' @export
#'
#' @examples
simulate_edge_model <- function(model_type, aggregated, location_effect=TRUE, age_diff_effect=TRUE) {
  num_nodes <- 10
  num_locations <- 5
  max_obs <- 20
  edge_weights <- matrix(rnorm(num_nodes^2, 0, 1), num_nodes, num_nodes)
  edge_weights <- edge_weights * upper.tri(edge_weights)
  ages <- rnorm(num_nodes, 20, 3)
  locations <- rnorm(num_locations, 0, 1)

  if (model_type %in% c("binary", "count")) {
    df_sim <- data.frame(event=numeric(), node_1_id=numeric(), node_2_id=numeric(), age_diff=numeric(), age_1=numeric(), age_2=numeric(), location=numeric(), duration=numeric())
    df_true <- data.frame(node_1_id=numeric(), node_2_id=numeric(), edge_weight=numeric(), age_diff=numeric())
  }
  if (model_type == "duration") {
    df_sim <- data.frame(event=numeric(), node_1_id=numeric(), node_2_id=numeric(), age_diff=numeric(), location=numeric())
    df_sim_agg <- data.frame(event_count=numeric(), node_1_id=numeric(), node_2_id=numeric())
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
            predictor <- predictor + age_diff
          }
          if (location_effect) {
            predictor <- predictor + locations[location_id]
          }
          if (model_type == "binary") {
            event <- rbinom(1, 1, plogis(predictor))
            df_sim[nrow(df_sim) + 1, ] <- list(event=event, node_1_id=i, node_2_id=j, age_diff=age_diff, age_1=ages[i], age_2=ages[j], location=location_id, duration=1)
          }
          if (model_type == "count") {
            duration <- runif(1, 1, max_obs)
            event <- rpois(1, exp(predictor) * duration)
            df_sim[nrow(df_sim) + 1, ] <- list(event=event, node_1_id=i, node_2_id=j, age_diff=age_diff, age_1=ages[i], age_2=ages[j], location=location_id, duration=duration)
          }
        }
        if (model_type == "duration") {
          # Draw a lambda, sample K, and loop through K
          duration <- runif(1, 100, 1000) # Hours observed
          lambda <- runif(1, 0, 0.1) # Events per hour
          K <- rpois(1, lambda * duration)
          for (k in 1:K) {
            location_id <- sample.int(num_locations, 1)
            predictor <- -5 + edge_weights[i, j] + 0.25 * age_diff + locations[location_id]
            event <- min(rexp(1, lambda/plogis(predictor)), 1)
            df_sim[nrow(df_sim) + 1, ] <- list(event=event, node_1_id=i, node_2_id=j, age_diff=age_diff, location=location_id)
          }
          df_sim_agg[nrow(df_sim_agg) + 1, ] <- list(event_count=K, node_1_id=i, node_2_id=j)
        }

        # Set true dataframe
        df_true[nrow(df_true) + 1, ] <- list(node_1_id=i, node_2_id=j, edge_weight=edge_weights[i, j], age_diff=age_diff)
      }
    }
  }
  if (aggregated) {
    df_sim <- dplyr::summarise(dplyr::group_by(df_sim, node_1_id, node_2_id), event=sum(event), duration=sum(duration))
  }
  df_sim$node_1_id <- factor(df_sim$node_1_id, levels=1:num_nodes)
  df_sim$node_2_id <- factor(df_sim$node_2_id, levels=1:num_nodes)
  df_true$node_1_id <- factor(df_true$node_1_id, levels=1:num_nodes)
  df_true$node_2_id <- factor(df_true$node_2_id, levels=1:num_nodes)
  if (model_type %in% c("binary", "count")) {
    return(list(df_sim=df_sim, df_true=df_true))
  }
  if (model_type == "duration") {
    df_sim_agg$node_1_id <- factor(df_sim_agg$node_1_id, levels=1:num_nodes)
    df_sim_agg$node_2_id <- factor(df_sim_agg$node_2_id, levels=1:num_nodes)
    return(list(df=df_sim, df_agg=df_sim_agg))
  }

  stop("Model type not supported")
}
