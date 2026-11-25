#' Generate a mice mids object
#'
#' @param edgemodel_list A bisonR edge model (or a list of them)
#' @param data_list A dataframe of regression variables compatible with a brms model (or a list of them).
#' @param num_draws Number of draws from the network posterior
#' @param param_names List of column names for network features (e.g. node_1_id, node_2_id or just node_id)
#' @param target_name Name of network feature (edge, node, global)
#' @param metric_name Name of network metric to calculate
#' @param z_score Whether to Z-score bison variable or not.
#'
#' @return A mice mids object
#' @examples
#' \donttest{
#' sim_data <- simulate_bison_model("binary", aggregated = TRUE)
#' df <- sim_data$df_sim
#' fit_edge <- bison_model(
#'   (event | duration) ~ dyad(node_1_id, node_2_id),
#'   data = df,
#'   model_type = "binary_conjugate",
#'   priors = get_default_priors("binary_conjugate")
#' )
#' # One imputed dataset per posterior draw of node strength
#' df_nodal <- data.frame(node = levels(df$node_1_id))
#' imputed <- bison_mice(fit_edge, df_nodal, "node", "node", "strength", num_draws = 5)
#' head(mice::complete(imputed, 1))
#' }
#' @export
bison_mice <- function(edgemodel_list, data_list, param_names, target_name, metric_name, num_draws=100, z_score=FALSE) {
  check_suggested("mice", "bison_mice")
  if (target_name == "global" && is.data.frame(data_list)) {
    data_list <- lapply(seq_len(nrow(data_list)), function(i) data_list[i, , drop = FALSE])
  }

  # If the data aren't lists yet, convert them
  if (inherits(edgemodel_list, "bison_model")) edgemodel_list <- list(edgemodel_list)
  if (is.data.frame(data_list)) data_list <- list(data_list)

  new_bison_term <- bison_term_name(target_name, metric_name)

  posterior_samples_list <- list()

  if (length(edgemodel_list) != length(data_list) && target_name != "global") {
    stop(
      "Got ", length(edgemodel_list), " edge model(s) but ", length(data_list),
      " data frame(s). Provide one data frame per edge model.",
      call. = FALSE
    )
  }

  for (i in seq_along(edgemodel_list)) {
    model <- edgemodel_list[[i]]
    metric_string <- paste0(target_name, "_", metric_name)
    if (target_name == "edge") {
      node_1_ids <- lookup_node_ids(model, data_list[[i]], param_names[1], i)
      node_2_ids <- lookup_node_ids(model, data_list[[i]], param_names[2], i)
      dyad_ids <- get_dyad_ids(node_1_ids, node_2_ids, model$dyad_to_idx, model$directed)
      if (anyNA(dyad_ids)) {
        bad <- which(is.na(dyad_ids))[1]
        stop(
          "Row ", bad, " of data frame ", i, " refers to a dyad that isn't in edge model ", i, ".",
          call. = FALSE
        )
      }
      if (metric_name != "weight") {
        stop("Only edge_weight is available as an edge-level metric in bison_brm().", call. = FALSE)
      }
      # Edge weights are used on the link scale (logit or log), which is
      # unbounded and usually better suited to regression.
      posterior_samples <- extract_metric(model, metric_string, num_draws, transform = FALSE)[, dyad_ids, drop = FALSE]
    } else if (target_name == "node") {
      node_ids <- lookup_node_ids(model, data_list[[i]], param_names[1], i)
      posterior_samples <- extract_metric(model, metric_string, num_draws)[, node_ids, drop = FALSE]
    } else if (target_name == "global") {
      posterior_samples <- extract_metric(model, metric_string, num_draws)
    } else {
      stop("Network feature \"", target_name, "\" not supported. Use edge, node or global.", call. = FALSE)
    }

    # Make sure posterior samples are in matrix format (for num_draws=1).
    posterior_samples_list[[i]] <- matrix(posterior_samples, nrow = num_draws)
  }

  # Generate list of dataframes.
  imputed_dataframes <- lapply(1:num_draws, function(i) {
    # For each posterior draw, combine dataframes from different edge models and data
    new_data_list <- lapply(1:length(edgemodel_list), function(j) {
      new_data <- data_list[[j]]
      if (z_score) {
        new_data[new_bison_term] <- (posterior_samples_list[[j]][i, ] - mean(posterior_samples_list[[j]][i, ]))/sd(posterior_samples_list[[j]][i, ])
      } else {
        new_data[new_bison_term] <- posterior_samples_list[[j]][i, ]
      }
      new_data["bison_network"] <- j
      new_data
    })
    new_data <- dplyr::bind_rows(new_data_list)
    new_data$bison_network <- factor(new_data$bison_network)
    new_data
  })

  imputed_dataframes <- dplyr::bind_rows(imputed_dataframes, .id=".imp")
  imputed_dataframes <- dplyr::mutate(imputed_dataframes, .imp=as.integer(.imp))

  original_dataframes <- lapply(1:length(edgemodel_list), function(j) {
    new_data <- data_list[[j]]
    new_data["bison_network"] <- j
    new_data$bison_network <- factor(new_data$bison_network)
    new_data
  })
  original_dataframes <- dplyr::bind_rows(original_dataframes)
  original_dataframes[".imp"] <- 0

  combined_dataframes <- dplyr::bind_rows(list(original_dataframes, imputed_dataframes))

  # mice warns about "logged events" when imputed columns are constant within an
  # imputation, which is expected here, so that warning is muffled.
  withCallingHandlers(
    mice::as.mids(as.data.frame(combined_dataframes)),
    warning = function(w) {
      if (grepl("logged events", conditionMessage(w))) invokeRestart("muffleWarning")
    }
  )
}

# Find the indices (in the edge model) of the nodes named in a column of a
# regression data frame, with a clear error if any are missing.
lookup_node_ids <- function(model, data, column, model_number) {
  if (!column %in% colnames(data)) {
    stop("Column `", column, "` not found in data frame ", model_number, ".", call. = FALSE)
  }
  node_names <- as.character(data[[column]])
  ids <- match(node_names, names(model$node_to_idx))
  if (anyNA(ids)) {
    missing <- unique(node_names[is.na(ids)])
    stop(
      "Node(s) ", paste(utils::head(missing, 5), collapse = ", "),
      if (length(missing) > 5) ", ..." else "",
      " in column `", column, "` of data frame ", model_number,
      " aren't in edge model ", model_number, ".",
      call. = FALSE
    )
  }
  ids
}
