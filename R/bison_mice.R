#' Prepare network data for multiple imputation with the mice package
#'
#' Combines a regression data frame with posterior draws of a network metric,
#' producing one completed copy of the data for each draw. This is the step
#' [bison_brm()] uses to carry the uncertainty in the network into a
#' regression. Use `bison_mice()` directly to fit the regression with software
#' other than brms, or to inspect the data that [bison_brm()] would use.
#'
#' @details
#' **Multiple imputation** is a method for missing data: the analysis is
#' repeated on several plausible versions of the data and the results are
#' combined, so the final estimates include the uncertainty about the missing
#' values. bisonR treats a network metric as a variable whose true values are
#' unknown but described by the posterior distribution of the edge model (the
#' distribution of values consistent with the data and priors). Each posterior
#' draw, a single sample from that distribution, gives one version of the data.
#'
#' The [mice](https://cran.r-project.org/package=mice) package is an R package
#' for multiple imputation. It stores the versions of a dataset in a
#' `mids` object ("multiply imputed data set"). Regression packages that
#' support mice, such as brms ([brms::brm_multiple()]), and functions in mice
#' itself ([with()][mice::with.mids()] and [mice::pool()]) can fit a model to each
#' version and combine the results. Use [mice::complete()] to get the data
#' frames back.
#'
#' The network metric is added as a column named after the metric, such as
#' `bison_node_strength` or `bison_edge_weight`. A factor `bison_network` gives
#' the position of the edge model in `edgemodel_list` for each row. Edge weights
#' are on the link scale (log-odds for binary models, log rates for count
#' models), which is unbounded and usually better suited to regression. Node
#' and global metrics are calculated from edge weights on the natural scale, as
#' in [extract_metric()].
#'
#' @param edgemodel_list A fitted `bison_model` object, or a list of them for
#'   an analysis of several networks.
#' @param data_list A data frame of regression variables, or a list of data
#'   frames with one per edge model, in the same order. For edge and node
#'   metrics, each row refers to a dyad or node in the corresponding edge model.
#'   For global metrics, give one data frame with one row per edge model.
#' @param num_draws Number of posterior draws to use, which is the number of
#'   imputed datasets. 100 is a reasonable value for a final analysis; small
#'   values (such as 5) are useful for testing code.
#' @param param_names Names of the columns in the data that identify the
#'   targets of the metric: two columns (e.g. `c("node_1_id", "node_2_id")`) for
#'   edge weights and one (e.g. `"node_id"`) for node metrics. The values must
#'   match the node names in the edge model. Not used for global metrics.
#' @param target_name The level of the metric: `"edge"`, `"node"` or
#'   `"global"`.
#' @param metric_name Name of the metric, without the target prefix, e.g.
#'   `"weight"` (the only edge metric), `"strength"`, `"degree[0.2]"` or `"cv"`.
#'   See [extract_metric()] for the metrics available.
#' @param z_score If `TRUE`, standardise the network variable within each
#'   posterior draw. Equivalent to `scale = "z"`.
#' @param scale How to transform the network variable within each posterior
#'   draw: `"none"` (default), `"centre"` (subtract the mean) or `"z"`
#'   (subtract the mean and divide by the standard deviation). For edge and node
#'   metrics this is done across the rows of each data frame separately; global
#'   metrics are transformed across networks.
#'
#' @return A `mids` object from the mice package, holding the original data
#'   (with the network variable missing) and one completed dataset per
#'   posterior draw.
#'
#' @seealso [bison_brm()], which calls this function and fits the regression;
#'   [extract_metric()] for posterior draws of a metric on its own.
#' @examplesIf requireNamespace("brms", quietly = TRUE) && requireNamespace("mice", quietly = TRUE) && !identical(Sys.getenv("BISONR_TEST_BRMS"), "false")
#' sim_data <- simulate_bison_model("binary", aggregated = TRUE)
#' df <- sim_data$df_sim
#' fit_edge <- bison_model(
#'   (event | effort) ~ dyad(node_1_id, node_2_id),
#'   data = df,
#'   model_type = "binary_conjugate",
#'   priors = get_default_priors("binary_conjugate")
#' )
#' # One imputed dataset per posterior draw of node strength
#' df_nodal <- data.frame(node = levels(df$node_1_id))
#' imputed <- bison_mice(fit_edge, df_nodal, "node", "node", "strength", num_draws = 5)
#' head(mice::complete(imputed, 1))
#' @export
bison_mice <- function(edgemodel_list, data_list, param_names, target_name, metric_name, num_draws=100, z_score=FALSE,
                       scale = c("none", "centre", "z")) {
  scale <- match.arg(scale)
  if (isTRUE(z_score)) {
    if (!scale %in% c("none", "z")) {
      warning("`z_score = TRUE` overrides `scale = \"", scale, "\"`; using scale = \"z\".", call. = FALSE)
    }
    scale <- "z"
  }
  check_suggested("mice", "bison_mice")
  if (target_name == "global" && is.data.frame(data_list)) {
    data_list <- lapply(seq_len(nrow(data_list)), function(i) data_list[i, , drop = FALSE])
  }

  # If the data aren't lists yet, convert them
  if (inherits(edgemodel_list, "bison_model")) edgemodel_list <- list(edgemodel_list)
  if (is.data.frame(data_list)) data_list <- list(data_list)

  new_bison_term <- bison_term_name(target_name, metric_name)

  posterior_samples_list <- list()

  if (length(edgemodel_list) != length(data_list)) {
    if (target_name == "global") {
      stop(
        "Got ", length(edgemodel_list), " edge model(s) but ", length(data_list),
        " row(s) of data. Global metrics need one row of data per edge model, in the same order.",
        call. = FALSE
      )
    }
    stop(
      "Got ", length(edgemodel_list), " edge model(s) but ", length(data_list),
      " data frame(s). Provide one data frame per edge model.",
      call. = FALSE
    )
  }

  # bison_network is set from the position of each edge model. A column with
  # this name in the data must agree with that order.
  for (j in seq_along(data_list)) {
    if ("bison_network" %in% colnames(data_list[[j]]) &&
        any(as_node_name(data_list[[j]]$bison_network) != as.character(j))) {
      stop(
        "The `bison_network` column in data frame ", j, " doesn't match the position of its edge ",
        "model (", j, "). `bison_network` is set from the order of the edge models, so put the ",
        "data in the same order or remove the column.",
        call. = FALSE
      )
    }
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

    if (any(!is.finite(posterior_samples))) {
      stop(
        "The network variable ", metric_string, " isn't finite for some rows of data frame ", i,
        " (for example, closeness of a node with no dyads). Remove those rows or use another metric.",
        call. = FALSE
      )
    }

    # Make sure posterior samples are in matrix format (for num_draws=1).
    posterior_samples_list[[i]] <- matrix(posterior_samples, nrow = num_draws)
  }

  # Global metrics have one value per network, so scaling is done across the
  # networks within each draw rather than within each network.
  if (target_name == "global" && scale != "none") {
    if (length(posterior_samples_list) < 2) {
      stop("Scaling a global metric needs at least two networks.", call. = FALSE)
    }
    combined <- do.call(cbind, posterior_samples_list)
    combined <- t(apply(combined, 1, scale_draw, scale = scale))
    combined <- matrix(combined, nrow = num_draws)
    posterior_samples_list <- lapply(seq_along(posterior_samples_list), function(j) combined[, j, drop = FALSE])
    scale <- "none"
  }

  # Generate list of dataframes.
  imputed_dataframes <- lapply(1:num_draws, function(i) {
    # For each posterior draw, combine dataframes from different edge models and data
    new_data_list <- lapply(1:length(edgemodel_list), function(j) {
      new_data <- data_list[[j]]
      new_data[new_bison_term] <- scale_draw(posterior_samples_list[[j]][i, ], scale)
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
  node_names <- as_node_name(data[[column]])
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

# Transform one posterior draw of a network variable.
scale_draw <- function(x, scale) {
  if (scale == "centre") {
    return(x - mean(x))
  }
  if (scale == "z") {
    if (length(x) < 2 || stats::sd(x) == 0) {
      stop("Can't standardise a network variable with fewer than two distinct values.", call. = FALSE)
    }
    return((x - mean(x)) / stats::sd(x))
  }
  x
}
