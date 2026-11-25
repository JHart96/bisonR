# Transform edge samples from the link scale to the natural scale
# (probabilities for binary models, rates for count models).
edge_transform <- function(obj, edge_samples = obj$edge_samples) {
  if (obj$model_type == "binary") {
    return(stats::plogis(edge_samples))
  }
  if (obj$model_type == "count") {
    return(exp(edge_samples))
  }
  stop("Unknown model type \"", obj$model_type, "\".", call. = FALSE)
}

# Build an igraph network with every node in the model (in node_to_idx order)
# and one edge per dyad (in dyad order). Edge weights are added later.
bison_base_igraph <- function(obj, dyad_rows = seq_len(obj$num_dyads)) {
  node_names <- names(obj$node_to_idx)
  edges <- data.frame(
    from = node_names[obj$dyad_to_idx[dyad_rows, 1]],
    to = node_names[obj$dyad_to_idx[dyad_rows, 2]],
    stringsAsFactors = FALSE
  )
  igraph::graph_from_data_frame(
    edges,
    directed = obj$directed,
    vertices = data.frame(name = node_names, stringsAsFactors = FALSE)
  )
}

#' Convert a fitted edge model into a list of igraph networks
#'
#' Each network is one draw from the posterior distribution of networks, with
#' edge weights on the natural scale (probabilities or rates). The vertices are
#' always all the nodes in the model, in the same order as the model's node
#' levels, so results can be matched back to node names.
#'
#' @param obj A fitted `bison_model` object.
#' @param num_draws Number of networks to draw from the posterior.
#'
#' @return A list of `igraph` objects.
#'
#' @examples
#' sim_data <- simulate_bison_model("binary", aggregated = TRUE)
#' fit_edge <- bison_model(
#'   (event | effort) ~ dyad(node_1_id, node_2_id),
#'   data = sim_data$df_sim,
#'   model_type = "binary_conjugate",
#'   priors = get_default_priors("binary_conjugate")
#' )
#' nets <- bison_to_igraph(fit_edge, num_draws = 10)
#' # Posterior distribution of the strength of each node
#' strengths <- sapply(nets, igraph::strength)
#' @export
bison_to_igraph <- function(obj, num_draws = 1000) {
  if (is.null(obj$dyad_to_idx)) {
    stop("This model has no edge weights (it has no dyad() term).", call. = FALSE)
  }
  base_net <- bison_base_igraph(obj)
  edge_samples <- edge_transform(obj)
  draw_ids <- sample_draw_ids(nrow(edge_samples), num_draws)
  lapply(draw_ids, function(draw_id) {
    net <- base_net
    igraph::E(net)$weight <- edge_samples[draw_id, ]
    net
  })
}

#' Convert duration data to binary
#'
#' @description Uses a naive method for converting duration data to binary data:
#' durations are divided by `max_event_length` and rounded. Use with caution, and
#' make sure all durations are in the same units.
#'
#' @param df Original data frame, ready for use with [bison_model()].
#' @param event_name Name of the column with event durations.
#' @param divisor_name Name of the column with observation durations.
#' @param max_event_length The maximum length of time a single social event could
#'   last.
#'
#' @return The input data frame with the `event_name` and `divisor_name` columns
#'   replaced by counts suitable for a binary model.
#'
#' @examples
#' df <- data.frame(
#'   node_1 = c("A", "A", "B"), node_2 = c("B", "C", "C"),
#'   time_together = c(30, 0, 10), time_observed = c(60, 60, 60)
#' )
#' convert_duration_to_binary(df, "time_together", "time_observed", max_event_length = 10)
#' @export
convert_duration_to_binary <- function(df, event_name, divisor_name, max_event_length) {
  new_df <- df
  new_df[[event_name]] <- round(df[[event_name]] / max_event_length)
  new_df[[divisor_name]] <- round(df[[divisor_name]] / max_event_length)
  new_df
}

#' Convert a group-by-individual matrix to bisonR format
#'
#' Group-by-individual (GBI) data record which individuals were seen together in
#' each group or sampling event. This function produces one row for each group
#' and each pair of individuals where *at least one* of the pair was present.
#' The `event` column is 1 if both were present and 0 if only one was. Groups
#' where neither individual was seen carry no information about the pair and are
#' left out.
#'
#' @param gbi A matrix with R rows (groups) and C columns (individuals) of 0s and
#'   1s. If it has column names they are used as node names; otherwise nodes are
#'   numbered `1:C`.
#' @param group_properties Optional vector of length R with a property of each
#'   group (e.g. location), added as column `group_property`.
#' @param individual_properties Optional vector of length C with a property of
#'   each individual, added as columns `node_1_property` and `node_2_property`.
#' @param individual_constraints Optional vector of length C. Only pairs of
#'   individuals with the same value are included (for example, individuals that
#'   could have been seen together because they were in the study area at the
#'   same time).
#'
#' @return A data frame with columns `node_1`, `node_2`, `event` and `group_id`,
#'   plus any property columns, ready for use with a binary [bison_model()] with
#'   an effort of 1, e.g. `(event | 1) ~ dyad(node_1, node_2)`.
#'
#' @examples
#' gbi <- matrix(rbinom(20 * 5, 1, 0.4), 20, 5, dimnames = list(NULL, LETTERS[1:5]))
#' df <- convert_gbi_to_bison(gbi)
#' head(df)
#' @export
convert_gbi_to_bison <- function(gbi, group_properties = NULL, individual_properties = NULL,
                                 individual_constraints = NULL) {
  gbi <- as.matrix(gbi)
  num_individuals <- ncol(gbi)
  if (is.null(individual_constraints)) {
    individual_constraints <- rep(0, num_individuals)
  }

  node_names <- colnames(gbi)
  if (is.null(node_names)) {
    node_names <- seq_len(num_individuals)
  }

  # All pairs i < j that satisfy the constraints, ordered by i then j.
  pairs <- which(
    upper.tri(matrix(0, num_individuals, num_individuals)) &
      outer(individual_constraints, individual_constraints, "=="),
    arr.ind = TRUE
  )
  pairs <- pairs[order(pairs[, 1], pairs[, 2]), , drop = FALSE]

  rows <- lapply(seq_len(nrow(pairs)), function(k) {
    i <- pairs[k, 1]
    j <- pairs[k, 2]
    present <- which(gbi[, i] == 1 | gbi[, j] == 1)
    if (length(present) == 0) {
      return(NULL)
    }
    df_new <- data.frame(
      node_1 = node_names[i],
      node_2 = node_names[j],
      event = gbi[present, i] * gbi[present, j],
      group_id = present
    )
    if (!is.null(group_properties)) {
      df_new$group_property <- group_properties[present]
    }
    if (!is.null(individual_properties)) {
      df_new$node_1_property <- individual_properties[i]
      df_new$node_2_property <- individual_properties[j]
    }
    df_new
  })
  df <- do.call(rbind, rows)
  if (is.null(df)) {
    df <- data.frame(node_1 = numeric(), node_2 = numeric(), event = numeric(), group_id = numeric())
  }
  rownames(df) <- NULL
  df
}

# Directory where compiled Stan models are cached between R sessions. Set
# options(bisonR.cache_dir = ...) to change it.
bison_cache_dir <- function() {
  cache_dir <- getOption("bisonR.cache_dir", tools::R_user_dir("bisonR", which = "cache"))
  ok <- dir.exists(cache_dir) || dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  if (!ok || file.access(cache_dir, 2) != 0) {
    cache_dir <- file.path(tempdir(), "bisonR")
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }
  cache_dir
}

# Compile (or load from the cache) one of the package's Stan models. Compiled
# models are keyed on the model source and the CmdStan version, so they are only
# recompiled when either changes.
build_stan_model <- function(model_name) {
  cmdstan_ok <- !is.null(tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL))
  if (!cmdstan_ok) {
    stop(
      "CmdStan isn't installed, so this model can't be fitted. ",
      "Install it with cmdstanr::install_cmdstan(), or use a conjugate model ",
      "(model_type = \"binary_conjugate\" or \"count_conjugate\").",
      call. = FALSE
    )
  }
  source_file <- system.file("stan", paste0(model_name, ".stan"), package = "bisonR")
  key <- paste0(
    substr(unname(tools::md5sum(source_file)), 1, 12), "_",
    gsub("[^0-9]", "", as.character(cmdstanr::cmdstan_version()))
  )
  stan_file <- file.path(bison_cache_dir(), paste0(model_name, "_", key, ".stan"))
  if (!file.exists(stan_file)) {
    file.copy(source_file, stan_file)
  }
  cmdstanr::cmdstan_model(stan_file, stanc_options = list("O1"))
}

bison_colors <- c("#2b9392", "#5aabaa", "#89c3c2", "#b8dbda")
bison_highlight <- "#d1495b"

col2rgba <- function(col, alpha) {
  rgba <- grDevices::col2rgb(col, alpha = TRUE) / 255
  grDevices::rgb(rgba[1], rgba[2], rgba[3], alpha)
}

# Stop with an informative message if an optional (Suggests) package is missing.
check_suggested <- function(pkg, fn) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Package '", pkg, "' is needed for ", fn, "() but is not installed. ",
      "Install it with install.packages(\"", pkg, "\").",
      call. = FALSE
    )
  }
  invisible(TRUE)
}
