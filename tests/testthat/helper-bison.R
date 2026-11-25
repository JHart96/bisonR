# Shared helpers for the bisonR test suite.

skip_if_no_cmdstan <- function() {
  testthat::skip_if(
    is.null(tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL)),
    "CmdStan is not installed"
  )
}

# Silence MCMC progress output in tests.
quiet <- function(expr) {
  suppressMessages(utils::capture.output(result <- expr))
  result
}

conjugate_priors <- function(model_type) {
  get_default_priors(paste0(model_type, "_conjugate"))
}

# A minimal bison_model object with fixed (known) edge weights, for testing
# metric calculations without fitting a model.
fake_bison_model <- function(nodes, dyads, weights, directed = FALSE, num_draws = 5) {
  node_to_idx <- seq_along(nodes)
  names(node_to_idx) <- nodes
  dyad_to_idx <- cbind(match(dyads[, 1], nodes), match(dyads[, 2], nodes))
  edge_samples <- matrix(stats::qlogis(weights), num_draws, length(weights), byrow = TRUE)
  obj <- list(
    chain = edge_samples, edge_samples = edge_samples,
    num_nodes = length(nodes), num_dyads = nrow(dyad_to_idx),
    node_to_idx = node_to_idx, dyad_to_idx = dyad_to_idx,
    dyad_names = paste(dyads[, 1], dyads[, 2], sep = if (directed) " -> " else " <-> "),
    directed = directed, model_type = "binary", conjugate = TRUE
  )
  class(obj) <- "bison_model"
  obj
}

# brms models are slow to compile, so CI runs the brms tests on one platform
# only (BISONR_TEST_BRMS = "false" elsewhere).
skip_if_brms_disabled <- function() {
  testthat::skip_if(identical(Sys.getenv("BISONR_TEST_BRMS"), "false"), "brms tests disabled on this platform")
}
