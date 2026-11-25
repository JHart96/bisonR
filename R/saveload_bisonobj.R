#' Save a fitted bisonR model to disk
#'
#' Fitted models contain a CmdStan fit object whose draws live in temporary CSV
#' files, so `saveRDS()` on its own does not produce a file that can be reloaded
#' in a new R session. `bison_save()` reads all draws into memory first and then
#' saves the complete object as a single `.rds` file.
#'
#' @param bisonobj A fitted model, such as the output of [bison_model()].
#' @param file Path of the `.rds` file to write. `.rds` is appended if missing.
#'
#' @return The path to the saved file, invisibly.
#' @seealso [bison_load()]
#'
#' @examples
#' \donttest{
#' sim_data <- simulate_bison_model("binary", aggregated = TRUE)
#' fit_edge <- bison_model(
#'   (event | duration) ~ dyad(node_1_id, node_2_id),
#'   data = sim_data$df_sim,
#'   model_type = "binary_conjugate",
#'   priors = get_default_priors("binary_conjugate")
#' )
#' path <- tempfile(fileext = ".rds")
#' bison_save(fit_edge, path)
#' fit_reloaded <- bison_load(path)
#' }
#' @export
bison_save <- function(bisonobj, file) {
  if (!grepl("\\.rds$", file, ignore.case = TRUE)) {
    file <- paste0(file, ".rds")
  }
  obj <- bisonobj
  if (!is.null(obj$fit) && inherits(obj$fit, "CmdStanFit")) {
    # Load draws, sampler diagnostics etc. into the R object so it is self-contained.
    obj$fit$draws()
    try(obj$fit$sampler_diagnostics(), silent = TRUE)
    try(obj$fit$init(), silent = TRUE)
    try(obj$fit$profiles(), silent = TRUE)
  }
  # The compiled model object only points at a temporary executable, so drop it.
  obj$stan_model <- NULL
  obj$bisonR_version <- as.character(utils::packageVersion("bisonR"))
  saveRDS(obj, file)
  invisible(file)
}

#' Load a bisonR model saved with `bison_save()`
#'
#' @param file Path of an `.rds` file written by [bison_save()].
#'
#' @return The fitted model object.
#' @seealso [bison_save()]
#' @inherit bison_save examples
#' @export
bison_load <- function(file) {
  if (!file.exists(file) && file.exists(paste0(file, ".rds"))) {
    file <- paste0(file, ".rds")
  }
  readRDS(file)
}
