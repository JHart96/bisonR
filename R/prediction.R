#' Prior predictive check
#'
#' Fits the model using only the priors (ignoring the observed events) and plots
#' the data the model would predict. This shows whether the priors imply
#' plausible data before any data are seen.
#'
#' @param formula,data,model_type,priors As for [bison_model()].
#' @param plot_type Type of predictive plot, passed to [plot_predictions()] as
#'   `type`: `"density"` or `"stats"`.
#' @param options A list of extra options. Currently only `directed`
#'   (`TRUE`/`FALSE`) is used; prefer the `directed` argument.
#' @param directed Whether the network is directed.
#' @param num_draws Number of datasets to simulate from the priors.
#' @param ... Further arguments passed to [bison_model()].
#'
#' @return The prior-only model, invisibly.
#'
#' @examples
#' sim_data <- simulate_bison_model("binary", aggregated = TRUE)
#' prior_predictive_check(
#'   (event | duration) ~ dyad(node_1_id, node_2_id),
#'   data = sim_data$df_sim,
#'   model_type = "binary_conjugate",
#'   priors = get_default_priors("binary_conjugate")
#' )
#' @export
prior_predictive_check <- function(formula, data, model_type, plot_type = "density",
                                   options = NULL, priors = NULL,
                                   directed = isTRUE(options$directed), num_draws = 50, ...) {
  fit <- bison_model(
    formula = formula,
    data = data,
    model_type = model_type,
    directed = directed,
    priors = priors,
    priors_only = TRUE,
    ...
  )
  plot_predictions(fit, num_draws = num_draws, draw_data = FALSE, type = plot_type)
  invisible(fit)
}
