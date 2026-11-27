#' Prior predictive check
#'
#' Fits the model using only the priors (ignoring the observed events) and plots
#' the data the model would predict. This shows whether the priors imply
#' plausible data before any data are seen. A prior describes which values of
#' a parameter are plausible before seeing the data. Use it before fitting
#' with [bison_model()]: if the predicted rates are implausible for the
#' species and sampling method, change the priors.
#'
#' @param formula,data,model_type,priors As for [bison_model()]. The data are
#'   used for the structure of the model (dyads and observation effort), not
#'   for the events.
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
#' @seealso [prior_check()] to plot the priors themselves,
#'   [get_default_priors()], [plot_predictions()]
#'
#' @examples
#' sim_data <- simulate_bison_model("binary", aggregated = TRUE)
#' prior_predictive_check(
#'   (event | effort) ~ dyad(node_1_id, node_2_id),
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
