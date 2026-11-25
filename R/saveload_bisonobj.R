#' Function for Saving Bison Models as RDS
#' 
#' Saves Bison Objects. Saves as two RDS objects- a cmdstanr object and a list 
#' of other information.
#' 
#' Objects will save with the format: "filename_cmdstanfit.RDS" and "filename.RDS".
#' 
#' @examples
#' 
#' library(bisonR)
#' sim_data <- simulate_bison_model("binary", aggregated = TRUE)
#' df <- sim_data$df_sim
#' priors <- get_default_priors("binary")
#'
#' fit_edge <- bison_model(
#' (event | duration) ~ dyad(node_1_id, node_2_id), 
#' data=df, 
#' model_type="binary",
#' priors=priors
#' )
#' 
#' bison_save(fit_edge, filename = "testsave")
#' 
#' @seealso [bison_load()]

bison_save = function(bisonobj, filename){
  
  fitname = paste(filename, "_cmdstanfit.RDS", sep ="")
  bisonobj$fit$save_object(fitname)
  
  therest = bisonobj
  therest$fit = NULL
  therestname = paste(filename, ".RDS", sep = "")
  saveRDS(therest, therestname)
  
}

#' Function for Loading Saved Bison Models as RDS
#' 
#' Load Bison Objects saved by the \code{\link{bison_save}}. Saves as two RDS objects- a cmdstanr object and a list 
#' of other information.
#' 
#' Objects will save with the format: "filename_cmdstanfit.RDS" and "filename.RDS".
#' 
#' #' @examples
#' 
#' library(bisonR)
#' sim_data <- simulate_bison_model("binary", aggregated = TRUE)
#' df <- sim_data$df_sim
#' priors <- get_default_priors("binary")
#'
#' fit_edge <- bison_model(
#' (event | duration) ~ dyad(node_1_id, node_2_id), 
#' data=df, 
#' model_type="binary",
#' priors=priors
#' )
#' 
#' bison_save(fit_edge, filename = "testsave")
#' rm(fit_edge)
#' fit_edge = bison_load(filename = "testsave)
#' 
#' @seealso [bison_save()]


bison_load = function(filename){
  
  fitname = paste(filename, "_cmdstanfit.RDS", sep ="")
  bisonfit = readRDS(fitname)
  
  therestname = paste(filename, ".RDS", sep = "")
  therest = readRDS(therestname)
  
  bisonobj = therest
  bisonobj$fit = bisonfit
  return(bisonobj)
}
# loaded.bison.mod = bison_load("testsave")


