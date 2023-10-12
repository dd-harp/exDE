# specialized methods for the a eir-based travel_malaria model

#' @title A model for travel malaria based on the eir in destinations
#' @description Implements [travel_malaria] through a model for an eir-based model
#' @inheritParams travel_malaria
#' @return a [numeric]
#' @export
travel_malaria.eir <- function(t, pars) {
  with(pars$TRAVEL,{
    eir = travel_eir_scale*travel_eir_t(t, pars)*frac_time_spent_traveling
    return(F_foi(pars$MYZpar$b*eir, pars))
})}

#' @title Set up parameters and function for an eir-based travel_malaria model
#' @param pars a [list]
#' @param travel_Opts a list, overwrites default values
#' @param frac_time_spent_traveling is the fraction of time spent travel
#' @param travel_eir_scale a non-negative numeric value to scale the mean eir experienced while traveling
#' @param travel_eir_t the temporal pattern for eir experienced while traveling
#' @return none
#' @export
setup_travel_eir <- function(pars, travel_Opts = list(), frac_time_spent_traveling = 0.01, travel_eir_scale=0, travel_eir_t=NULL) {with(travel_Opts,{

  TRAVEL <- list()
  class(TRAVEL) <- 'eir'
  pars$TRAVEL <- TRAVEL

  pars$TRAVEL$travel_eir_scale = travel_eir_scale
  if(is.null(travel_eir_t)) travel_eir_t = function(t, pars){1}
  pars$TRAVEL$travel_eir_t = travel_eir_t
  pars$TRAVEL$frac_time_spent_traveling = frac_time_spent_traveling

  return(pars)
})}
