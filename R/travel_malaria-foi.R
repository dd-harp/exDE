# specialized methods for the a FoI-based travel_malaria model

#' @title A model for the travel FoI
#' @description Implements [travel_malaria] through a model for the travel FoI
#' @inheritParams travel_malaria
#' @return a [numeric]
#' @export
travel_malaria.foi <- function(t, pars) {
  with(pars$TRVpar,{
    return(delta_scale*delta_t(t, pars))
})}

#' @title Set up parameters and function for the FoI-based travel_malaria model
#' @param pars a [list]
#' @param travel_Opts a list, overwrites default values
#' @param delta_scale a non-negative numeric value to scale the mean travel_malaria FoI
#' @param delta_t the temporal pattern for travel_malaria FoI
#' @return none
#' @export
setup_travel_foi <- function(pars, travel_Opts = list(), delta_scale=0, delta_t=NULL) {with(travel_Opts,{

  TRVpar <- list()
  class(TRVpar) <- 'foi'
  pars$TRVpar <- TRVpar

  pars$TRVpar$delta_scale = delta_scale
  if(is.null(delta_t)) delta_t = function(t, pars){1}
  pars$TRVpar$delta_t = delta_t

  return(pars)
})}
