# specialized methods for a static model of travel malaria

#' @title travel, a static model
#' @description Implements [travel_malaria] for the static model (do nothing)
#' @inheritParams travel_malaria
#' @return a [numeric]
#' @export
travel_malaria.static <- function(t, pars) {
  return(pars$vars$delta)
}

#' @title A function to set up malaria importation
#' @description Setup a static model for travel malaria
#' @param pars a [list]
#' @param delta the travel FoI
#' @return a [list]
#' @export
setup_travel_static = function(pars, delta=0){
  TRAVEL <- list()
  class(TRAVEL) <- 'static'
  pars$TRAVEL <- TRAVEL
  pars$vars$delta = delta
  return(pars)
}
