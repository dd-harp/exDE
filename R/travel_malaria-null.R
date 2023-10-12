# specialized methods for a null model of travel malaria

#' @title travel, a null model
#' @description Implements [travel_malaria] for the null model (do nothing)
#' @inheritParams travel_malaria
#' @return a [numeric]
#' @export
travel_malaria.null <- function(t, pars) {
  return(0)
}

#' @title A function to set up malaria importation
#' @description Setup a null model for travel malaria
#' @param pars a [list]
#' @return a [list]
#' @export
setup_travel_null = function(pars){

  TRAVEL <- list()
  class(TRAVEL) <- 'null'
  pars$TRAVEL <- TRAVEL

  return(pars)
}
