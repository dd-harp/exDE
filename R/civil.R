 # generic methods for exogenous forcing by civil factors

 #' @title Set up exogenous variables for civil forcing
 #' @description This method dispatches on the type of `pars$CIVIL`.
 #' @param t current simulation time
 #' @param pars a [list]
 #' @return [list]
 #' @export
 Civil <- function(t, pars) {
   UseMethod("Civil", pars$CIVIL)
 }

 #' @title Set up exogenous variables for civil forcing
 #' @description Implements [Civil] for the null model of exogenous forcing (do nothing)
 #' @inheritParams Civil
 #' @return [list]
 #' @export
 Civil.null <- function(t, pars) {pars}

 #' @title Set up the null model for exogenous forcing (do nothing)
 #' @param pars a [list]
 #' @return [list]
 #' @export
 setup_civil_null <- function(pars) {
   CIVIL <- list()
   class(CIVIL) <- 'null'
   pars$CIVIL <- CIVIL
   return(pars)
 }

 #' @title Set up exogenous variables for civil forcing
 #' @description Implements [Civil] for civil forcing
 #' @param t current simulation time
 #' @param pars a [list]
 #' @return pars a [list]
 #' @export
 Civil.forced <- function(t, pars) {
   pars = Weather(t, pars)
   pars = Hydrology(t, pars)
   return(pars)
 }

 #' @title Make parameters for the null model of civil forcing (do nothing)
 #' @param pars a [list]
 #' @return [list]
 #' @export
 setup_civil_forced <- function(pars) {
   CIVIL <- list()
   class(CIVIL) <- 'forced'
   pars$CIVIL <- CIVIL
   pars = setup_shock_null(pars)
   pars = setup_development_null(pars)
   return(pars)
 }

 #' @title Check civil
 #' @param pars a [list]
 #' @return [list]
 #' @export
 check_civil <- function(pars) {
   UseMethod("check_civil", pars$EfSz)
 }

 #' @title Check civil
 #' @param pars a [list]
 #' @return [list]
 #' @export
 check_civil.null <- function(pars) {
   setup_civil_forced(pars)
 }

 #' @title Check civil
 #' @param pars a [list]
 #' @return [list]
 #' @export
 check_civil.forced<- function(pars) {pars}

