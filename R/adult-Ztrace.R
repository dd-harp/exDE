# specialized methods for the adult mosquito Ztrace model

#' @title Number of infective adults in each patch
#' @description Implements [F_Z] for the Ztrace  model.
#' @inheritParams F_Z
#' @return a [numeric] vector of length `nHabitats`
#' @export
F_Z.Ztrace <- function(t, y, pars) {
  with(pars$MYZpar, Zm*Zf(t, pars))
}

#' @title Number of infective adults in each patch
#' @description Implements [F_Z] for the Ztrace  model.
#' @inheritParams F_Z_lag
#' @return a [numeric] vector of length `nHabitats`
#' @export
F_Z_lag.Ztrace <- function(t, y, pars, lag) {
  pars$MYZpar$Zf(t-lag)
}

#' @title Derivatives for aquatic stage mosquitoes
#' @description Implements [dMYZdt] for the Ztrace (forced emergence) model.
#' @inheritParams dMYZdt
#' @return a [numeric] vector
#' @export
dMYZdt.Ztrace <- function(t, y, pars, Lambda, kappa, MosyBehavior){
  numeric(0)
}

#' @title Add indices for aquatic stage mosquitoes to parameter list
#' @description Implements [make_indices_MYZ] for Ztrace (forced emergence) model.
#' @inheritParams make_indices_MYZ
#' @return none
#' @export
make_indices_MYZ.Ztrace <- function(pars) {
  pars$MYZ_ix <- integer(0)
  return(pars)
}

#' @title Make parameters for Ztrace aquatic mosquito model
#' @param pars an [environment]
#' @param Zm a vector of mean mosquito densities
#' @param Zf a [function] of the form Zf(t, pars) that computes temporal fluctuations
#' @return none
#' @export
make_parameters_MYZ_Ztrace <- function(pars, Zm=1, Zf=NULL) {
  stopifnot(is.numeric(Zm))
  MYZpar <- list()
  class(MYZpar) <- 'Ztrace'
  if(is.null(Zf)) {
    MYZpar$Zf = function(t, pars){1}
  }
  pars$MYZpar <- MYZpar
  return(pars)
}

#' @title Make parameters for Ztrace aquatic mosquito model
#' @param pars an [environment]
#' @param MYZ0 is set to NULL for the Ztrace model
#' @return none
#' @export
make_inits_MYZ_Ztrace<- function(pars, MYZ0=NULL) {
  pars$MYZinits = numeric(0)
  return(pars)
}

#' @title Return initial values as a vector
#' @description Implements [get_inits_MYZ] for the GeRM model.
#' @inheritParams get_inits_MYZ
#' @return none
#' @export
get_inits_MYZ.Ztrace <- function(pars){
  numeric(0)
}

