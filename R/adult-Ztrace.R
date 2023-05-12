# specialized methods for the adult mosquito Ztrace model

#' @title Compute bloodfeeding and mortality rates
#' @description Implements [MosquitoBehavior] for the Ztrace model.
#' @inheritParams MosquitoBehavior
#' @return a named [list]
#' @export
MosquitoBehavior.Ztrace <- function(t, y, pars) {
  pars$MYZpar$f <- pars$MYZpar$f0
  pars$MYZpar$q <- pars$MYZpar$q0
  return(pars)
}

#' @title Number of infective adults in each patch
#' @description Implements [F_Z] for the Ztrace  model.
#' @inheritParams F_Z
#' @return a [numeric] vector of length `nHabitats`
#' @export
F_Z.Ztrace <- function(t, y, pars) {
  with(pars$MYZpar, return(Zm*Zf(t, pars)))
}

#' @title Number of eggs laid by adult mosquitoes
#' @description Implements [F_eggs] for the Ztrace model.
#' @inheritParams F_eggs
#' @return a [numeric] vector of length `nPatches`
#' @export
F_eggs.Ztrace <- function(t, y, pars) {
  return(numeric(0))
}

#' @title Derivatives for aquatic stage mosquitoes
#' @description Implements [dMYZdt] for the Ztrace (forced emergence) model.
#' @inheritParams dMYZdt
#' @return a [numeric] vector
#' @export
dMYZdt.Ztrace <- function(t, y, pars, Lambda, kappa){
  numeric(0)
}

#' @title Add indices for aquatic stage mosquitoes to parameter list
#' @description Implements [make_indices_MYZ] for Ztrace (forced emergence) model.
#' @inheritParams make_indices_MYZ
#' @return none
#' @export
make_indices_MYZ.Ztrace <- function(pars) {
  pars$MYZpar$MYZ_ix <- integer(0)
  return(pars)
}

#' @title Make parameters for Ztrace aquatic mosquito model
#' @param pars a [list]
#' @param Zm a vector of mean mosquito densities
#' @param f the blood feeding rate
#' @param q the human fraction
#' @param Zf a [function] of the form Zf(t, pars) that computes temporal fluctuations
#' @return none
#' @export
make_parameters_MYZ_Ztrace <- function(pars, Zm, f, q, Zf) {
  stopifnot(is.numeric(Zm))
  MYZpar <- list()
  class(MYZpar) <- 'Ztrace'
  xde <- "trace"
  class(xde) <- "trace"
  MYZpar$xde <- xde
  MYZpar$Zm <- Zm
  MYZpar$f0 <- f
  MYZpar$q0 <- q
  MYZpar$Zf = Zf
  pars$MYZpar <- MYZpar
  pars = MosquitoBehavior(pars)
  return(pars)
}

#' @title Make parameters for Ztrace aquatic mosquito model
#' @param pars a [list]
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

