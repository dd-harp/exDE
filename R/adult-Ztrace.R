# specialized methods for the adult mosquito Ztrace model

#' @title Compute bloodfeeding and mortality rates
#' @description Implements [MBionomics] for the Ztrace model.
#' @inheritParams MBionomics
#' @return a named [list]
#' @export
MBionomics.Ztrace <- function(t, y, pars, s) {

  with(pars$MYZpar[[s]],{
    pars$MYZpar[[s]]$f <- f0
    pars$MYZpar[[s]]$q <- q0
    return(pars)
})}

#' @title Blood feeding rate of the infective mosquito population
#' @description Implements [F_fqZ] for the Ztrace model.
#' @inheritParams F_fqZ
#' @return a [numeric] vector of length `nHabitats`
#' @export
F_fqZ.Ztrace <- function(t, y, pars, s) {
  with(pars$MYZpar[[s]], return(f*q*scale*Zf(t)))
}

#' @title Blood feeding rate of the infective mosquito population
#' @description Implements [F_fqM] for the Ztrace model.
#' @inheritParams F_fqM
#' @return a [numeric] vector of length `nHabitats`
#' @export
F_fqM.Ztrace <- function(t, y, pars, s) {
  numeric(0)
}

#' @title Number of eggs laid by adult mosquitoes
#' @description Implements [F_eggs] for the Ztrace model.
#' @inheritParams F_eggs
#' @return a [numeric] vector of length `nPatches`
#' @export
F_eggs.Ztrace <- function(t, y, pars, s) {
  return(numeric(0))
}

#' @title Derivatives for aquatic stage mosquitoes
#' @description Implements [dMYZdt] for the Ztrace (forced emergence) model.
#' @inheritParams dMYZdt
#' @return a [numeric] vector
#' @export
dMYZdt.Ztrace <- function(t, y, pars, s){
  numeric(0)
}

#' @title Setup the Ztrace model
#' @description Implements [setup_MYZpar] for the Ztrace model
#' @inheritParams setup_MYZpar
#' @return a [list] vector
#' @export
setup_MYZpar.Ztrace = function(MYZname, pars, s, MYZopts=NULL, EIPmod=NULL, calK=NULL){
  pars$MYZpar[[s]] = make_MYZpar_Ztrace(pars, MYZopts)
  return(pars)
}


#' @title Make parameters for Ztrace aquatic mosquito model
#' @param nPatches an integer
#' @param MYZopts a [list] of values that overwrites the defaults
#' @param f the blood feeding rate
#' @param q the human fraction
#' @param Zm a vector of mean mosquito densities
#' @param Zf a [function] of the form Zf(t, pars) that computes temporal fluctuations
#' @return none
#' @export
make_MYZpar_Ztrace = function(nPatches, MYZopts,
                              f = 0.3, q = 0.95,
                              Zm = 1, Zf = NULL){

  with(MYZopts,{
    MYZpar <- list()
    class(MYZpar) <- "Ztrace"

    MYZpar$f <- checkIt(f, pars$nPatches)
    MYZpar$q <- checkIt(q, pars$nPatches)
    MYZpar$f0 <- MYZpar$f
    MYZpar$q0 <- MYZpar$q

    MYZpar$scale <- checkIt(Zm, pars$nPatches)
    if(is.null(Zf)) Zf = function(t){return(1)}
    MYZpar$Zf <- Zf

    return(MYZpar)
  })}


#' @title Setup the Ztrace model
#' @description Implements [setup_MYZinits] for the Ztrace model
#' @inheritParams setup_MYZinits
#' @return a [list] vector
#' @export
setup_MYZinits.Ztrace = function(pars, s, MYZopts=NULL){
  return(numeric(0))
}

#' @title Parse the output of deSolve and return variables for the Ztrace model
#' @description Implements [parse_deout_MYZ] for Ztrace
#' @inheritParams parse_deout_MYZ
#' @return a [list]
#' @export
parse_deout_MYZ.Ztrace <- function(deout, pars, s) {
  return(list())
}

#' @title Add indices for aquatic stage mosquitoes to parameter list
#' @description Implements [make_indices_MYZ] for Ztrace (forced emergence) model.
#' @inheritParams make_indices_MYZ
#' @return none
#' @export
make_indices_MYZ.Ztrace <- function(pars, s) {
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
make_parameters_MYZ_Ztrace <- function(pars, Zm, f, q, Zf=NULL) {
  stopifnot(is.numeric(Zm))
  MYZpar <- list()
  class(MYZpar) <- 'Ztrace'
  xde <- "trace"
  class(xde) <- "trace"
  MYZpar$xde <- xde
  MYZpar$f0 <- f
  MYZpar$f <- f
  MYZpar$q0 <- q
  MYZpar$q <- q
  MYZpar$scale <- checkIt(Zm, pars$nPatches)
  if(is.null(Zf)) Zf = function(t){return(1)}
  MYZpar$Zf = Zf
  pars$MYZpar <- MYZpar
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

#' @title Update inits for Ztrace
#' @inheritParams update_inits_MYZ
#' @return none
#' @export
update_inits_MYZ.Ztrace <- function(pars, y0, s) {
  return(pars)
}

#' @title Return initial values as a vector
#' @description Implements [get_inits_MYZ] for the GeRM model.
#' @inheritParams get_inits_MYZ
#' @return none
#' @export
get_inits_MYZ.Ztrace <- function(pars, s){
  numeric(0)
}

