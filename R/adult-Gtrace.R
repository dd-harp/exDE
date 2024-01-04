# specialized methods for the adult mosquito Gtrace model

#' @title Compute bloodfeeding and mortality rates
#' @description Implements [MBionomics] for the Gtrace model.
#' @inheritParams MBionomics
#' @return a named [list]
#' @export
MBionomics.Gtrace <- function(t, y, pars, s) {
  return(pars)
}

#' @title Blood feeding rate of the infective mosquito population
#' @description Implements [F_fqZ] for the Gtrace model.
#' @inheritParams F_fqZ
#' @return a [numeric] vector of length `nHabitats`
#' @export
F_fqZ.Gtrace <- function(t, y, pars, s) {
  numeric(0)
}

#' @title Number of eggs laid by adult mosquitoes
#' @description Implements [F_eggs] for the Gtrace model.
#' @inheritParams F_eggs
#' @return a [numeric] vector of length `nPatches`
#' @export
F_eggs.Gtrace <- function(t, y, pars, s) {
  with(pars$MYZpar[[s]], return(Gm*Gf(t, pars$MYZpar[[s]])))
}

#' @title Derivatives for aquatic stage mosquitoes
#' @description Implements [dMYZdt] for the Gtrace (forced emergence) model.
#' @inheritParams dMYZdt
#' @return a [numeric] vector
#' @export
dMYZdt.Gtrace <- function(t, y, pars, s){
  numeric(0)
}


#' @title Setup the Gtrace
#' @description Implements [setup_MYZpar] for the Gtrace model
#' @inheritParams setup_MYZpar
#' @return a [list] vector
#' @export
setup_MYZpar.Gtrace = function(MYZname, pars, s, MYZopts=NULL, EIPmod=NULL, calK=NULL){
  pars$MYZpar[[s]] = make_MYZpar_Gtrace(pars$nPatches, MYZopts)
  return(pars)
}


#' @title Make parameters for Gtrace aquatic mosquito model
#' @param nPatches the number of patches in the model
#' @param MYZopts a [list] to overwrite the defaults
#' @param Gm a vector of mean mosquito densities
#' @param Gf a [function] of the form Gf(t, pars) that computes temporal fluctuations
#' @return none
#' @export
make_MYZpar_Gtrace = function(nPatches, MYZopts, Gm = 1, Gf=NULL){
  with(MYZopts,{
    MYZpar <- list()
    class(MYZpar) <- "Gtrace"

    MYZpar$Gm <- checkIt(Gm, nPatches)
    if(is.null(Gf)) Gf = function(t, MYZpar){return(1)}
    MYZpar$Gf = Gf

    return(MYZpar)
})}

#' @title Setup the Gtrace model
#' @description Implements [setup_MYZinits] for the Gtrace model
#' @inheritParams setup_MYZinits
#' @return a [list] vector
#' @export
setup_MYZinits.Gtrace = function(pars, s, MYZopts=NULL){
  return(pars)
}

#' @title Add indices for aquatic stage mosquitoes to parameter list
#' @description Implements [make_indices_MYZ] for Gtrace (forced emergence) model.
#' @inheritParams make_indices_MYZ
#' @return none
#' @export
make_indices_MYZ.Gtrace <- function(pars, s) {
  return(pars)
}

#' @title Parse the output of deSolve and return variables for the Gtrace model
#' @description Implements [parse_deout_MYZ] for Gtrace
#' @inheritParams parse_deout_MYZ
#' @return [list]
#' @export
parse_deout_MYZ.Gtrace <- function(deout, pars, s) {
  return(list())
}

#' @title Make parameters for Gtrace aquatic mosquito model
#' @param pars a [list]
#' @param Gm a vector of mean mosquito densities
#' @param Gf a [function] of the form Gf(t, pars) that computes temporal fluctuations
#' @return none
#' @export
make_parameters_MYZ_Gtrace <- function(pars, Gm, Gf) {
  stopifnot(is.numeric(Gm))
  MYZpar <- list()
  class(MYZpar) <- 'Gtrace'
  xde <- "trace"
  class(xde) <- "trace"
  MYZpar$xde <- xde
  MYZpar$Gm <- Gm
  MYZpar$Gf = Gf
  pars$MYZpar[[1]] <- MYZpar
  return(pars)
}

#' @title Make parameters for Gtrace aquatic mosquito model
#' @param pars a [list]
#' @param MYZ0 is set to NULL for the Gtrace model
#' @return none
#' @export
make_inits_MYZ_Gtrace<- function(pars, MYZ0=NULL) {
  pars$MYZinits[[1]] = numeric(0)
  return(pars)
}

#' @title Update inits for Gtrace
#' @inheritParams update_inits_MYZ
#' @return none
#' @export
update_inits_MYZ.Gtrace <- function(pars, y0, s) {
  return(pars)
}

#' @title Return initial values as a vector
#' @description Implements [get_inits_MYZ] for the GeRM model.
#' @inheritParams get_inits_MYZ
#' @return none
#' @export
get_inits_MYZ.Gtrace <- function(pars, s){
  return(c())
}

