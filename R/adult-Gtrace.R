# specialized methods for the adult mosquito Gtrace model

#' @title Compute bloodfeeding and mortality rates
#' @description Implements [MosquitoBehavior] for the Gtrace model.
#' @inheritParams MosquitoBehavior
#' @return a named [list]
#' @export
MosquitoBehavior.Gtrace <- function(t, y, pars) {
  return(pars)
}

#' @title Number of infective adults in each patch
#' @description Implements [F_Z] for the Gtrace  model.
#' @inheritParams F_Z
#' @return a [numeric] vector of length `nHabitats`
#' @export
F_Z.Gtrace <- function(t, y, pars) {
  return(numeric(0))
}

#' @title Number of eggs laid by adult mosquitoes
#' @description Implements [F_eggs] for the Gtrace model.
#' @inheritParams F_eggs
#' @return a [numeric] vector of length `nPatches`
#' @export
F_eggs.Gtrace <- function(t, y, pars) {
  with(pars$MYZpar, return(Gm*Gf(t, pars)))
}

#' @title Derivatives for aquatic stage mosquitoes
#' @description Implements [dMYZdt] for the Gtrace (forced emergence) model.
#' @inheritParams dMYZdt
#' @return a [numeric] vector
#' @export
dMYZdt.Gtrace <- function(t, y, pars, Lambda, kappa){
  numeric(0)
}


#' @title Setup the Gtrace
#' @description Implements [setup_MYZ] for the Gtrace model
#' @inheritParams setup_MYZ
#' @return a [list] vector
#' @export
setup_MYZ.Gtrace = function(pars, MYZname,
                               nPatches=1, MYZopts=NULL,
                               calK=diag(1)){

  pars$MYZname = "Gtrace"
  pars$nPatches = checkIt(nPatches, 1, "integer")

  pars = make_MYZpar_Gtrace(pars, MYZopts)
  pars$MYZinits = numeric(0)

  return(pars)
}


#' @title Make parameters for Gtrace aquatic mosquito model
#' @param pars a [list]
#' @param MYZopts a [list] to overwrite the defaults
#' @param Gm a vector of mean mosquito densities
#' @param Gf a [function] of the form Gf(t, pars) that computes temporal fluctuations
#' @return none
#' @export
make_MYZpar_Gtrace = function(pars, MYZopts,
                              Gm = 1, Gf=NULL){
  with(MYZopts,{
    MYZpar <- list()
    class(MYZpar) <- "Gtrace"

    MYZpar$Gm <- checkIt(Gm, pars$nPatches)
    if(is.null(Gf)) Gf = function(t, y, pars){return(1)}
    MYZpar$Gf = Gf

    pars$MYZpar = MYZpar
    return(pars)
})}

#' @title Add indices for aquatic stage mosquitoes to parameter list
#' @description Implements [make_indices_MYZ] for Gtrace (forced emergence) model.
#' @inheritParams make_indices_MYZ
#' @return none
#' @export
make_indices_MYZ.Gtrace <- function(pars) {
  pars$MYZpar$MYZ_ix <- integer(0)
  return(pars)
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
  pars$MYZpar <- MYZpar
  pars = MosquitoBehavior(pars)
  return(pars)
}

#' @title Make parameters for Gtrace aquatic mosquito model
#' @param pars a [list]
#' @param MYZ0 is set to NULL for the Gtrace model
#' @return none
#' @export
make_inits_MYZ_Gtrace<- function(pars, MYZ0=NULL) {
  pars$MYZinits = numeric(0)
  return(pars)
}

#' @title Return initial values as a vector
#' @description Implements [get_inits_MYZ] for the GeRM model.
#' @inheritParams get_inits_MYZ
#' @return none
#' @export
get_inits_MYZ.Gtrace <- function(pars){
  numeric(0)
}

