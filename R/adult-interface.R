# generic methods for adult component

#' @title Set bloodfeeding and mortality rates to baseline
#' @description This method dispatches on the type of `pars$MYZpar`. It should
#' set the values of the bionomic parameters to baseline values.
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @return a [list]
#' @export
MosquitoBehavior <- function(t, y, pars) {
  UseMethod("MosquitoBehavior", pars$MYZpar)
}

#' @title Time spent host seeking/feeding and resting/ovipositing
#' @description This method dispatches on the type of `pars$MYZpar`.
#' @param t current simulation time
#' @param pars a [list]
#' @return either a [numeric] vector if the model supports this feature, or [NULL]
#' @export
F_tau <- function(t, pars) {
  UseMethod("F_tau", pars$MYZpar)
}

#' @title Density of infectious mosquitoes
#' @description This method dispatches on the type of `pars$MYZpar`.
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @return a [numeric] vector of length `nPatches`
#' @export
F_Z <- function(t, y, pars) {
  UseMethod("F_Z", pars$MYZpar)
}

#' @title Number of eggs laid by adult mosquitoes
#' @description This method dispatches on the type of `pars$MYZpar`.
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @return a [numeric] vector of length `nPatches`
#' @export
F_eggs <- function(t, y, pars) {
  UseMethod("F_eggs", pars$MYZpar)
}


#' @title Derivatives for adult mosquitoes
#' @description This method dispatches on the type of `pars$MYZpar`.
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @param Lambda emergence rate of adult mosquitoes
#' @param kappa net infectiousness of human population
#' @return a [numeric] vector
#' @export
dMYZdt <- function(t, y, pars, Lambda, kappa) {
  UseMethod("dMYZdt", pars$MYZpar)
}

#' @title A function to set up adult mosquito models
#' @description This method dispatches on `MYZname`.
#' @param pars a [list]
#' @param MYZname a [character] string
#' @param nPatches a [numeric] value (an integer), the number of patches
#' @param MYZopts a [list]
#' @param calK is a [matrix]
#' @return none
#' @export
setup_MYZ = function(pars, MYZname, nPatches=1, MYZopts=list(), calK=diag(1)){
  class(MYZopts) <- MYZname
  UseMethod("setup_MYZ", MYZopts)
}

#' @title Add indices for adult mosquitoes to parameter list
#' @description This method dispatches on the type of `pars$MYZpar`.
#' @param pars a [list]
#' @return none
#' @export
make_indices_MYZ <- function(pars) {
  UseMethod("make_indices_MYZ", pars$MYZpar)
}

#' @title Parse the output of deSolve and return the variables by name in a list
#' @description This method dispatches on the type of `pars$MYZpar`. Adds the variables
#' from the MYZ model to varslist and returns it
#' @param varslist a [list] the object to be returned
#' @param deout a [matrix] of outputs from deSolve
#' @param pars a [list] that defines a model
#' @export
parse_deout_MYZ <- function(varslist, deout, pars) {
  UseMethod("parse_deout_MYZ", pars$MYZpar)
}

#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$MYZpar`.
#' @param pars a [list]
#' @return none
#' @export
get_inits_MYZ <- function(pars) {
  UseMethod("get_inits_MYZ", pars$MYZpar)
}

#' @title Make the mosquito demography matrix
#' @param g mortality rate
#' @param sigma emigration  rate
#' @param K mosquito dispersal matrix
#' @param nPatches number of patches
#' @return a [matrix] of dimensions `nPatches` by `nPatches`
#' @export
make_Omega <- function(g, sigma, K, nPatches) {
  diag(g, nPatches) + ((diag(nPatches) - K) %*% diag(sigma, nPatches))
}

