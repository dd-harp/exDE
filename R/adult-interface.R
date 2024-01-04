# generic methods for adult component

#' @title Set bloodfeeding and mortality rates to baseline
#' @description This method dispatches on the type of `pars$MYZpar`. It should
#' set the values of the bionomic parameters to baseline values.
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @param s the species index
#' @return a [list]
#' @export
MBionomics <- function(t, y, pars, s) {
  UseMethod("MBionomics", pars$MYZpar[[s]])
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

#' @title Blood feeding rate of the infective mosquito population
#' @description This method dispatches on the type of `pars$MYZpar`.
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @param s the species index
#' @return a [numeric] vector of length `nPatches`
#' @export
F_fqZ <- function(t, y, pars, s) {
  UseMethod("F_fqZ", pars$MYZpar[[s]])
}

#' @title Number of eggs laid by adult mosquitoes
#' @description This method dispatches on the type of `pars$MYZpar`.
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @param s the species index
#' @return a [numeric] vector of length `nPatches`
#' @export
F_eggs <- function(t, y, pars, s) {
  UseMethod("F_eggs", pars$MYZpar[[s]])
}


#' @title Derivatives for adult mosquitoes
#' @description This method dispatches on the type of `pars$MYZpar`.
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @param s the species index
#' @return the derivatives a [vector]
#' @export
dMYZdt <- function(t, y, pars, s) {
  UseMethod("dMYZdt", pars$MYZpar[[s]])
}

#' @title A function to set up adult mosquito models
#' @description This method dispatches on `MYZname`.
#' @param MYZname the name of the model
#' @param pars a [list]
#' @param s the species index
#' @param MYZopts a [list]
#' @param EIPmod is a [list]
#' @param calK is a [matrix]
#' @return [list]
#' @export
setup_MYZpar = function(MYZname, pars, s, MYZopts=list(), EIPmod=setup_eip_static(), calK=diag(1)){
  class(MYZname) <- MYZname
  UseMethod("setup_MYZpar", MYZname)
}

#' @title A function to set up adult mosquito models
#' @description This method dispatches on `MYZname`.
#' @param pars a [list]
#' @param s the species index
#' @param MYZopts a [list]
#' @return [list]
#' @export
setup_MYZinits = function(pars, s, MYZopts=list()){
  UseMethod("setup_MYZinits", pars$MYZpar[[s]])
}


#' @title Add indices for adult mosquitoes to parameter list
#' @description This method dispatches on the type of `pars$MYZpar`.
#' @param pars a [list]
#' @param s the species index
#' @return [list]
#' @export
make_indices_MYZ <- function(pars, s) {
  UseMethod("make_indices_MYZ", pars$MYZpar[[s]])
}

#' @title Parse the output of deSolve and return the variables by name in a list
#' @description This method dispatches on the type of `pars$MYZpar`.
#' It computes the variables by name and returns a named list.
#' @param deout a [matrix] of outputs from deSolve
#' @param pars a [list] that defines a model
#' @param s the species index
#' @return [list]
#' @export
parse_deout_MYZ <- function(deout, pars, s) {
  UseMethod("parse_deout_MYZ", pars$MYZpar[[s]])
}

#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$MYZpar`.
#' @param pars a [list]
#' @param s the species index
#' @return [numeric]
#' @export
get_inits_MYZ <- function(pars, s) {
  UseMethod("get_inits_MYZ", pars$MYZpar[[s]])
}

#' @title Set the initial values as a vector
#' @description This method dispatches on the type of `pars$MYZpar`.
#' @param pars a [list]
#' @param y0 a vector of variable values from a simulation
#' @param s the species index
#' @return a [list]
#' @export
update_inits_MYZ <- function(pars, y0, s) {
  UseMethod("update_inits_MYZ", pars$MYZpar[[s]])
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


#' @title Convert a model from dde to the corresponding ode
#' @description This method dispatches on the type of `pars$MYZpar$xde`
#' @param pars a [list]
#' @return a [list]
#' @export
dde2ode_MYZ = function(pars){
  UseMethod("dde2ode_MYZ", pars$MYZpar$xde)
}

#' @title Convert a model from dde to the corresponding ode
#' @description If it is already an ode, return pars unchanged.
#' @param pars a [list]
#' @return a [list]
#' @export
dde2ode_MYZ.ode = function(pars){pars}

#' @title Convert a model from dde to the corresponding ode
#' @description If it is a dde, return the corresponding ode
#' @param pars a [list]
#' @return a [list]
#' @export
dde2ode_MYZ.dde = function(pars){
  pars$MYZpar$xde <- "ode"
  pars$MYZpar$solve_as <- "ode"
  pars <- make_MYZpar_RM(pars, MYZopts<- pars$MYZpar,
                         calK=pars$MYZpar$calK)
  pars <- make_indices(pars)
  return(pars)
}
