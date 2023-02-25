# generic methods for adult component

#' @title Compute bloodfeeding and mortality rates
#' @description This method dispatches on the type of `pars$MYZpar`. It should,
#' at a minimum return the values `f`, `q`, `g` (blood feeding rate, human feeding
#' proportion, and mortality rate) at the current time, although it may return
#' vectors of these values at multiple times for models with delay. These baseline
#' values will be modified by the vector control component. The return type is a
#' named list with those 3 values, and `f`  should have an [attr] labeled `time`
#' giving the time(s) in the simulation that these bionomic values correspond to.
#' @param t current simulation time
#' @param y state vector
#' @param pars an [environment]
#' @return a [list]
#' @export
MosquitoBehavior <- function(t, y, pars) {
  UseMethod("MosquitoBehavior", pars$MYZpar)
}

#' @title Time spent host seeking/feeding and resting/ovipositing
#' @description This method dispatches on the type of `pars$MYZpar`.
#' @param t current simulation time
#' @param y state vector
#' @param pars an [environment]
#' @return either a [numeric] vector if the model supports this feature, or [NULL]
#' @export
F_tau <- function(t, y, pars) {
  UseMethod("F_tau", pars$MYZpar)
}


#' @title Density of infectious mosquitoes
#' @description This method dispatches on the type of `pars$MYZpar`.
#' @param t current simulation time
#' @param y state vector
#' @param pars an [environment]
#' @return a [numeric] vector of length `nPatches`
#' @export
F_Z <- function(t, y, pars) {
  UseMethod("F_Z", pars$MYZpar)
}

#' @title Density of lagged infectious mosquitoes
#' @description This method dispatches on the type of `pars$MYZpar`.
#' @param t current simulation time
#' @param y state vector
#' @param pars an [environment]
#' @param lag duration of lag `t-lag`
#' @return a [numeric] vector of length `nPatches`
#' @export
F_Z_lag <- function(t, y, pars, lag) {
  UseMethod("F_Z_lag", pars$MYZpar)
}

#' @title Number of eggs laid by adult mosquitoes
#' @description This method dispatches on the type of `pars$MYZpar`.
#' @param t current simulation time
#' @param y state vector
#' @param pars an [environment]
#' @return a [numeric] vector of length `nPatches`
#' @export
F_eggs <- function(t, y, pars) {
  UseMethod("F_eggs", pars$MYZpar)
}

#' @title Derivatives for adult mosquitoes
#' @description This method dispatches on the type of `pars$MYZpar`.
#' @param t current simulation time
#' @param y state vector
#' @param pars an [environment]
#' @param Lambda emergence rate of adult mosquitoes
#' @param kappa net infectiousness of human population
#' @param MosyBehavior values returned by [exDE::MosquitoBehavior], potentially modified by control [exDE::VectorControl]
#' @return a [numeric] vector
#' @export
dMYZdt <- function(t, y, pars, Lambda, kappa, MosyBehavior) {
  UseMethod("dMYZdt", pars$MYZpar)
}

#' @title Add indices for adult mosquitoes to parameter list
#' @description This method dispatches on the type of `pars$MYZpar`.
#' @param pars an [environment]
#' @return none
#' @export
make_indices_MYZ <- function(pars) {
  UseMethod("make_indices_MYZ", pars$MYZpar)
}

#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$MYZpar`.
#' @param pars an [environment]
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
