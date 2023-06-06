# Functions that exogenously force the availability of resources

#' @title Compute host availability
#' @description Compute the availability for the pathogen's hosts for blood feeding
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @return a [numeric] vector of length `nPatches`
#' @export
compute_W <- function(t, y, pars){
  H = F_H(t, y, pars)
  W = with(pars$Hpar, TaR %*% (wts_f*H))
  return(W)
}

#' @title Compute total availability of aquatic habitats
#' @description Computes the availability of aquatic habitats
#' @param pars a [list]
#' @return a [numeric] vector of length `nPatches`
#' @export
compute_Q <- function(pars){
  with(pars,
    return(calN %*% searchQ)
)}

#' @title Compute the total availability of blood hosts
#' @description Compute the availability for the pathogen's hosts for blood feeding
#' @param pars a [list]
#' @return a [numeric] vector of length `nPatches`
#' @export
compute_B <- function(pars){
  B = with(pars, W + OtherBloodHosts + Visitors)
  return(B)
}

#' @title Compute the fraction of bites occurring on local hosts
#' @description Computes the local_frac
#' @param pars a [list]
#' @return a [numeric] vector of length `nPatches`
#' @export
compute_local_frac <- function(pars){
  with(pars,
       return(W/(W+Visitors))
)}


