
#' @title Update the availability of blood hosts
#' @description This method dispatches on the type of `pars$BFpar`.
#' @param t current simulation time
#' @param y vector of state variables
#' @param pars an [list]
#' @return a [list]
#' @export
AvailableBlood <- function(t, y, pars) {
  UseMethod("AvailableBlood", pars$BFpar)
}

#' @title Compute availability of local humans for blood feeding
#' @description Compute the availability for the pathogen's hosts for blood feeding
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @return a [numeric] vector of length `nPatches`
#' @export
AvailableBlood.static <- function(t, y, pars){
  return(pars)
}

#' @title Compute availability of local humans for blood feeding
#' @description Compute the availability for the pathogen's hosts for blood feeding
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @return a [numeric] vector of length `nPatches`
#' @export
AvailableBlood.setup <- function(t, y, pars){
  for(s in 1:pars$nVectors)
    pars = compute_AvailableHosts(t, y, pars, s)

  return(pars)
}

#' @title Compute availability of local humans for blood feeding
#' @description Compute the availability for the pathogen's hosts for blood feeding
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @return a [numeric] vector of length `nPatches`
#' @export
AvailableBlood.forced <- function(t, y, pars){
  for(s in 1:pars$nVectors)
    pars = compute_AvailableHosts(t, y, pars, s)

  return(pars)
}

#' @title Compute availability blood hosts of the i^th species
#' @description Compute the availability for the pathogen's hosts for blood feeding
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @param s the vector species index
#' @return a [numeric] vector of length `nPatches`
#' @export
compute_AvailableHosts <- function(t, y, pars, s){

  H = F_H(t, y, pars, 1)

  pars$vars$Wi[[1]][[s]] = as.vector(with(pars$BFpar, TaR[[1]][[s]] %*% (searchWts[[1]][[s]]*H)))
  pars$vars$W[[s]] = pars$vars$Wi[[1]][[s]]
  if(pars$nHosts > 1){
    for(i in 2:pars$nHosts){
      H = F_H(t, y, pars, i)
      pars$vars$Wi[[i]][[s]] = as.vector(with(pars$BFpar, TaR[[i]][[s]] %*% (searchWts[[i]][[s]]*H)))
      pars$vars$W[[s]] =  pars$vars$W[[s]] + pars$vars$Wi[[i]][[s]]
    }
  }
  pars$vars$W[[s]] = as.vector(pars$vars$W[[s]])
  pars$vars$B[[s]] = as.vector(with(pars$vars, W[[s]] + Visitors[[s]] + Other[[s]]))

  return(pars)
}
