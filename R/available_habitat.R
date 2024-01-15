
#' @title Update the availability of aquatic habitats
#' @description This method dispatches on the type of `pars$EGGpar`.
#' @param pars an [list]
#' @return a [list]
#' @export
AvailableHabitat <- function(pars) {
  UseMethod("AvailableHabitat", pars$EGGpar)
}

#' @title Compute total availability of aquatic habitats
#' @description Computes the availability of aquatic habitats for the static model (do nothing)
#' @param pars a [list]
#' @return a [numeric] vector of length `nPatches`
#' @export
AvailableHabitat.static <- function(pars){
  return(pars)
}

#' @title Compute total availability of aquatic habitats
#' @description Computes the availability of aquatic habitats for the static model (do nothing)
#' @param pars a [list]
#' @return a [numeric] vector of length `nPatches`
#' @export
AvailableHabitat.simple <- function(pars){
  return(pars)
}

#' @title Compute total availability of aquatic habitats
#' @description Computes the availability of aquatic habitats
#' @param pars a [list]
#' @return a [numeric] vector of length `nPatches`
#' @export
AvailableHabitat.forced <- function(pars){
  for(s in 1:pars$nVectors) pars = compute_AvailableHabitat(pars, s)
  return(pars)
}

#' @title Compute total availability of aquatic habitats
#' @description Computes the availability of aquatic habitats
#' @param pars a [list]
#' @param s the vector species index
#' @return a [numeric] vector of length `nPatches`
#' @export
compute_AvailableHabitat <- function(pars, s){
  habitats = with(pars, calN %*% EGGpar$searchWts[[s]])
  Q = habitats + pars$vars$ovitraps + pars$vars$non_habitats
  Qfrac = as.vector(habitats/Q)
  Qfrac[which(Q==0)]=0
  pars$vars$Q[[s]] = as.vector(Q)
  pars$vars$Qfrac[[s]] =  as.vector(Qfrac)
  return(pars)
}
