# Methods to compute mixing and parasite / pathogen transmission during bloood feeding

#' @title Compute the mixing matrix and transmission terms
#' @description This method dispatches on the type of `pars$BFpar`
#' @param t current simulation time
#' @param y state vector
#' @param pars, a [list]
#' @return [list]
#' @export
Transmission <- function(t, y, pars){
  UseMethod("Transmission", pars$BFpar)
}

#' @title Compute transmission terms dynamically, no update required
#' @description This method dispatches on the type of `pars$BFpar`
#' @inheritParams Transmission
#' @return [list]
#' @export
Transmission.static <- function(t, y, pars){

  pars = compute_EIR(t, y, pars)
  pars = compute_kappa(t, y, pars)

  return(pars)
}

#' @title Compute transmission terms, the dynamic case
#' @description This method dispatches on the type of `pars$BFpar`
#' @inheritParams Transmission
#' @return [list]
#' @export
Transmission.dynamic <- function(t, y, pars){

  pars = LocalFrac(pars)
  pars = compute_beta(t, y, pars)
  pars = compute_EIR(t, y, pars)
  pars = compute_kappa(t, y, pars)

  return(pars)
}

#' @title Compute transmission terms, the dynamic case
#' @description This method dispatches on the type of `pars$BFpar`
#' @inheritParams Transmission
#' @return [list]
#' @export
Transmission.setup <- function(t, y, pars){

  pars = setup_local_fraction_simple(pars)
  for(s in 1:pars$nVectors)
    pars = compute_AvailableHosts(t, y, pars, s)
  pars = compute_beta(t, y, pars)
  pars = compute_EIR(t, y, pars)
  pars = compute_kappa(t, y, pars)

  class(pars$BFpar) <- "static"
  return(pars)
}

#' @title Compute beta
#' @description This function computes the mixing matrix, beta
#' @param t current simulation time
#' @param y state vector
#' @param pars, a [list]
#' @return [list]
#' @export
compute_beta <- function(t, y, pars){

  for(i in 1:pars$nHosts){
    H = F_H(t, y, pars, i)
    for(s in 1:pars$nVectors){
      W = pars$vars$W[[s]]
      wts = pars$BFpar$searchWts[[i]][[s]]
      TaR = pars$BFpar$TaR[[i]][[s]]
      pars$beta[[i]][[s]] <- F_beta(H, W, wts, TaR)
    }
  }

  return(pars)
}

#' @title Compute beta, the biting distribution matrix
#' @param H human / host population density
#' @param W human / host availability in the patches
#' @param wts_f the blood feeding search weights
#' @param TaR (time at risk), a [matrix]  dimensions `nPatches` by `nStrata`
#' @return a [matrix] of dimensions `nStrata` by `nPatches`
#' @export
F_beta = function(H, W, wts_f, TaR){
  beta = diag(wts_f, length(H)) %*% t(TaR) %*% diag(1/W, length(W))
  return(beta)
}

#' @title Compute EIR
#' @description This function computes the EIR for each stratum of each host species
#' @param t current simulation time
#' @param y state vector
#' @param pars, a [list]
#' @return [list]
#' @export
compute_EIR <- function(t, y, pars){

  for(i in 1:pars$nHosts){
    pars$EIR[[i]] <- F_EIR(t, y, pars, i, 1)
    if(pars$nVectors > 1)
      for(s in 2:pars$nVectors){
        pars$EIR[[i]] <- pars$EIR[[i]] + F_EIR(t, y, pars, i, s)
    }
  }

  return(pars)
}

#' @title Compute EIR for each vector-host pair
#' @description This function computes the EIR from each vector species for each stratum of each host species
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @return [list]
#' @export
compute_EIR_full <- function(t, y, pars){

  for(i in 1:pars$nHosts){
    pars$eir[[i]][[1]] <- F_EIR(t, y, pars, i, 1)
    pars$EIR[[i]] <- pars$EIR[[i]] + pars$eir[[i]][[1]]
    if(s>1)
      for(s in 2:pars$nVectors){
        pars$eir[[i]][[s]] <- F_EIR(t, y, pars, i, s)
        pars$EIR[[i]] <- pars$EIR[[i]] + pars$eir[[i]][[s]]
    }
  }

  return(pars)
}

#' @title Entomological inoculation rate on human strata
#' @description Compute the daily EIR for all the strata at time `t`
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @param i the host species index
#' @param s the vector species index
#' @return [numeric] vector of length `nStrata`
#' @export
F_EIR <- function(t, y, pars, i, s) {

  fqZ <- F_fqZ(t, y, pars, s)
  beta <- pars$beta[[i]][[s]]
  lf <- pars$vars$local_frac[[s]]
  eir = beta %*% (fqZ*lf)

  return(as.vector(eir))
}

#' @title Compute kappa
#' @description This function computes kappa for each vector species in each patch
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @return [list]
#' @export
compute_kappa <- function(t, y, pars){
  for(s in 1:pars$nVectors){
    kappa <- F_kappa(t, y, pars, 1, s)

    if(pars$nHosts>1)
      for(i in 2:pars$nHosts)
        kappa <- kappa + F_kappa(t, y, pars, i, s)

    kappa = with(pars$vars, local_frac[[s]]*kappa + (1-local_frac[[s]])*x_visitors[[s]])
    pars$kappa[[s]] = kappa
  }

  return(pars)
}

#' @title Net infectiousness of human population to mosquitoes
#' @description Compute the net infectiousness of humans in each patch at time `t`
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @param i the host species index
#' @param s the vector species index
#' @return a [numeric] vector of length `nPatches`
#' @export
F_kappa <- function(t, y, pars, i, s) {

  beta = pars$beta[[i]][[s]]
  Wi = pars$vars$Wi[[i]][[s]]
  W = pars$vars$W[[s]]
  kappa = Wi/W*(as.vector(t(beta) %*% F_X(t, y, pars, i)))

  return(as.vector(kappa))
}
