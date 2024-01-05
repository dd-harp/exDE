# Compute the Transmission Terms

#' @title Compute the biting distribution matrix, beta
#' @description This method dispatches on the type of `pars$xde`
#' @param t current simulation time
#' @param y state vector
#' @param pars, a [list]
#' @return [list]
#' @export
Transmission <- function(t, y, pars){

  # beta: mixing
  for(i in 1:pars$nHosts)
    for(s in 1:pars$nVectors)
      pars$beta[[i]][[s]] <- F_beta(t, y, pars, s, i)

  # EIR: entomological inoculation rate
  for(i in 1:pars$nHosts){
    pars$EIR[[i]] <- 0*F_EIR(t, y, pars, pars$beta[[i]][[1]], 1)
    for(s in 1:pars$nVectors){
#      pars$eir[[i]][[s]] <- F_EIR(t, y, pars, pars$beta[[i]][[s]], s)
#      pars$EIR[[i]] <- pars$EIR[[i]] + pars$eir[[i]][[s]]
      pars$EIR[[i]] <- pars$EIR[[i]] + F_EIR(t, y, pars, pars$beta[[i]][[s]], s)
    }
    pars$EIR[[i]] = pars$EIR[[i]]*pars$local_frac
  }

  for(s in 1:pars$nVectors){
    pars$kappa[[s]] <- 0* F_kappa(t, y, pars, pars$beta[[1]][[s]], 1)
    for(i in 1:pars$nHosts){
#      pars$ni[[i]][[s]] <- F_kappa(t, y, pars, pars$beta[[i]][[s]], i)
#      pars$kappa[[s]] <- pars$kappa[[s]] + W[[i]]/Wtot*pars$ni[[i]][[s]]
      pars$kappa[[s]] <- pars$kappa[[s]] + F_kappa(t, y, pars, pars$beta[[i]][[s]], i)
    }
    pars$kappa[[s]] = with(pars, local_frac*kappa[[s]] + (1-local_frac)*x_visitors)
  }
  return(pars)
}

#' @title Compute the biting distribution matrix, beta
#' @description This method dispatches on the type of `pars$xde`
#' @param t current simulation time
#' @param y state vector
#' @param pars, a [list]
#' @param s the vector species index
#' @param i the host species index
#' @return [list]
#' @export
F_beta <- function(t, y, pars, s, i){
  H <- F_H(t, y, pars, i)
  beta = compute_beta(H, pars$Hpar[[i]]$wts_f, pars$Hpar[[i]]$TaR)
  return(beta)
}

#' @title Entomological inoculation rate on human strata
#' @description Compute the daily EIR for all the strata at time `t`
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @param beta, a [matrix]
#' @param s the vector species index
#' @return [numeric] vector of length `nStrata`
#' @export
F_EIR <- function(t, y, pars, beta, s) {
  fqZ <- F_fqZ(t, y, pars, s)
  return(as.vector(beta %*% fqZ))
}

#' @title Net infectiousness of human population to mosquitoes
#' @description Compute the net infectiousness of humans in each patch at time `t`
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @param beta, a [matrix]
#' @param i the host species index
#' @return a [numeric] vector of length `nPatches`
#' @export
F_kappa <- function(t, y, pars, beta, i) {
  if(i>1) browser()
  kappa = as.vector(t(beta) %*% F_X(t, y, pars, i))
#  kappa = pars$HostAvailability[[i]]/pars$Wtot*as.vector(t(beta) %*% F_X(t, y, pars, i))
  return(kappa)
}

#' @title Compute beta, the biting distribution matrix
#' @param H the human population size
#' @param wts_f the blood feeding search weights
#' @param TaR (time at risk), a [matrix]  dimensions `nPatches` by `nStrata`
#' @return a [matrix] of dimensions `nStrata` by `nPatches`
#' @export
compute_beta = function(H, wts_f, TaR){
  W <- as.vector(TaR %*% (wts_f*H))
  beta <- diag(wts_f, length(H)) %*% t(TaR) %*% diag(1/W, length(W))
  return(beta)
}
