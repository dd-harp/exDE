#' @title Compute the EIR
#' @description Using the output of deSolve,
#' compute the EIR at every point in time
#' @param deout a [matrix], the output of deSolve
#' @param pars a [list]
#' @return [matrix]
#' @export
compute_EIR <- function(deout, pars) {
  UseMethod("compute_EIR", pars$xde)
}

#' @title Compute the EIR
#' @description Using the output of deSolve,
#' compute the EIR at every point in time
#' @param deout a [matrix], the output of deSolve
#' @param pars a [list]
#' @return [matrix]
#' @export
compute_EIR_default <- function(deout, pars) {
  ix = 1:length(deout[,1])
  eir = sapply(ix, compute_EIR_i, deout=deout, pars=pars)
  return(eir)
}

#' @title Compute the EIR for the ith point in time
#' @description Using the output of deSolve,
#' compute the EIR at one point in time
#' @param i an [integer]
#' @param deout a matrix, the output of deSolve
#' @param pars a [list]
#' @return [numeric]
#' @export
compute_EIR_i <- function(i, deout, pars) {
  t = deout[i,1]
  y = deout[i,-1]
  beta <- F_beta(t, y, pars)
  eir = F_EIR(t, y, pars, beta)
  return(eir)
}

#' @title Compute the EIR
#' @description Implements [compute_EIR] for the full DDEs
#' @inheritParams compute_EIR
#' @return a [matrix]
#' @export
compute_EIR.dde <- function(deout, pars) {
  compute_EIR_default(deout, pars)
}

#' @title Compute the EIR
#' @description Implements [compute_EIR] for the full ODEs
#' @inheritParams compute_EIR
#' @return a [matrix]
#' @export
compute_EIR.ode <- function(deout, pars) {
  compute_EIR_default(deout, pars)
}

#' @title Compute the EIR
#' @description Implements [compute_EIR] for aquatic ODEs
#' @inheritParams compute_EIR
#' @return a [matrix]
#' @export
compute_EIR.aqua <- function(deout, pars) {
  return(numeric(0))
}

#' @title Compute the EIR
#' @description Implements [compute_EIR] for aquatic DDEs
#' @inheritParams compute_EIR
#' @return a [vector]
#' @export
compute_EIR.aqua_dde <- function(deout, pars) {
  return(numeric(0))
}

#' @title Compute the EIR
#' @description Implements [compute_EIR] for mosquito ODEs
#' @inheritParams compute_EIR
#' @return a [vector]
#' @export
compute_EIR.mosy <- function(deout, pars) {
  return(numeric(0))
}

#' @title Compute the EIR
#' @description Implements [compute_EIR] for mosquito DDEs
#' @inheritParams compute_EIR
#' @return a [vector]
#' @export
compute_EIR.mosy_dde <- function(deout, pars) {
  return(numeric(0))
}

#' @title Compute the EIR
#' @description Implements [compute_EIR] for human ODEs
#' @inheritParams compute_EIR
#' @return a [vector]
#' @export
compute_EIR.human<- function(deout, pars) {
  compute_EIR_default(deout, pars)
}

#' @title Compute the EIR
#' @description Implements [compute_EIR] for cohort ODEs
#' @inheritParams compute_EIR
#' @return a [vector]
#' @export
compute_EIR.cohort<- function(deout, pars) {
  t = deout[,1]
  eir = pars$F_eir(t, pars)
}

#' @title Compute the fqZ
#' @description Using the output of [deSolve::ode] or [deSolve::dede],
#' compute the fqZ at every point in time
#' @param deout a matrix, the output of deSolve
#' @param pars a [list]
#' @return [vector]
#' @export
compute_fqZ <- function(deout, pars) {
  UseMethod("compute_fqZ", pars$xde)
}

#' @title Compute the fqZ
#' @description Using the output of [deSolve::ode] or [deSolve::dede],
#' compute the fqZ at every point in time
#' @param deout a matrix, the output of deSolve
#' @param pars a [list]
#' @return [vector]
#' @export
compute_fqZ_default <- function(deout, pars) {
  ix = 1:length(deout[,1])
  fqZ = sapply(ix, compute_fqZ_i, deout=deout, pars=pars)
  return(fqZ)
}

#' @title Compute the fqZ for the ith
#' @description Using the output of [deSolve::ode] or [deSolve::dede],
#' compute the fqZ at one point in time
#' @param i an [integer]
#' @param deout a matrix, the output of deSolve
#' @param pars a [list]
#' @return [numeric]
#' @export
compute_fqZ_i <- function(i, deout, pars) {
  t = deout[i,1]
  y = deout[i,-1]
  fqZ = F_fqZ(t, y, pars)
  return(fqZ)
}

#' @title Compute the fqZ
#' @description Implements [compute_fqZ] for class(pars$xde)=dde
#' @inheritParams compute_fqZ
#' @return a [vector]
#' @export
compute_fqZ.dde <- function(deout, pars) {
  compute_fqZ_default(deout, pars)
}

#' @title Compute the fqZ
#' @description Implements [compute_fqZ] for class(pars$xde)=ode
#' @inheritParams compute_fqZ
#' @return a [vector]
#' @export
compute_fqZ.ode <- function(deout, pars) {
  compute_fqZ_default(deout, pars)
}

#' @title Compute the fqZ
#' @description Implements [compute_fqZ] for class(pars$xde)=aqua
#' @inheritParams compute_fqZ
#' @return a [vector]
#' @export
compute_fqZ.aqua <- function(deout, pars) {
  return(numeric(0))
}

#' @title Compute the fqZ
#' @description Implements [compute_fqZ] for class(pars$xde)=aqua_dde
#' @inheritParams compute_fqZ
#' @return a [vector]
#' @export
compute_fqZ.aqua_dde <- function(deout, pars) {
  return(numeric(0))
}

#' @title Compute the fqZ
#' @description Implements [compute_fqZ] for class(pars$xde)=mosy
#' @inheritParams compute_fqZ
#' @return a [vector]
#' @export
compute_fqZ.mosy <- function(deout, pars) {
  return(numeric(0))
}

#' @title Compute the fqZ
#' @description Implements [compute_fqZ] for class(pars$xde)=mosy_dde
#' @inheritParams compute_fqZ
#' @return a [vector]
#' @export
compute_fqZ.mosy_dde <- function(deout, pars) {
  return(numeric(0))
}

#' @title Compute the fqZ
#' @description Implements [compute_fqZ] for class(pars$xde)=human
#' @inheritParams compute_fqZ
#' @return a [vector]
#' @export
compute_fqZ.human<- function(deout, pars) {
  compute_fqZ_default(deout, pars)
}

#' @title Compute the fqZ
#' @description Implements [compute_fqZ] for class(pars$xde)=cohort
#' @inheritParams compute_fqZ
#' @return a [vector]
#' @export
compute_fqZ.cohort<- function(deout, pars) {
  return(numeric(0))
}

#' @title Compute the NI
#' @description Using the output of [deSolve::ode] or [deSolve::dede],
#' compute the NI for each stratum
#' @param deout a matrix, the output of deSolve
#' @param pars a [list]
#' @return [numeric] containing the NI
#' @export
compute_NI <- function(deout, pars) {
  UseMethod("compute_NI", pars$xde)
}

#' @title Compute the NI
#' @description Using the output of [deSolve::ode] or [deSolve::dede],
#' compute the NI for each stratum
#' @param deout a matrix, the output of deSolve
#' @param pars a [list]
#' @return [numeric] containing the NI
#' @export
compute_NI_default <- function(deout, pars) {
  ix = 1:length(deout[,1])
  NI = sapply(ix, compute_NI_i, deout=deout, pars=pars)
  return(t(NI))
}

#' @title Compute NI once
#' @description Using the output of [deSolve::ode] or [deSolve::dede],
#' compute the NI at a point in time
#' @param i an [integer]
#' @param deout a matrix, the output of deSolve
#' @param pars a [list]
#' @return a [list] containing the NI
#' @export
compute_NI_i <- function(i, deout, pars) {
  t = deout[i,1]
  y = deout[i,-1]
  X <- F_X(t, y, pars)
  H = F_H(t, y, pars)
  NI = X/H
  return(NI)
}

#' @title Compute the NI
#' @description Implements [compute_NI] for class(pars$xde)=ode
#' @inheritParams compute_NI
#' @return a [vector]
#' @export
compute_NI.ode <- function(deout, pars) {
  compute_NI_default(deout, pars)
}

#' @title Compute the NI
#' @description Implements [compute_NI] for class(pars$xde)=dde
#' @inheritParams compute_NI
#' @return a [vector]
#' @export
compute_NI.dde <- function(deout, pars) {
  compute_NI_default(deout, pars)
}

#' @title Compute the NI
#' @description Implements [compute_NI] for class(pars$xde)=aqua
#' @inheritParams compute_NI
#' @return a [vector]
#' @export
compute_NI.aqua <- function(deout, pars) {
  return(numeric(0))
}

#' @title Compute the NI
#' @description Implements [compute_NI] for class(pars$xde)=aqua_dde
#' @inheritParams compute_NI
#' @return a [vector]
#' @export
compute_NI.aqua_dde <- function(deout, pars) {
  return(numeric(0))
}

#' @title Compute the NI
#' @description Implements [compute_NI] for class(pars$xde)=mosy
#' @inheritParams compute_NI
#' @return a [vector]
#' @export
compute_NI.mosy <- function(deout, pars) {
  return(numeric(0))
}

#' @title Compute the NI
#' @description Implements [compute_NI] for class(pars$xde)=mosy_dde
#' @inheritParams compute_NI
#' @return a [vector]
#' @export
compute_NI.mosy_dde <- function(deout, pars) {
  return(numeric(0))
}

#' @title Compute the NI
#' @description Implements [compute_NI] for class(pars$xde)=human
#' @inheritParams compute_NI
#' @return a [vector]
#' @export
compute_NI.human<- function(deout, pars) {
  compute_NI_default(deout, pars)
}

#' @title Compute the NI
#' @description Implements [compute_NI] for class(pars$xde)=cohort
#' @inheritParams compute_NI
#' @return a [vector]
#' @export
compute_NI.cohort<- function(deout, pars) {
  compute_NI_default(deout, pars)
}

#' @title Compute kappa
#' @description Using the output of [deSolve::ode] or [deSolve::dede],
#' compute the kappa for each stratum
#' @param deout a matrix, the output of deSolve
#' @param pars a [list]
#' @return [numeric] containing the kappa
#' @export
compute_kappa <- function(deout, pars) {
  UseMethod("compute_kappa", pars$xde)
}

#' @title Compute the kappa
#' @description Using the output of [deSolve::ode] or [deSolve::dede],
#' compute the kappa for each stratum
#' @param deout a matrix, the output of deSolve
#' @param pars a [list]
#' @return [numeric] containing the kappa
#' @export
compute_kappa_default <- function(deout, pars) {
  ix = 1:length(deout[,1])
  kappa = sapply(ix, compute_kappa_i, deout=deout, pars=pars)
  return(t(kappa))
}

#' @title Compute kappa once
#' @description Using the output of [deSolve::ode] or [deSolve::dede],
#' compute the kappa at a point in time
#' @param i an [integer]
#' @param deout a matrix, the output of deSolve
#' @param pars a [list]
#' @return a [list] containing the kappa
#' @export
compute_kappa_i <- function(i, deout, pars) {
  t = deout[i,1]
  y = deout[i,-1]
  beta <- F_beta(t, y, pars)
  kappa = F_kappa(t, y, pars, beta)
  return(kappa)
}

#' @title Compute the kappa
#' @description Implements [compute_kappa] for class(pars$xde)=dde
#' @inheritParams compute_kappa
#' @return a [vector]
#' @export
compute_kappa.dde <- function(deout, pars) {
  compute_kappa_default(deout, pars)
}

#' @title Compute the kappa
#' @description Implements [compute_kappa] for class(pars$xde)=ode
#' @inheritParams compute_kappa
#' @return a [vector]
#' @export
compute_kappa.ode <- function(deout, pars) {
  compute_kappa_default(deout, pars)
}

#' @title Compute the kappa
#' @description Implements [compute_kappa] for class(pars$xde)=aqua
#' @inheritParams compute_kappa
#' @return a [vector]
#' @export
compute_kappa.aqua <- function(deout, pars) {
  return(numeric(0))
}

#' @title Compute the kappa
#' @description Implements [compute_kappa] for class(pars$xde)=aqua_dde
#' @inheritParams compute_kappa
#' @return a [vector]
#' @export
compute_kappa.aqua_dde <- function(deout, pars) {
  return(numeric(0))
}

#' @title Compute the kappa
#' @description Implements [compute_kappa] for class(pars$xde)=mosy
#' @inheritParams compute_kappa
#' @return a [vector]
#' @export
compute_kappa.mosy <- function(deout, pars) {
  return(numeric(0))
}

#' @title Compute the kappa
#' @description Implements [compute_kappa] for class(pars$xde)=mosy_dde
#' @inheritParams compute_kappa
#' @return a [vector]
#' @export
compute_kappa.mosy_dde <- function(deout, pars) {
  return(numeric(0))
}

#' @title Compute the kappa
#' @description Implements [compute_kappa] for class(pars$xde)=human
#' @inheritParams compute_kappa
#' @return a [vector]
#' @export
compute_kappa.human<- function(deout, pars) {
  compute_kappa_default(deout, pars)
}

#' @title Compute the kappa
#' @description Implements [compute_kappa] for class(pars$xde)=cohort
#' @inheritParams compute_kappa
#' @return a [vector]
#' @export
compute_kappa.cohort<- function(deout, pars) {
  return(numeric(0))
}
