# specialized methods for a human trace model

#' @title Size of effective infectious human population
#' @description Implements [F_X] for the trace model.
#' @inheritParams F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_X.trace <- function(t, y, pars, i) {
  with(pars$Xpar[[i]],  pars$Hpar[[i]]$H*Kf(t))
}

#' @title Size of the human population
#' @description Implements [F_H] for the trace model.
#' @inheritParams F_H
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.trace <- function(t, y, pars, i) {
  with(pars$Hpar[[i]], H)
}

#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_pr] for the trace model.
#' @inheritParams F_pr
#' @return a [numeric] vector numeric(0)
#' @export
F_pr.trace <- function(varslist, pars, i) {
  return(numeric(0))
}

#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_b] for the trace model.
#' @inheritParams F_b
#' @return a [numeric] vector of length `nStrata`
#' @export
F_b.trace <- function(y, pars, i) {
  return(1)
}

#' @title Derivatives for human population
#' @description Implements [dXdt] for the trace model.
#' @inheritParams dXdt
#' @return a [numeric] vector
#' @export
dXdt.trace <- function(t, y, pars, i) {
  numeric(0)
}

#' @title Setup Xpar.trace
#' @description Implements [setup_Xpar] for the trace model
#' @inheritParams setup_Xpar
#' @return a [list] vector
#' @export
setup_Xpar.trace = function(Xname, pars, i, Xopts=list()){
  nStrata= pars$nPatches

  pars$Hpar[[i]]$nStrata = nStrata
  pars$Hpar[[i]]$H = checkIt(pars$Hpar[[i]]$H, nStrata)
  pars$BFpar$TimeSpent[[i]] = diag(1, nStrata)
  pars = make_TaR(0, pars, i, 1)
  pars$BFpar$searchWts[[i]][[1]] = checkIt(pars$BFpar$searchWts[[i]][[1]], nStrata)
  pars$BFpar$residence[[i]] = checkIt(pars$BFpar$residence[[i]], nStrata)
  pars$Xpar[[i]] = make_Xpar_trace(nStrata, Xopts)
  return(pars)
}

#' @title Setup Xinits.trace
#' @description Implements [setup_Xinits] for the trace model
#' @inheritParams setup_Xinits
#' @return a [list] vector
#' @export
setup_Xinits.trace = function(pars, i, Xopts=list()){
  pars$Xinits[[i]] = list()
  return(pars)
}

#' @title Make parameters for human null model
#' @param nStrata the number of population strata in the model
#' @param Xopts a [list] that could overwrite defaults
#' @param kappa value
#' @param Kf a function
#' @return a [list]
#' @export
make_Xpar_trace = function(nStrata, Xopts=list(),
                         kappa = 0.1,
                         Kf = NULL){
  with(Xopts,{
    Xpar = list()
    class(Xpar) <- "trace"
    kappa = checkIt(kappa, nStrata)
    if(is.null(Kf)) Kf = function(t){return(kappa)}
    Xpar$Kf = Kf
    return(Xpar)
})}

#' @title Add indices for human population to parameter list
#' @description Implements [make_indices_X] for the trace model.
#' @inheritParams make_indices_X
#' @return none
#' @importFrom utils tail
#' @export
make_indices_X.trace <- function(pars, i) {
  return(pars)
}

#' @title Parse the output of deSolve and return variables for the trace model
#' @description Implements [parse_deout_X] for the trace model
#' @inheritParams parse_deout_X
#' @return none
#' @export
parse_deout_X.trace <- function(deout, pars,i) {
  return(list())
}

#' @title Make parameters for trace human model
#' @param pars a [list]
#' @param kappa net infectiousness
#' @return a [list]
#' @export
make_parameters_X_trace <- function(pars, kappa) {
  Xpar <- list()
  class(Xpar) <- c('trace')
  Xpar$scale = checkIt(kappa, pars$nPatches)
  Xpar$Kf = function(t){return(1)}
  pars$Xpar[[1]] <- Xpar
  return(pars)
}

#' @title Make inits for trace human model
#' @param pars a [list]
#' @return none
#' @export
make_inits_X_trace <- function(pars) {
  pars$Xinits[[1]] <- numeric(0)
  return(pars)
}

#' @title Update inits for the trace human model from a vector of states
#' @inheritParams update_inits_X
#' @return none
#' @export
update_inits_X.trace <- function(pars, y0, i) {
  return(pars)
}

#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Xpar`.
#' @inheritParams get_inits_X
#' @return none
#' @export
get_inits_X.trace <- function(pars, i){
  numeric(0)
}
