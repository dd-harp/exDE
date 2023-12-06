# specialized methods for a human trace model

#' @title Size of effective infectious human population
#' @description Implements [F_X] for the trace model.
#' @inheritParams F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_X.trace <- function(t, y, pars) {
  with(pars$Xpar, kappa)
}

#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_b] for the trace model.
#' @inheritParams F_b
#' @return a [numeric] vector of length `nStrata`
#' @export
F_b.trace <- function(y, pars) {
  numeric(0)
}

#' @title Derivatives for human population
#' @description Implements [dXdt] for the trace model.
#' @inheritParams dXdt
#' @return a [numeric] vector
#' @export
dXdt.trace <- function(t, y, pars, FoI) {
  numeric(0)
}


#' @title Setup Xpar.trace
#' @description Implements [setup_X] for the SIS model
#' @inheritParams setup_X
#' @return a [list] vector
#' @export
setup_X.trace = function(pars, Xname, Xopts=list()){

  pars$Xname = "trace"
  pars = make_Xpar_trace(pars, Xopts)
  pars$Xinits = numeric(0)

  return(pars)
}

#' @title Make parameters for human null model
#' @param pars a [list]
#' @param Xopts a [list] that could overwrite defaults
#' @param kappa value
#' @param Kf a function
#' @return a [list]
#' @export
make_Xpar_trace = function(pars, Xopts=list(),
                         kappa = 0.1,
                         Kf = NULL){
  with(Xopts,{
    Xpar = list()
    class(Xpar) <- "trace"

    Xpar$kappa = checkIt(kappa, pars$nPatches)
    if(is.null(Kf)) Kf = function(t, y, pars){return(1)}
    Xpar$Kf = Kf

    pars$Xpar = Xpar
    pars$Hpar$wts_f = rep(1, pars$nPatches)
    pars$Hpar$H = rep(1, pars$nPatches)
    pars$Hpar$residence = c(1:pars$nPatches)
    pars$Hpar$TaR = diag(1, pars$nPatches)
    return(pars)
  })}

#' @title Add indices for human population to parameter list
#' @description Implements [make_indices_X] for the trace model.
#' @inheritParams make_indices_X
#' @return none
#' @importFrom utils tail
#' @export
make_indices_X.trace <- function(pars) {
  pars$X_ix <- integer(0)
  return(pars)
}

#' @title Parse the output of deSolve and return variables for the trace model
#' @description Implements [parse_deout_X] for the trace model
#' @inheritParams parse_deout_X
#' @return none
#' @export
parse_deout_X.trace <- function(deout, pars) {
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
  Xpar$kappa <- kappa
  pars$Xpar <- Xpar
  return(pars)
}

#' @title Make inits for trace human model
#' @param pars a [list]
#' @return none
#' @export
make_inits_X_trace <- function(pars) {
  pars$Xinits <- numeric(0)
  return(pars)
}

#' @title Update inits for the trace human model from a vector of states
#' @param pars a [list]
#' @param y0 a vector of initial values
#' @return none
#' @export
update_inits_X.trace <- function(pars, y0) {
  return(pars)
}

#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Xpar`.
#' @param pars a [list]
#' @return none
#' @export
get_inits_X.trace <- function(pars){
  numeric(0)
}
