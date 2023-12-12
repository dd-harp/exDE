# specialized methods for the human SIS model

#' @title Size of effective infectious human population
#' @description Implements [F_X] for the SIS model.
#' @inheritParams F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_X.SIS <- function(t, y, pars) {
  with(pars$Xpar, y[X_ix]*c)
}

#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_pr] for the SIS model.
#' @inheritParams F_pr
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pr.SIS <- function(varslist, pars) {
  pr = with(varslist$XH, X/H)
  return(pr)
}


#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_b] for the SIS model.
#' @inheritParams F_b
#' @return a [numeric] vector of length `nStrata`
#' @export
F_b.SIS <- function(y, pars) {
  with(pars$Xpar, b)
}

#' @title Derivatives for human population
#' @description Implements [dXdt] for the SIS model, no demography.
#' @inheritParams dXdt
#' @return a [numeric] vector
#' @export
dXdt.SISdX <- function(t, y, pars, FoI) {
  with(pars$Xpar, {

    X <- y[X_ix]
    H <- F_H(t, y, pars)

    dX <- FoI*(H - X) - r*X

    return(c(dX))
  })
}

#' @title Derivatives for human population
#' @description Implements [dXdt] for the SIS model with demography.
#' @inheritParams dXdt
#' @return a [numeric] vector
#' @export
dXdt.SISdXdH <- function(t, y, pars, FoI) {
  with(pars$Xpar, {

    H <- F_H(t, y, pars)
    X <- y[X_ix]

    dX <- FoI*(H - X) - r*X + dHdt(t, X, pars)
    dH <- Births(t, H, pars) + dHdt(t, H, pars)

    return(c(dX, dH))
  })
}

#' @title Setup Xpar.SIS
#' @description Implements [setup_X] for the SIS model
#' @inheritParams setup_X
#' @return a [list] vector
#' @export
setup_X.SIS = function(pars, Xname, Xopts=list()){

  pars$Xname = "SIS"
  pars = make_Xpar_SIS(pars, Xopts)
  pars = make_Xinits_SIS(pars, Xopts)

  return(pars)
}

#' @title Make parameters for SIS human model, with defaults
#' @param pars a [list]
#' @param Xopts a [list] that could overwrite defaults
#' @param b transmission probability (efficiency) from mosquito to human
#' @param c transmission probability (efficiency) from human to mosquito
#' @param r recovery rate
#' @return a [list]
#' @export
make_Xpar_SIS = function(pars, Xopts=list(),
                         b=0.55, r=1/180, c=0.15){
  with(Xopts,{
    Xpar = list()
    class(Xpar) <- c("SIS", "SISdX")

    Xpar$b = checkIt(b, pars$nStrata)
    Xpar$c = checkIt(c, pars$nStrata)
    Xpar$r = checkIt(r, pars$nStrata)

    pars$Xpar = Xpar
    return(pars)
  })}

#' @title Make initial values for the SIS human model, with defaults
#' @param pars a [list]
#' @param Xopts a [list] to overwrite defaults
#' @param X0 the initial values of the parameter X
#' @return a [list]
#' @export
make_Xinits_SIS = function(pars, Xopts = list(), X0=1){with(Xopts,{
  inits = list()
  inits$X0 = checkIt(X0, pars$nStrata)
  pars$Xinits = inits
  return(pars)
})}

#' @title Parse the output of deSolve and return variables for the SIS model
#' @description Implements [parse_deout_X] for the SIS model
#' @inheritParams parse_deout_X
#' @return none
#' @export
parse_deout_X.SIS <- function(deout, pars) {
  time = deout[,1]
  Hlist <- parse_deout_H(deout, pars)
  with(Hlist,{
    X = deout[,pars$Xpar$X_ix+1]
    return(list(time=time, X=X, H=H))
})}

#' @title Compute the HTC for the SIS model
#' @description Implements [HTC] for the SIS model with demography.
#' @inheritParams HTC
#' @return a [numeric] vector
#' @export
HTC.SIS <- function(pars) {
  with(pars$Xpar,
    return(c/r)
  )
}

#' @title Add indices for human population to parameter list
#' @description Implements [make_indices_X] for the SIS model.
#' @inheritParams make_indices_X
#' @return none
#' @importFrom utils tail
#' @export
make_indices_X.SIS <- function(pars) {
  pars$Xpar$X_ix <- seq(from = pars$max_ix+1, length.out = pars$nStrata)
  pars$max_ix <- tail(pars$Xpar$X_ix, 1)
  return(pars)
}

#' @title Make parameters for SIS human model
#' @param pars a [list]
#' @param b transmission probability (efficiency) from mosquito to human
#' @param c transmission probability (efficiency) from human to mosquito
#' @param r recovery rate
#' @return a [list]
#' @export
make_parameters_X_SIS <- function(pars, b, c, r) {
  stopifnot(is.numeric(b), is.numeric(c), is.numeric(r))
  Xpar <- list()
  class(Xpar) <- c('SIS', 'SISdX')
  Xpar$b <- b
  Xpar$c <- c
  Xpar$r <- r
  pars$Xpar <- Xpar
  return(pars)
}

#' @title Make inits for SIS human model
#' @param pars a [list]
#' @param X0 size of infected population in each strata
#' @return none
#' @export
make_inits_X_SIS <- function(pars, X0) {
  stopifnot(is.numeric(X0))
  pars$Xinits <- list(X0=X0)
  return(pars)
}

#' @title Update inits for the SIS human model from a vector of states
#' @param pars a [list]
#' @param y0 a vector of initial values
#' @return none
#' @export
update_inits_X.SIS <- function(pars, y0) {
  X0 = y0[pars$Xpar$X_ix]
  pars = make_inits_X_SIS(pars, X0)
  return(pars)
}


#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Xpar`.
#' @param pars a [list]
#' @return none
#' @export
get_inits_X.SIS <- function(pars){
  pars$Xinits$X0
}
