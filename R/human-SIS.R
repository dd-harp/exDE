# specialized methods for the human SIS model

#' @title Size of effective infectious human population
#' @description Implements [F_X] for the SIS model.
#' @inheritParams F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_X.SIS <- function(t, y, pars) {
  with(pars$Xpar, y[X_ix]*c)
}

#' @title Derivatives for human population
#' @description Implements [dXdt] for the SIS model, no demography.
#' @inheritParams dXdt
#' @return a [numeric] vector
#' @export
dXdt.SISdX <- function(t, y, pars, EIR) {
  with(pars$Xpar, {
    X <- y[X_ix]
    H <- F_H(t, y, pars)
    foi = F_foi(b*EIR, pars)
    dX <- foi*(H - X) - r*X
    return(c(dX))
  })
}


#' @title Derivatives for human population
#' @description Implements [dXdt] for the SIS model with demography.
#' @inheritParams dXdt
#' @return a [numeric] vector
#' @export
dXdt.SISdXdH <- function(t, y, pars, EIR) {

  with(pars$Xpar, {
    H <- F_H(t, y, pars)
    X <- y[X_ix]
    foi = F_foi(b*EIR, pars)
    dX <- foi*(H - X) - r*X + dHdt(t, X, pars)
    dH <- Births(t, H, pars) + dHdt(t, H, pars)

    return(c(dX, dH))
  })
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
#' @param pars an [list]
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
#' @param pars an [environment]
#' @param X0 size of infected population in each strata
#' @return none
#' @export
make_inits_X_SIS <- function(pars, X0) {
  stopifnot(is.numeric(X0))
  pars$Xinits <- list(X0=X0)
  return(pars)
}

#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Xpar`.
#' @param pars an [environment]
#' @return none
#' @export
get_inits_X.SIS <- function(pars){
  pars$Xinits$X0
}
