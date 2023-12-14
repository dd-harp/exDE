# specialized methods for the human SIP model

#' @title Size of effective infectious human population
#' @description Implements [F_X] for the SIP model.
#' @inheritParams F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_X.SIP <- function(t, y, pars) {
  with(pars$Xpar, y[X_ix] * c)
}

#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_pr] for the SIP model.
#' @inheritParams F_pr
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pr.SIP <- function(varslist, pars) {
  pr = with(varslist$XH, X/H)
  return(pr)
}

#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_b] for the SIP model.
#' @inheritParams F_b
#' @return a [numeric] vector of length `nStrata`
#' @export
F_b.SIP <- function(y, pars) {
  with(pars$Xpar, b)
}

#' @title Derivatives for human population
#' @description Implements [dXdt] for the SIP model.
#' @inheritParams dXdt
#' @return a [numeric] vector
#' @export
dXdt.SIPdX <- function(t, y, pars, FoI) {

  with(pars$Xpar, {
    X <- y[X_ix]
    P <- y[P_ix]
    H <- F_H(t, y, pars)

    dX <- (1-rho)*FoI*(H - X - P) - (r+xi)*X
    dP <- rho*FoI*(H - X - P) + xi*(H-P) - eta*P

    return(c(dX, dP))
  })
}


#' @title Compute the HTC for the SIP model
#' @description Implements [HTC] for the SIP model with demography.
#' @inheritParams HTC
#' @return a [numeric] vector
#' @export
HTC.SIP <- function(pars) {
  with(pars$Xpar,
       return((1-rho)*b/(r+xi)*xi/(eta+xi))
  )
}

#' @title Derivatives for human population
#' @description Implements [dXdt] for the SIP model.
#' @inheritParams dXdt
#' @return a [numeric] vector
#' @export
dXdt.SIPdXdH <- function(t, y, pars, FoI) {

  with(pars$Xpar, {

    X <- y[X_ix]
    P <- y[P_ix]
    H <- F_H(t, y, pars)

    dX <- (1-rho)*FoI*(H - X - P) - (r+xi)*X + dHdt(t, X, pars)
    dP <- rho*FoI*(H - X - P) + xi*(H-P) - eta*P + dHdt(t, P, pars)
    dH <- Births(t, H, pars) + dHdt(t, H, pars)

    return(c(dX, dP, dH))
  })
}

#' @title Setup Xpar.SIP
#' @description Implements [setup_X] for the SIP model
#' @inheritParams setup_X
#' @return a [list] vector
#' @export
setup_X.SIP = function(pars, Xname, Xopts=list()){

  pars$Xname = "SIP"
    pars = make_Xpar_SIP(pars, Xopts)
  pars = make_Xinits_SIP(pars, Xopts)

  return(pars)
}

#' @title Make parameters for SIP human model, with defaults
#' @param pars a [list]
#' @param Xopts a [list] that could overwrite defaults
#' @param b transmission probability (efficiency) from mosquito to human
#' @param c transmission probability (efficiency) from human to mosquito
#' @param r recovery rate
#' @param rho probability of successful treatment upon infection
#' @param eta prophylaxis waning rate
#' @param xi background treatment rate
#' @return a [list]
#' @export
make_Xpar_SIP = function(pars, Xopts=list(),
                         b=0.55, r=1/180, c=0.15,
                         rho=.1, eta=1/25, xi=1/365){
  with(Xopts,{
    Xpar = list()
    class(Xpar) <- c("SIP", "SIPdX")

    Xpar$b = checkIt(b, pars$nStrata)
    Xpar$c = checkIt(c, pars$nStrata)
    Xpar$r = checkIt(r, pars$nStrata)
    Xpar$rho = checkIt(rho, pars$nStrata)
    Xpar$eta = checkIt(eta, pars$nStrata)
    Xpar$xi = checkIt(xi, pars$nStrata)

    pars$Xpar = Xpar
    return(pars)
})}


#' @title Make initial values for the SIP human model, with defaults
#' @param pars a [list]
#' @param Xopts a [list] that could overwrite defaults
#' @param X0 the initial values of the parameter X
#' @param P0 the initial values of the parameter P
#' @return a [list]
#' @export
make_Xinits_SIP = function(pars, Xopts = list(),
                           X0=1, P0=0){with(Xopts,{
  inits = list()
  inits$X0 = checkIt(X0, pars$nStrata)
  inits$P0 = checkIt(P0, pars$nStrata)
  pars$Xinits = inits
  return(pars)
})}

#' @title Parse the output of deSolve and return variables for the SIP model
#' @description Implements [parse_deout_X] for the SIP model
#' @inheritParams parse_deout_X
#' @return none
#' @export
parse_deout_X.SIP <- function(deout, pars) {
  time = deout[,1]
  Hlist <- parse_deout_H(deout, pars)
  with(Hlist,{
    X = deout[,pars$Xpar$X_ix+1]
    P = deout[,pars$Xpar$P_ix+1]
  return(list(time=time, X=X,P=P,H=H))
})}

#' @title Add indices for human population to parameter list
#' @description Implements [make_indices_X] for the SIP model.
#' @inheritParams make_indices_X
#' @return none
#' @importFrom utils tail
#' @export
make_indices_X.SIP <- function(pars) {

  pars$Xpar$X_ix <- seq(from = pars$max_ix+1, length.out = pars$nStrata)
  pars$max_ix <- tail(pars$Xpar$X_ix, 1)

  pars$Xpar$P_ix <- seq(from = pars$max_ix+1, length.out = pars$nStrata)
  pars$max_ix <- tail(pars$Xpar$P_ix, 1)

  return(pars)
}

#' @title Make parameters for SIP human model
#' @param pars a [list]
#' @param b transmission probability (efficiency) from mosquito to human
#' @param c transmission probability (efficiency) from human to mosquito
#' @param r recovery rate
#' @param rho probability of successful treatment upon infection
#' @param eta prophylaxis waning rate
#' @param xi background treatment rate
#' @return none
#' @export
make_parameters_X_SIP <- function(pars, b, c, r, rho, eta, xi){
  stopifnot(is.numeric(b), is.numeric(c), is.numeric(r), is.numeric(rho), is.numeric(eta), is.numeric(xi))
  Xpar <- list()
  class(Xpar) <- c('SIP', 'SIPdX')
  Xpar$b <- b
  Xpar$c <- c
  Xpar$r <- r
  Xpar$rho <- rho
  Xpar$eta <- eta
  Xpar$xi <- xi
  pars$Xpar <- Xpar
  return(pars)
}

#' @title Make inits for SIP human model
#' @param pars a [list]
#' @param X0 size of infected population in each strata
#' @param P0 size of population protected by prophylaxis in each strata
#' @return none
#' @export
make_inits_X_SIP <- function(pars, X0, P0) {
  stopifnot(is.numeric(X0), is.numeric(P0))
  pars$Xinits = list(X0 = X0, P0 = P0)
  return(pars)
}

#' @title Update inits for the SIP human model from a vector of states
#' @param pars a [list]
#' @param y0 a vector of initial values
#' @return none
#' @export
update_inits_X.SIP <- function(pars, y0) {
  X0 = y0[pars$Xpar$X_ix]
  P0 = y0[pars$Xpar$P_ix]
  pars = make_inits_X_SIP(pars, X0, P0)
  return(pars)
}


#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Xpar`.
#' @param pars a [list]
#' @return none
#' @export
get_inits_X.SIP <- function(pars){with(pars$Xinits,{
  c(X0, P0)
})}


#' Plot the density of infected individuals for the SIP model
#'
#' @inheritParams xde_plot_X
#' @export
xde_plot_X.SIP = function(pars, clrs=c("black", "darkgreen"), llty=1, stable=FALSE, add_axes=TRUE){
  vars=with(pars$outputs,if(stable==TRUE){stable_orbits}else{orbits})

  if(add_axes==TRUE)
    with(vars$XH,
         plot(time, 0*time, type = "n", ylim = c(0, max(H)),
              ylab = "# Infected", xlab = "Time"))

  xde_lines_X(vars$XH, pars, clrs, llty)
}


#' Add lines for the density of infected individuals for the SIP model
#'
#' @inheritParams xde_lines_X
#'
#' @export
xde_lines_X.SIP = function(XH, pars, clrs=c("black", "darkgreen"), llty=1){
  with(XH,{
    if(pars$nStrata==1) {
      lines(time, X, col=clrs[1], lty = llty[1])
      lines(time, P, col=clrs[2], lty = llty[1])
    }
    if(pars$nStrata>1){
      if (length(clrs)==1) clrs=matrix(clrs, 2, pars$nStrata)
      if (length(llty)==1) llty=rep(llty, pars$nStrata)
      for(i in 1:pars$nStrata){
        lines(time, X[,i], col=clrs[1,i], lty = llty[i])
        lines(time, P[,i], col=clrs[2,i], lty = llty[i])
      }
    }
  })}

