# specialized methods for the aquatic mosquito basic competition model

#' @title Reset aquatic parameters to baseline
#' @description Implements [LBionomics] for the RM model
#' @inheritParams LBionomics
#' @return a named [list]
#' @export
LBionomics.basic <- function(t, y, pars, s) {with(pars$Lpar[[s]],{
  pars$Lpar[[s]]$psi <- psi0
  pars$Lpar[[s]]$phi <- phi0
  pars$Lpar[[s]]$theta <- theta0

  return(pars)
})}


#' @title Number of newly emerging adults from each larval habitat
#' @description Implements [F_alpha] for the basic competition model.
#' @inheritParams F_alpha
#' @return a [numeric] vector of length `nHabitats`
#' @export
F_alpha.basic <- function(t, y, pars, s) {
  L <- y[pars$ix$L[[s]]$L_ix]
  with(pars$Lpar[[s]],{
    return(psi*L)
  })
}

#' @title Derivatives for aquatic stage mosquitoes
#' @description Implements [dLdt] for the basic competition model.
#' @inheritParams dLdt
#' @return a [numeric] vector
#' @export
dLdt.basic <- function(t, y, pars, s) {
  eta <- pars$eta[[s]]
  with(pars$ix$L[[s]],{
    L <- y[L_ix]
    with(pars$Lpar[[s]], {
      dL = eta - (psi + phi + (theta*L))*L
      return(dL)
    })
  })
}

#' @title Setup Lpar for the basic model
#' @description Implements [setup_Lpar] for the basic model
#' @inheritParams setup_Lpar
#' @return a [list] vector
#' @export
setup_Lpar.basic = function(Lname, pars, s, Lopts=list()){
  pars$Lpar[[s]] = make_Lpar_basic(pars$nHabitats, Lopts)
  return(pars)
}

#' @title Setup the basic model
#' @description Implements [setup_Linits] for the basic model
#' @inheritParams setup_Linits
#' @return a [list]
#' @export
setup_Linits.basic = function(pars, s, Lopts=list()){
  pars$Linits[[s]] = make_Linits_basic(pars$nHabitats, Lopts)
  return(pars)
}

#' @title Make parameters for basic competition aquatic mosquito model
#' @param nHabitats the number of habitats in the model
#' @param Lopts a [list] that overwrites default values
#' @param psi maturation rates for each aquatic habitat
#' @param phi density-independent mortality rates for each aquatic habitat
#' @param theta density-dependent mortality terms for each aquatic habitat
#' @return a [list] with Lpar added
#' @export
make_Lpar_basic = function(nHabitats, Lopts=list(), psi=1/8, phi=1/8, theta=1/100){with(Lopts,{
  Lpar = list()
  class(Lpar) <- "basic"
  Lpar$psi0 = checkIt(psi, nHabitats)
  Lpar$phi0 = checkIt(phi, nHabitats)
  Lpar$theta0 = checkIt(theta, nHabitats)

  return(Lpar)
})}

#' @title Make inits for basic competition aquatic mosquito model
#' @param nHabitats the number of habitats in the model
#' @param Lopts a [list] that overwrites default values
#' @param L0 initial conditions
#' @return a [list] with Linits added
#' @export
make_Linits_basic = function(nHabitats, Lopts=list(), L0=1){with(Lopts,{
  L0 = checkIt(L0, nHabitats)
  return(list(L=L0))
})}

#' @title Add indices for aquatic stage mosquitoes to parameter list
#' @description Implements [make_indices_L] for basic competition model.
#' @inheritParams make_indices_L
#' @return none
#' @importFrom utils tail
#' @export
make_indices_L.basic <- function(pars, s) {with(pars,{

  L_ix <- seq(from = max_ix+1, length.out=nHabitats)
  max_ix <- tail(L_ix, 1)

  pars$max_ix = max_ix
  pars$ix$L[[s]] = list(L_ix=L_ix)
  return(pars)
})}

#' @title Parse the variable names for the basic model
#' @description Implements [parse_deout_L] for basic competition model.
#' @inheritParams parse_deout_L
#' @return [list]
#' @export
parse_deout_L.basic <- function(deout, pars, s) {
  time = deout[,1]
  L = deout[,pars$ix$L[[s]]$L_ix+1]
  return(list(time=time, L=L))
}


#' @title Make parameters for basic competition aquatic mosquito model
#' @param pars a [list]
#' @param psi maturation rates for each aquatic habitat
#' @param phi density-independent mortality rates for each aquatic habitat
#' @param theta density-dependent mortality terms for each aquatic habitat
#' @return a [list] with Lpar added
#' @export
make_parameters_L_basic <- function(pars, psi, phi, theta) {
  stopifnot(is.numeric(psi), is.numeric(phi), is.numeric(theta))
  Lpar <- list()
  class(Lpar) <- 'basic'
  Lpar$psi <- psi
  Lpar$phi <- phi
  Lpar$theta <- theta

  Lpar$psi0 <- psi
  Lpar$phi0 <- phi
  Lpar$theta0 <- theta

  pars$Lpar[[1]] <- Lpar
  return(pars)
}

#' @title Make inits for basic competition aquatic mosquito model
#' @param pars a [list]
#' @param L0 initial conditions
#' @return a [list] with Linits added
#' @export
make_inits_L_basic <- function(pars, L0){
  stopifnot(is.numeric(L0))
  pars$Linits[[1]] <- list(L=L0)
  return(pars)
}

#' @title Update inits for the basic aquatic mosquito competition model
#' @inheritParams update_inits_L
#' @return none
#' @export
update_inits_L.basic <- function(pars, y0, s) {
  L = y0[pars$ix$L[[s]]$L_ix]
  pars = make_Linits_basic(pars$nHabitats, L0=L)
  return(pars)
}

#' @title Return initial values as a vector
#' @description Implements [get_inits_L] for the GeRM model.
#' @inheritParams get_inits_L
#' @return none
#' @export
get_inits_L.basic <- function(pars, s){
  pars$Linits[[s]]$L
}

