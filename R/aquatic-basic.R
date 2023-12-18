# specialized methods for the aquatic mosquito basic competition model

#' @title Reset aquatic parameters to baseline
#' @description Implements [LBionomics] for the RM model
#' @inheritParams LBionomics
#' @return a named [list]
#' @export
LBionomics.basic <- function(t, y, pars) {
  pars$Lpar$psi <- pars$Lpar$psi0
  pars$Lpar$phi <- pars$Lpar$phi0
  pars$Lpar$theta <- pars$Lpar$theta0
  return(pars)
}


#' @title Number of newly emerging adults from each larval habitat
#' @description Implements [F_alpha] for the basic competition model.
#' @inheritParams F_alpha
#' @return a [numeric] vector of length `nHabitats`
#' @export
F_alpha.basic <- function(t, y, pars) {
  L <- y[pars$ix$L$L_ix]
  with(pars$Lpar,{
    psi*L
  })
}

#' @title Derivatives for aquatic stage mosquitoes
#' @description Implements [dLdt] for the basic competition model.
#' @inheritParams dLdt
#' @return a [numeric] vector
#' @export
dLdt.basic <- function(t, y, pars, eta) {
  L <- y[pars$ix$L$L_ix]
  with(pars$Lpar, {
    dL = eta - (psi + phi + (theta*L))*L
    return(dL)
  })
}

#' @title Setup Lpar.basic
#' @description Implements [setup_L] for the basic model
#' @inheritParams setup_L
#' @return a [list] vector
#' @export
setup_L.basic = function(pars, Lname,
                            membership=1, searchQ=1,
                            Lopts=list()){

  pars$nHabitats=length(membership)
  pars$Lname = "basic"
  pars$membership = membership
  pars$searchQ = checkIt(searchQ, length(membership), F)
  pars$calN = make_calN(pars$nPatches, membership)
  pars$calU = make_calU(pars$calN, pars$searchQ)

  with(Lopts,{

    pars <- make_Lpar_basic(pars, Lopts)
    pars <- LBionomics.basic(0, 0, pars)
    pars <- make_Linits_basic(pars, Lopts)

    return(pars)
 })}

#' @title Make parameters for basic competition aquatic mosquito model
#' @param pars a [list]
#' @param Lopts a [list] that overwrites default values
#' @param psi maturation rates for each aquatic habitat
#' @param phi density-independent mortality rates for each aquatic habitat
#' @param theta density-dependent mortality terms for each aquatic habitat
#' @return a [list] with Lpar added
#' @export
make_Lpar_basic = function(pars, Lopts=list(), psi=1/8, phi=1/8, theta=1/100){with(Lopts,{
  Lpar = list()
  class(Lpar) <- "basic"
  Lpar$psi0 = checkIt(psi, pars$nHabitats)
  Lpar$phi0 = checkIt(phi, pars$nHabitats)
  Lpar$theta0 = checkIt(theta, pars$nHabitats)

  pars$Lpar = Lpar
  return(pars)
})}

#' @title Make inits for basic competition aquatic mosquito model
#' @param pars a [list]
#' @param Lopts a [list] that overwrites default values
#' @param L0 initial conditions
#' @return a [list] with Linits added
#' @export
make_Linits_basic = function(pars, Lopts=list(), L0 = 1){with(Lopts,{
  inits = list()
  inits$L0 = checkIt(L0, pars$nHabitats)
  pars$Linits = inits
  return(pars)
})}

#' @title Add indices for aquatic stage mosquitoes to parameter list
#' @description Implements [make_indices_L] for basic competition model.
#' @inheritParams make_indices_L
#' @return none
#' @importFrom utils tail
#' @export
make_indices_L.basic <- function(pars) {
  pars$ix$L$L_ix <- seq(from = pars$max_ix+1, length.out = pars$nHabitats)
  pars$max_ix <- tail(pars$ix$L$L_ix, 1)
  return(pars)
}

#' @title Parse the variable names for the basic model
#' @description Implements [parse_deout_L] for basic competition model.
#' @inheritParams parse_deout_L
#' @return [list]
#' @export
parse_deout_L.basic <- function(deout, pars) {
  time = deout[,1]
  L = deout[,pars$ix$L$L_ix+1]
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
  Lpar$psi0 <- psi
  Lpar$phi0 <- phi
  Lpar$theta0 <- theta
  pars$Lpar <- Lpar
  pars <- LBionomics.basic(0, 0, pars)
  return(pars)
}

#' @title Make inits for basic competition aquatic mosquito model
#' @param pars a [list]
#' @param L0 initial conditions
#' @return a [list] with Linits added
#' @export
make_inits_L_basic <- function(pars, L0){
  stopifnot(is.numeric(L0))
  pars$Linits <- list(L0=L0)
  return(pars)
}

#' @title Update inits for the basic aquatic mosquito competition model
#' @param pars a [list]
#' @param y0 a vector of initial values
#' @return none
#' @export
update_inits_L.basic <- function(pars, y0) {
  L0 = y0[pars$ix$L$L_ix]
  pars = make_inits_L_basic(pars, L0)
  return(pars)
}

#' @title Return initial values as a vector
#' @description Implements [get_inits_L] for the GeRM model.
#' @inheritParams get_inits_L
#' @return none
#' @export
get_inits_L.basic <- function(pars){
  pars$Linits$L0
}

