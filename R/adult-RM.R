# specialized methods for the adult mosquito RM model

#' @title Reset bloodfeeding and mortality rates to baseline
#' @description Implements [MBionomics] for the RM model
#' @inheritParams MBionomics
#' @return a named [list]
#' @export
MBionomics.RM <- function(t, y, pars) {

  pars$MYZpar$f <- pars$MYZpar$f0
  pars$MYZpar$q <- pars$MYZpar$q0
  pars$MYZpar$g <- pars$MYZpar$g0
  pars$MYZpar$sigma <- pars$MYZpar$sigma0
  pars$MYZpar$nu <- pars$MYZpar$nu0

  return(pars)
}

#' @title Blood feeding rate of the infective mosquito population
#' @description Implements [F_fqZ] for the RM model.
#' @inheritParams F_fqZ
#' @return a [numeric] vector of length `nPatches`
#' @export
F_fqZ.RM <- function(t, y, pars) {
  with(pars$MYZpar, f*q)*y[pars$MYZpar$Z_ix]
}

#' @title Number of eggs laid by adult mosquitoes
#' @description Implements [F_eggs] for the RM model.
#' @inheritParams F_eggs
#' @return a [numeric] vector of length `nPatches`
#' @export
F_eggs.RM <- function(t, y, pars) {
  with(pars$MYZpar, {
    M <- y[M_ix]
    return(M*nu*eggsPerBatch)
  })
}

#' @title Derivatives for adult mosquitoes
#' @description Implements [dMYZdt] for the RM ODE model.
#' @inheritParams dMYZdt
#' @return a [numeric] vector
#' @export
dMYZdt.RM_ode <- function(t, y, pars, Lambda, kappa) {

  nPatches <- pars$nPatches

  with(pars$MYZpar,{

    M <- y[M_ix]
    P <- y[P_ix]
    Y <- y[Y_ix]
    Z <- y[Z_ix]

    Omega <- make_Omega(g, sigma, calK, nPatches)
    Upsilon <- expm(-Omega*eip)

    dMdt <- Lambda - (Omega %*% M)
    dPdt <- f*(M - P) - (Omega %*% P)
    dYdt <- f*q*kappa*(M - Y) - (Omega %*% Y)
    dZdt <- Upsilon %*% diag(f*q*kappa, nPatches) %*% (M - Y) - (Omega %*% Z)

    return(c(dMdt, dPdt, dYdt, dZdt))
  })
}

#' @title Derivatives for adult mosquitoes
#' @description Implements [dMYZdt] for the RM DDE model.
#' @inheritParams dMYZdt
#' @return a [numeric] vector
#' @importFrom deSolve lagvalue
#' @importFrom deSolve lagderiv
#' @export
dMYZdt.RM_dde <- function(t, y, pars, Lambda, kappa) {

  nPatches <- pars$nPatches
  eip <- pars$MYZpar$eip
  with(pars$MYZpar,{

    if (t < eip) {
      M_eip <- pars$MYZinits$M0
      Y_eip <- pars$MYZinits$Y0
      fqkappa_eip <- kappa*with(pars$MYZpar,f*q)
      g_eip <- g
      sigma_eip <- sigma
    } else {
      M_eip <- lagvalue(t = t - eip, nr = M_ix)
      Y_eip <- lagvalue(t = t - eip, nr = Y_ix)
      fqkappa_eip <- lagderiv(t = t-eip, nr = fqkappa_ix)
      g_eip <- lagderiv(t = t-eip, nr = g_ix)
      sigma_eip <- lagderiv(t = t-eip, nr = sigma_ix)
    }

    M <- y[M_ix]
    P <- y[P_ix]
    Y <- y[Y_ix]
    Z <- y[Z_ix]
    Upsilon <- matrix(data = y[Upsilon_ix], nrow = nPatches, ncol = nPatches)

    Omega <- make_Omega(g, sigma, calK, nPatches)
    Omega_eip <- make_Omega(g_eip, sigma_eip, calK, nPatches)

    dMdt <- Lambda - (Omega %*% M)
    dPdt <- f*(M - P) - (Omega %*% P)
    dYdt <- f*q*kappa*(M - Y) - (Omega %*% Y)
    dZdt <- Upsilon %*% (fqkappa_eip * (M_eip - Y_eip)) - (Omega %*% Z)
    dUdt <- as.vector((Omega_eip - Omega) %*% Upsilon)

    return(c(dMdt, dPdt, dYdt, dZdt, dUdt, f*q*kappa, g, sigma))
  })
}


#' @title Setup the RM model
#' @description Implements [setup_MYZ] for the RM model
#' @inheritParams setup_MYZ
#' @return a [list] vector
#' @export
setup_MYZ.RM = function(pars, MYZname, nPatches=1, MYZopts=list(), calK=diag(1)){

  pars$MYZname = "RM"
  pars$nPatches = nPatches

  pars = make_MYZpar_RM(pars, MYZopts, calK)
  pars = make_MYZinits_RM(pars, MYZopts)

  # extra setup for dde solving
  if(pars$MYZpar$xde == "dde"){
    pars$xde = "dde"
    class(pars$xde) = "dde"
    Omega <- with(pars$MYZpar, make_Omega(g, sigma, calK, nPatches))
    Upsilon <- expm::expm(-Omega*pars$MYZpar$eip)
    pars$MYZinits$Upsilon0 = as.vector(Upsilon)
    pars$MYZinits$dummy = rep(0, 3*pars$nPatches)
  }
  return(pars)
}


#' @title Make parameters for RM ODE adult mosquito model
#' @param pars a [list]
#' @param MYZopts a [list] of values that overwrites the defaults
#' @param calK mosquito dispersal matrix of dimensions `nPatches` by `nPatches`
#' @param solve_as is either `ode` to solve as an ode or `dde` to solve as a dde
#' @param g mosquito mortality rate
#' @param sigma emigration rate
#' @param f feeding rate
#' @param q human blood fraction
#' @param nu oviposition rate, per mosquito
#' @param eggsPerBatch eggs laid per oviposition
#' @param eip length of extrinsic incubation period
#' @return a [list]
#' @export
make_MYZpar_RM = function(pars, MYZopts=list(), calK,
                          solve_as = "dde",
                          g=1/12, sigma=1/8,
                          f=0.3, q=0.95, eip=11,
                          nu=1, eggsPerBatch=60){

  stopifnot(is.matrix(calK))
  stopifnot(dim(calK) == c(pars$nPatches, pars$nPatches))

  with(MYZopts,{
    MYZpar <- list()

    MYZpar$xde <- solve_as
    class(MYZpar$xde) <- solve_as
    if(solve_as == 'dde') class(MYZpar) <- c('RM', 'RM_dde')
    if(solve_as == 'ode') class(MYZpar) <- c('RM', 'RM_ode')

    MYZpar$g0      <- checkIt(g, pars$nPatches)
    MYZpar$sigma0  <- checkIt(sigma, pars$nPatches)
    MYZpar$f0      <- checkIt(f, pars$nPatches)
    MYZpar$q0      <- checkIt(q, pars$nPatches)
    MYZpar$nu0     <- checkIt(nu, pars$nPatches)
    MYZpar$eggsPerBatch <- eggsPerBatch
    MYZpar$eip <- eip
    MYZpar$calK <- calK

    pars$MYZpar = MYZpar
    pars = MBionomics.RM(0, 0, pars)

    return(pars)
  })}

#' @title Make inits for RM adult mosquito model
#' @param pars a [list]
#' @param MYZopts a [list] of values that overwrites the defaults
#' @param M0 total mosquito density at each patch
#' @param P0 total parous mosquito density at each patch
#' @param Y0 infected mosquito density at each patch
#' @param Z0 infectious mosquito density at each patch
#' @return a [list]
#' @export
make_MYZinits_RM = function(pars, MYZopts = list(),
                            M0=5, P0=1, Y0=1, Z0=1){
  with(MYZopts,{
    inits = list()
    inits$M0 = checkIt(M0, pars$nPatches)
    inits$P0 = checkIt(P0, pars$nPatches)
    inits$Y0 = checkIt(Y0, pars$nPatches)
    inits$Z0 = checkIt(Z0, pars$nPatches)

    pars$MYZinits = inits
    return(pars)
  })}

#' @title Add indices for adult mosquitoes to parameter list
#' @description Implements [make_indices_MYZ] for the RM model.
#' @inheritParams make_indices_MYZ
#' @return none
#' @importFrom utils tail
#' @export
make_indices_MYZ.RM_ode <- function(pars) {

  pars$MYZpar$M_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$MYZpar$M_ix, 1)

  pars$MYZpar$P_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$MYZpar$P_ix, 1)

  pars$MYZpar$Y_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$MYZpar$Y_ix, 1)

  pars$MYZpar$Z_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$MYZpar$Z_ix, 1)

  return(pars)
}

#' @title Add indices for adult mosquitoes to parameter list
#' @description Implements [make_indices_MYZ] for the RM model.
#' @inheritParams make_indices_MYZ
#' @return none
#' @importFrom utils tail
#' @export
make_indices_MYZ.RM_dde <- function(pars) {

  pars$MYZpar$M_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$MYZpar$M_ix, 1)

  pars$MYZpar$P_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$MYZpar$P_ix, 1)

  pars$MYZpar$Y_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$MYZpar$Y_ix, 1)

  pars$MYZpar$Z_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$MYZpar$Z_ix, 1)

  pars$MYZpar$Upsilon_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches^2)
  pars$max_ix <- tail(pars$MYZpar$Upsilon_ix, 1)

  pars$MYZpar$fqkappa_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$MYZpar$fqkappa_ix, 1)

  pars$MYZpar$g_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$MYZpar$g_ix, 1)

  pars$MYZpar$sigma_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$MYZpar$sigma_ix, 1)

  return(pars)
}

#' @title Make parameters for RM ODE adult mosquito model
#' @param pars a [list]
#' @param g mosquito mortality rate
#' @param sigma emigration rate
#' @param calK mosquito dispersal matrix of dimensions `nPatches` by `nPatches`
#' @param f feeding rate
#' @param q human blood fraction
#' @param nu oviposition rate, per mosquito
#' @param eggsPerBatch eggs laid per oviposition
#' @param eip length of extrinsic incubation period
#' @param solve_as is either `ode` to solve as an ode or `dde` to solve as a dde
#' @return none
#' @export
make_parameters_MYZ_RM <- function(pars, g, sigma, f, q, nu, eggsPerBatch, eip, calK, solve_as = 'dde') {
  stopifnot(is.numeric(g), is.numeric(sigma), is.numeric(f), is.numeric(q), is.numeric(nu), is.numeric(eggsPerBatch))

  MYZpar <- list()
  MYZpar$xde = solve_as
  class(MYZpar$xde) <- solve_as

  if(solve_as == 'dde') class(MYZpar) <- c('RM', 'RM_dde')
  if(solve_as == 'ode') class(MYZpar) <- c('RM', 'RM_ode')

  MYZpar$g0 <- g
  MYZpar$sigma0 <- sigma
  MYZpar$f0 <- f
  MYZpar$q0 <- q
  MYZpar$nu0 <- nu
  MYZpar$eggsPerBatch <- eggsPerBatch
  MYZpar$eip <- eip
  MYZpar$calK <- calK

  pars$MYZpar <- MYZpar
  pars = MBionomics.RM(0, 0, pars)
  return(pars)
}

#' @title Parse the output of deSolve and return variables for the RM model
#' @description Implements [parse_deout_MYZ] for the RM model
#' @inheritParams parse_deout_MYZ
#' @return a [list]
#' @export
parse_deout_MYZ.RM <- function(deout, pars) {
  M = deout[,pars$MYZpar$M_ix+1]
  P = deout[,pars$MYZpar$P_ix+1]
  Y = deout[,pars$MYZpar$Y_ix+1]
  Z = deout[,pars$MYZpar$Z_ix+1]
  y = Y/M
  z = Z/M
  parous = P/M
  return(list(M=M, P=P, Y=Y, Z=Z, y=y, z=z, parous))
}

#' @title Make inits for RM adult mosquito model
#' @param pars a [list]
#' @param y0 a vector of initial values
#' @return none
#' @export
update_inits_MYZ.RM_ode <- function(pars, y0) {
  M0 = y0[pars$MYZpar$M_ix]
  P0 = y0[pars$MYZpar$P_ix]
  Y0 = y0[pars$MYZpar$Y_ix]
  Z0 = y0[pars$MYZpar$Z_ix]
  pars = make_inits_MYZ_RM_ode(pars, M0, P0, Y0, Z0)
  return(pars)
}

#' @title Make inits for RM adult mosquito model
#' @param pars a [list]
#' @param M0 total mosquito density at each patch
#' @param P0 total parous mosquito density at each patch
#' @param Y0 infected mosquito density at each patch
#' @param Z0 infectious mosquito density at each patch
#' @return none
#' @export
make_inits_MYZ_RM_ode <- function(pars, M0, P0, Y0, Z0) {
  pars$MYZinits = list(M0=M0, P0=P0, Y0=Y0, Z0=Z0)
  return(pars)
}

#' @title Make inits for RM adult mosquito model
#' @param pars a [list]
#' @param M0 total mosquito density at each patch
#' @param P0 total parous mosquito density at each patch
#' @param Y0 infected mosquito density at each patch
#' @param Z0 infectious mosquito density at each patch
#' @param Upsilon0 the initial values of Upsilon
#' @return none
#' @export
make_inits_MYZ_RM_dde <- function(pars, M0, P0, Y0, Z0, Upsilon0) {
  pars$MYZinits = list(M0=M0, P0=P0, Y0=Y0, Z0=Z0, Upsilon0=Upsilon0, dummy=rep(0, 3*pars$nPatches))
  return(pars)
}

#' @title Update inits for RM adult mosquito model
#' @param pars a [list]
#' @param y0 a vector of initial values
#' @return none
#' @export
update_inits_MYZ.RM_dde <- function(pars, y0) {
  M0 = y0[pars$MYZpar$M_ix]
  P0 = y0[pars$MYZpar$P_ix]
  Y0 = y0[pars$MYZpar$Y_ix]
  Z0 = y0[pars$MYZpar$Z_ix]
  Upsilon0 = y0[pars$MYZpar$Upsilon_ix]
  pars = make_inits_MYZ_RM_dde(pars, M0, P0, Y0, Z0, Upsilon0)
  return(pars)
}

#' @title Return initial values as a vector
#' @description Implements [get_inits_MYZ] for the RM model.
#' @inheritParams get_inits_MYZ
#' @return none
#' @export
get_inits_MYZ.RM_ode <- function(pars) {with(pars$MYZinits,{
  c(M0, P0, Y0, Z0)
})}

#' @title Return initial values as a vector
#' @description Implements [get_inits_MYZ] for the RM model.
#' @inheritParams get_inits_MYZ
#' @return none
#' @export
get_inits_MYZ.RM_dde <- function(pars) {with(pars$MYZinits,{
  c(M0, P0, Y0, Z0, as.vector(Upsilon0), rep(0, 3*pars$nPatches))
})}
