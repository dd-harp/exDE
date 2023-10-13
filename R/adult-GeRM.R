# specialized methods for the adult mosquito GeRM model

#' @title Reset bloodfeeding and mortality rates to baseline
#' @description Implements [MBionomics] for the GeRM model
#' @inheritParams MBionomics
#' @return a named [list]
#' @export
MBionomics.GeRM <- function(t, y, pars) {
  with(pars,{
    pars$MYZpar$f = F_f(t, pars)
    pars$MYZpar$q = F_q(t, pars)
    pars$MYZpar$g = F_g(t, pars)
    pars$MYZpar$sigma = F_sigma(t, pars)
    pars$MYZpar$nu = F_nu(t, pars)
    return(pars)
})}


#' @title Blood feeding rate of the infective mosquito population
#' @description Implements [F_fqZ] for the GeRM model.
#' @inheritParams F_fqZ
#' @return a [numeric] vector of length `nPatches`
#' @export
F_fqZ.GeRM <- function(t, y, pars) {
  with(pars$MYZpar, f*q)*y[pars$MYZpar$Z_ix]
}

#' @title Number of eggs laid by adult mosquitoes
#' @description Implements [F_eggs] for the GeRM model.
#' @inheritParams F_eggs
#' @return a [numeric] vector of length `nPatches`
#' @export
F_eggs.GeRM <- function(t, y, pars) {
  with(pars$MYZpar, {
    G <- y[G_ix]
    return(G*nu*eggsPerBatch)
  })
}

#' @title Derivatives for adult mosquitoes
#' @description Implements [dMYZdt] for the GeRM ODE model.
#' @inheritParams dMYZdt
#' @return a [numeric] vector
#' @export
dMYZdt.GeRM_ode <- function(t, y, pars, Lambda, kappa) {

  nPatches <- pars$nPatches

  with(pars$MYZpar,{

    M <- y[M_ix]
    G <- y[G_ix]
    Y <- y[Y_ix]
    Z <- y[Z_ix]

    Omega <- make_Omega(g, sigma, calK, nPatches)
    Upsilon <- expm(-Omega*eip)

    dMdt <- Lambda - (Omega %*% M)
    dGdt <- f*(M - G) - nu*G - (Omega %*% G)
    dYdt <- f*q*kappa*(M - Y) - (Omega %*% Y)
    dZdt <- Upsilon %*% diag(f*q*kappa, nPatches) %*% (M - Y) - (Omega %*% Z)

    return(c(dMdt, dGdt, dYdt, dZdt))
  })
}

#' @title Derivatives for adult mosquitoes
#' @description Implements [dMYZdt] for the GeRM DDE model.
#' @inheritParams dMYZdt
#' @return a [numeric] vector
#' @importFrom deSolve lagvalue
#' @importFrom deSolve lagderiv
#' @export
dMYZdt.GeRM_dde <- function(t, y, pars, Lambda, kappa) {

  nPatches <- pars$nPatches
  eip <- pars$MYZpar$eip

  with(pars$MYZpar,{
    if (t < eip) {
      M_eip <- pars$MYZinits$M0
      Y_eip <- pars$MYZinits$Y0
      kappa_eip <- kappa
      f_eip <- f
      q_eip <- q
      g_eip <- g
      sigma_eip <- sigma
    } else {
      M_eip <- lagvalue(t = t - eip, nr = M_ix)
      Y_eip <- lagvalue(t = t - eip, nr = Y_ix)
      kappa_eip <- lagderiv(t = t-eip, nr = kappa_ix)
      f_eip <- lagderiv(t = t-eip, nr = f_ix)
      q_eip <- lagderiv(t = t-eip, nr = q_ix)
      g_eip <- lagderiv(t = t-eip, nr = g_ix)
      sigma_eip <- lagderiv(t = t-eip, nr = sigma_ix)
    }

    M <- y[M_ix]
    G <- y[G_ix]
    Y <- y[Y_ix]
    Z <- y[Z_ix]
    Upsilon <- matrix(data = y[Upsilon_ix], nrow = nPatches, ncol = nPatches)

    Omega <- make_Omega(g, sigma, calK, nPatches)
    Omega_eip <- make_Omega(g_eip, sigma_eip, calK, nPatches)

    dMdt <- Lambda - (Omega %*% M)
    dGdt <- f*(M - G) - nu*G - (Omega %*% G)
    dYdt <- f*q*kappa*(M - Y) - (Omega %*% Y)
    dZdt <- Upsilon %*% diag(f_eip*q_eip*kappa_eip, nPatches) %*% (M_eip - Y_eip) - (Omega %*% Z)
    dUdt <- as.vector((Omega_eip - Omega) %*% Upsilon)

    return(c(dMdt, dGdt, dYdt, dZdt, dUdt, kappa, f, q, g, sigma))
  })
}

#' @title Setup the GeRM model for adult mosquitoes
#' @description Implements [setup_MYZ] for the GeRM model
#' @inheritParams setup_MYZ
#' @return a [list] vector
#' @export
setup_MYZ.GeRM = function(pars, MYZname,
                             nPatches=1, MYZopts=list(),
                             calK=diag(1)){

  pars$MYZname = "GeRM"
  pars$nPatches = checkIt(nPatches, 1, "integer")

  pars = make_MYZpar_GeRM(pars, MYZopts, calK)
  pars = make_MYZinits_GeRM(pars, MYZopts)

  if(pars$MYZpar$xde == "dde"){
    pars$xde = "dde"
    class(pars$xde) = "dde"
    Omega <- with(pars$MYZpar, make_Omega(g, sigma, calK, nPatches))
    Upsilon <- expm::expm(-Omega*pars$MYZpar$eip)
    pars$MYZinits$Upsilon0 = as.vector(Upsilon)
  }
  return(pars)
}


#' @title Make parameters for a GeRM ODE adult mosquito model
#' @param pars a [list]
#' @param MYZopts a [list] to overwrite defaults
#' @param calK mosquito dispersal matrix of dimensions `nPatches` by `nPatches`
#' @param g mosquito mortality rate
#' @param setup_Fg a [list] to set up F_g
#' @param sigma emigration rate
#' @param setup_Fsigma a [list] to set up F_sigma
#' @param f feeding rate
#' @param setup_Ff a [list] to set up F_f
#' @param q human blood fraction
#' @param setup_Fq a [list] to set up F_q
#' @param nu oviposition rate, per mosquito
#' @param setup_Fnu a [list] to set up F_nu
#' @param eip length of extrinsic incubation period
#' @param setup_Feip a [list] to set up F_eip
#' @param eggsPerBatch eggs laid per oviposition
#' @param solve_as is either `ode` to solve as an ode or `dde` to solve as a dde
#' @return none
#' @export
make_MYZpar_GeRM = function(pars, MYZopts=list(), calK,
                            solve_as = "dde",
                            g=1/12, setup_Fg = list(),
                            sigma=1/8, setup_Fsigma = list(),
                            f=0.3, setup_Ff = list(),
                            q=0.95, setup_Fq = list(),
                            nu=1, setup_Fnu = list(),
                            eip=11, setup_Feip = list(),
                            eggsPerBatch=60){

  stopifnot(is.matrix(calK))
  stopifnot(dim(calK) == c(pars$nPatches, pars$nPatches))

  with(MYZopts,{
    MYZpar <- list()

    MYZpar$xde <- solve_as
    class(MYZpar$xde) <- solve_as
    if(solve_as == 'dde') class(MYZpar) <- c('GeRM', 'GeRM_dde')
    if(solve_as == 'ode') class(MYZpar) <- c('GeRM', 'GeRM_ode')

    if(length(setup_Fg) == 0){
      MYZpar$g_par <- list()
      class(MYZpar$g_par) <- "static"
      MYZpar$g0 <- checkIt(g, pars$nPatches)
    } else MYZpar$g_par <- setup_Fx(setup_Fg)

    if(length(setup_Fsigma) == 0){
      MYZpar$sigma_par <- list()
      class(MYZpar$sigma_par) <- "static"
      MYZpar$sigma0 <- checkIt(sigma, pars$nPatches)
    } else MYZpar$sigma_par <- setup_Fx(setup_Fsigma)

    if(length(setup_Ff) == 0){
      MYZpar$f_par <- list()
      class(MYZpar$f_par) <- "static"
      MYZpar$f0 <- checkIt(f, pars$nPatches)
    } else MYZpar$f_par = setup_Fx(setup_Ff)

    if(length(setup_Fq) == 0){
      MYZpar$q_par <- list()
      class(MYZpar$q_par) <- "static"
      MYZpar$q0 <- checkIt(q, pars$nPatches)
    } else MYZpar$q_par <- setup_Fx(setup_Fq)

    if(length(setup_Fnu) == 0){
      MYZpar$nu_par <- list()
      class(MYZpar$nu_par) <- "static"
      MYZpar$nu0 <- checkIt(nu, pars$nPatches)
    } else MYZpar$nu_par= setup_Fx(setup_Fnu)

    MYZpar$eip <- eip

    MYZpar$eggsPerBatch <- eggsPerBatch

    MYZpar$calK <- calK

    pars$MYZpar = MYZpar
    pars = MBionomics(0, 0, pars)

    return(pars)
})}

#' @title Make inits for GeRM adult mosquito model
#' @param pars a [list]
#' @param MYZopts a [list] that overwrites the defaults
#' @param M0 total mosquito density at each patch
#' @param G0 total parous mosquito density at each patch
#' @param Y0 infected mosquito density at each patch
#' @param Z0 infectious mosquito density at each patch
#' @return none
#' @export
make_MYZinits_GeRM = function(pars, MYZopts = list(),
                              M0=5, G0=1, Y0=1, Z0=1){
  with(MYZopts,{
    inits = list()
    inits$M0 = checkIt(M0, pars$nPatches)
    inits$G0 = checkIt(G0, pars$nPatches)
    inits$Y0 = checkIt(Y0, pars$nPatches)
    inits$Z0 = checkIt(Z0, pars$nPatches)

    pars$MYZinits = inits
    return(pars)
  })}

#' @title Add indices for adult mosquitoes to parameter list
#' @description Implements [make_indices_MYZ] for the GeRM model.
#' @inheritParams make_indices_MYZ
#' @return none
#' @importFrom utils tail
#' @export
make_indices_MYZ.GeRM_ode <- function(pars) {

  pars$MYZpar$M_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$MYZpar$M_ix, 1)

  pars$MYZpar$G_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$MYZpar$G_ix, 1)

  pars$MYZpar$Y_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$MYZpar$Y_ix, 1)

  pars$MYZpar$Z_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$MYZpar$Z_ix, 1)

  return(pars)
}

#' @title Add indices for adult mosquitoes to parameter list
#' @description Implements [make_indices_MYZ] for the GeRM model.
#' @inheritParams make_indices_MYZ
#' @return none
#' @importFrom utils tail
#' @export
make_indices_MYZ.GeRM_dde <- function(pars) {

  pars$MYZpar$M_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$MYZpar$M_ix, 1)

  pars$MYZpar$G_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$MYZpar$G_ix, 1)

  pars$MYZpar$Y_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$MYZpar$Y_ix, 1)

  pars$MYZpar$Z_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$MYZpar$Z_ix, 1)

  pars$MYZpar$Upsilon_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches^2)
  pars$max_ix <- tail(pars$MYZpar$Upsilon_ix, 1)

  pars$MYZpar$kappa_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$MYZpar$kappa_ix, 1)

  pars$MYZpar$f_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$MYZpar$f_ix, 1)

  pars$MYZpar$q_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$MYZpar$q_ix, 1)

  pars$MYZpar$g_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$MYZpar$g_ix, 1)

  pars$MYZpar$sigma_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$MYZpar$sigma_ix, 1)

  return(pars)
}

#' @title Parse the output of deSolve and return variables for the GeRM model
#' @description Implements [parse_deout_MYZ] for the GeRM model.
#' @inheritParams parse_deout_MYZ
#' @return none
#' @export
parse_deout_MYZ.GeRM <- function(varslist, deout, pars) {
  varslist$M = deout[,pars$MYZpar$M_ix+1]
  varslist$G = deout[,pars$MYZpar$G_ix+1]
  varslist$Y = deout[,pars$MYZpar$Y_ix+1]
  varslist$Z = deout[,pars$MYZpar$Z_ix+1]
  varslist$g = with(varslist, g/M)
  varslist$y = with(varslist, Y/M)
  varslist$z = with(varslist, Z/M)
  return(varslist)
}

#' @title Make parameters for a GeRM ODE adult mosquito model
#' @param pars a [list]
#' @param g mosquito mortality rate
#' @param sigma emigration rate
#' @param f feeding rate
#' @param q human blood fraction
#' @param nu oviposition rate, per mosquito
#' @param eggsPerBatch eggs laid per oviposition
#' @param eip length of extrinsic incubation period
#' @param calK mosquito dispersal matrix of dimensions `nPatches` by `nPatches`
#' @param solve_as is either `ode` to solve as an ode or `dde` to solve as a dde
#' @return none
#' @export
make_parameters_MYZ_GeRM_static <- function(pars, g, sigma, f, q, nu, eggsPerBatch, eip, calK, solve_as = 'dde') {
  stopifnot(is.numeric(g), is.numeric(sigma), is.numeric(f), is.numeric(q), is.numeric(nu), is.numeric(eggsPerBatch))

  MYZpar <- list()
  MYZpar$xde = solve_as
  class(MYZpar$xde) <- solve_as

  if(solve_as == 'dde') class(MYZpar) <- c('GeRM', 'GeRM_dde')
  if(solve_as == 'ode') class(MYZpar) <- c('GeRM', 'GeRM_ode')

  MYZpar$f_par <- list()
  class(MYZpar$f_par) <- "static"
  MYZpar$q_par <- list()
  class(MYZpar$q_par) <- "static"
  MYZpar$g_par <- list()
  class(MYZpar$g_par) <- "static"
  MYZpar$sigma_par <- list()
  class(MYZpar$sigma_par) <- "static"
  MYZpar$nu_par <- list()
  class(MYZpar$nu_par) <- "static"
  MYZpar$eip_par <- list()
  class(MYZpar$eip_par) <- "static"

  MYZpar$g0 <- g
  MYZpar$f0 <- f
  MYZpar$q0 <- q
  MYZpar$sigma0 <- sigma
  MYZpar$nu0 <- nu
  MYZpar$eggsPerBatch <- eggsPerBatch
  MYZpar$eip <- eip
  MYZpar$calK <- calK

  pars$MYZpar <- MYZpar
  pars = MBionomics(0, 0, pars)
  return(pars)
}

#' @title Make inits for GeRM adult mosquito model
#' @param pars a [list]
#' @param M0 total mosquito density at each patch
#' @param G0 total parous mosquito density at each patch
#' @param Y0 infected mosquito density at each patch
#' @param Z0 infectious mosquito density at each patch
#' @return none
#' @export
make_inits_MYZ_GeRM_ode <- function(pars, M0, G0, Y0, Z0) {
  pars$MYZinits = list(M0=M0, G0=G0, Y0=Y0, Z0=Z0)
  return(pars)
}

#' @title Make inits for GeRM adult mosquito model
#' @param pars a [list]
#' @param M0 total mosquito density at each patch
#' @param G0 total parous mosquito density at each patch
#' @param Y0 infected mosquito density at each patch
#' @param Z0 infectious mosquito density at each patch
#' @param Upsilon0 the initial values of Upsilon
#' @return none
#' @export
make_inits_MYZ_GeRM_dde <- function(pars, M0, G0, Y0, Z0, Upsilon0) {
  pars$MYZinits = list(M0=M0, G0=G0, Y0=Y0, Z0=Z0, Upsilon0=Upsilon0, rep(0, 5*pars$nPatches))
  return(pars)
}

#' @title Return initial values as a vector
#' @description Implements [get_inits_MYZ] for the GeRM model.
#' @inheritParams get_inits_MYZ
#' @return none
#' @export
get_inits_MYZ.GeRM_ode <- function(pars) {with(pars$MYZinits,{
  c(M0, G0, Y0, Z0)
})}

#' @title Return initial values as a vector
#' @description Implements [get_inits_MYZ] for the GeRM model.
#' @inheritParams get_inits_MYZ
#' @return none
#' @export
get_inits_MYZ.GeRM_dde <- function(pars) {with(pars$MYZinits,{
  c(M0, G0, Y0, Z0, as.vector(Upsilon0), rep(0, 5*pars$nPatches))
})}

