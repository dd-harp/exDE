# specialized methods for the adult mosquito GeRM model

#' @title Set the availability of resources
#' @description Implements [ResourceAvailability] for the GeRM model
#' @inheritParams ResourceAvailability
#' @return a named [list]
#' @export
ResourceAvailability.GeRM <- function(t, y, pars) {
  pars$O = F_other(t, pars)
  pars$S = F_sugar(t, pars)
  pars$W = computeW(t, y, pars)
  pars$B = computeB(t, pars)
  pars$local_frac = with(pars, W/(W + Visitors))
  pars$Q = computeQ(t, pars)
  return(pars)
}

#' @title Reset bloodfeeding and mortality rates to baseline
#' @description Implements [MosquitoBehavior] for the GeRM model
#' @inheritParams MosquitoBehavior
#' @return a named [list]
#' @export
MosquitoBehavior.GeRM <- function(t, y, pars) {
  with(pars,{
    pars$MYZpar$f = F_f(t, pars)
    pars$MYZpar$q = F_q(t, pars)
    pars$MYZpar$g = F_g(t, pars)
    pars$MYZpar$sigma = F_sigma(t, pars)
    pars$MYZpar$nu = F_nu(t, pars)
    pars$MYZpar$eip = F_eip(t, pars)
    return(pars)
})}

#' @title Density of infectious mosquitoes
#' @description Implements [F_Z] for the GeRM model.
#' @inheritParams F_Z
#' @return a [numeric] vector of length `nPatches`
#' @export
F_Z.GeRM <- function(t, y, pars) {
  y[pars$MYZpar$Z_ix]
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
  pars = MosquitoBehavior.RM(0, 0, pars)
  return(pars)
}

#' @title Set up a static exogenous forcing for the GeRM ODE adult mosquito model
#' @param pars a [list]
#' @param other is the availability of other blood hosts
#' @param W is the availbility of the pathogen's hosts
#' @param zeta is a shape parameter
#' @param Q is the availability of aquatic habitats
#' @param sugar is sugar availability
#' @return none
#' @export
setup_forcing_MYZ_GeRM_basic <- function(pars, other, sugar, W, zeta, Q) {
  pars = make_parameters_exogenous_forced(pars)
  RApar = list()
  class(RApar) = "GeRM"
  pars$RApar = RApar

  Wpar = list()
  class(Wpar) <- "static"
  pars$Wpar = Wpar
  pars$W = W
  #pars$W = pars$TaR %*% (pars$Hpar$wts_f * F_H(t, y, pars))

  OBpar = list()
  class(OBpar) <- "static"
  pars$OBpar = OBpar
  pars$other = other

  Bpar = list()
  class(Bpar) <- "static"
  Bpar$zeta = zeta
  pars$Bpar = Bpar
  pars$B =  pars$W + pars$other^zeta

  SGRpar = list()
  class(SGRpar) <- "static"
  pars$SGRpar = SGRpar
  pars$sugar = sugar

  Qpar = list()
  class(Qpar) <- "static"
  pars$Qpar = Qpar
  pars$Q = Q

  pars = ResourceAvailability(0, 0, pars)
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
#' @param pars an [environment]
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

