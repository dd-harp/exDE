# specialized methods for the adult mosquito RM model

#' @title Reset bloodfeeding and mortality rates to baseline
#' @description Implements [MBionomics] for the RM model
#' @inheritParams MBionomics
#' @return the model as a [list]
#' @export
MBionomics.RM <- function(t, y, pars, s) {
  with(pars$MYZpar[[s]],{
    pars$MYZpar[[s]]$f <- f0
    pars$MYZpar[[s]]$q <- q0
    pars$MYZpar[[s]]$g <- g0
    pars$MYZpar[[s]]$sigma <- sigma0
    pars$MYZpar[[s]]$nu <- nu0
    pars$MYZpar[[s]]$eip<- EIP(t, EIPmod)

    return(pars)
})}

#' @title The net blood feeding rate of the infective mosquito population in a patch
#' @description Implements [F_fqZ] for the RM model.
#' @inheritParams F_fqZ
#' @return a [numeric] vector of length `nPatches`
#' @export
F_fqZ.RM <- function(t, y, pars, s) {
  with(pars$MYZpar[[s]], f*q)*y[pars$ix$MYZ[[s]]$Z_ix]
}

#' @title Number of eggs laid by adult mosquitoes
#' @description Implements [F_eggs] for the RM model.
#' @inheritParams F_eggs
#' @return a [numeric] vector of length `nPatches`
#' @export
F_eggs.RM <- function(t, y, pars, s) {
  M <- y[pars$ix$MYZ[[s]]$M_ix]
  with(pars$MYZpar[[s]],{
    return(M*nu*eggsPerBatch)
  })
}

#' @title Derivatives for adult mosquitoes
#' @description Implements [dMYZdt] for the RM ODE model.
#' @inheritParams dMYZdt
#' @return a [numeric] vector
#' @export
dMYZdt.RM_ode <- function(t, y, pars, s) {
  Lambda = pars$Lambda[[s]]
  kappa = pars$kappa [[s]]

  with(pars$ix$MYZ[[s]],{
    M <- y[M_ix]
    P <- y[P_ix]
    Y <- y[Y_ix]
    Z <- y[Z_ix]

    with(pars$MYZpar[[s]],{
      dMdt <- Lambda - (Omega %*% M)
      dPdt <- f*(M - P) - (Omega %*% P)
      dYdt <- f*q*kappa*(M - Y) - (Omega %*% Y)
      dZdt <- Upsilon %*% diag(f*q*kappa, nPatches) %*% (M - Y) - (Omega %*% Z)

      return(c(dMdt, dPdt, dYdt, dZdt))
    })
  })
}

#' @title Derivatives for adult mosquitoes
#' @description Implements [dMYZdt] for the RM DDE model.
#' @inheritParams dMYZdt
#' @return a [numeric] vector
#' @importFrom deSolve lagvalue
#' @importFrom deSolve lagderiv
#' @export
dMYZdt.RM_dde <- function(t, y, pars, s){

  Lambda = pars$Lambda[[s]]
  kappa = pars$kappa [[s]]

  with(pars$ix$MYZ[[s]],{
    M <- y[M_ix]
    P <- y[P_ix]
    Y <- y[Y_ix]
    Z <- y[Z_ix]
    Upsilon <- y[Upsilon_ix]

    with(pars$MYZpar[[s]],{

      if (t <= eip) {
        M_eip <- pars$MYZinits[[s]]$M
        Y_eip <- pars$MYZinits[[s]]$Y
        fqkappa_eip <- kappa*f*q
        g_eip <- g
        sigma_eip <- sigma
      } else {
        M_eip <- lagvalue(t = t - eip, nr = M_ix)
        Y_eip <- lagvalue(t = t - eip, nr = Y_ix)
        fqkappa_eip <- lagderiv(t = t-eip, nr = fqkappa_ix)
        g_eip <- lagderiv(t = t-eip, nr = g_ix)
        sigma_eip <- lagderiv(t = t-eip, nr = sigma_ix)
      }

      Omega <- make_Omega(g, sigma, calK, nPatches)
      Omega_eip <- make_Omega(g_eip, sigma_eip, calK, nPatches)
      dMdt <- Lambda - (Omega %*% M)
      dPdt <- f*(M - P) - (Omega %*% P)
      dYdt <- f*q*kappa*(M - Y) - (Omega %*% Y)
      dZdt <- Upsilon %*% (fqkappa_eip * (M_eip - Y_eip)) - (Omega %*% Z)
      dUdt <- as.vector(((1-dEIPdt(t,EIPmod))*Omega_eip - Omega) %*% Upsilon)

      return(c(dMdt, dPdt, dYdt, dZdt, dUdt, f*q*kappa, g, sigma))
    })
  })
}


#' @title Setup MYZpar for the RM model
#' @description Implements [setup_MYZpar] for the RM model
#' @inheritParams setup_MYZpar
#' @return a [list] vector
#' @export
setup_MYZpar.RM = function(MYZname, pars, s, MYZopts=list(), EIPmod, calK){
  pars$MYZpar[[s]] = make_MYZpar_RM(pars$nPatches, MYZopts, EIPmod, calK)
  return(pars)
}

#' @title Make parameters for RM ODE adult mosquito model
#' @param nPatches is the number of patches, an integer
#' @param MYZopts a [list] of values that overwrites the defaults
#' @param EIPmod a [list] that defines the EIP model
#' @param calK a mosquito dispersal matrix of dimensions `nPatches` by `nPatches`
#' @param g mosquito mortality rate
#' @param sigma emigration rate
#' @param f feeding rate
#' @param q human blood fraction
#' @param nu oviposition rate, per mosquito
#' @param eggsPerBatch eggs laid per oviposition
#' @return a [list]
#' @export
make_MYZpar_RM = function(nPatches, MYZopts=list(), EIPmod, calK,
                          g=1/12, sigma=1/8, f=0.3, q=0.95,
                          nu=1, eggsPerBatch=60){

  stopifnot(is.matrix(calK))
  stopifnot(dim(calK) == c(nPatches, nPatches))

  with(MYZopts,{
    MYZpar <- list()
    if(!exists("solve_as")) solve_as = "dde"
    MYZpar$xde <- solve_as
    class(MYZpar$xde) <- solve_as

    if(solve_as == "ode") class(MYZpar) <- c('RM', 'RM_ode')
    if(solve_as == "dde") class(MYZpar) <- c('RM', 'RM_dde')

    MYZpar$nPatches <- nPatches

    MYZpar$g      <- checkIt(g, nPatches)
    MYZpar$sigma  <- checkIt(sigma, nPatches)
    MYZpar$f      <- checkIt(f, nPatches)
    MYZpar$q      <- checkIt(q, nPatches)
    MYZpar$nu     <- checkIt(nu, nPatches)
    MYZpar$eggsPerBatch <- eggsPerBatch

    # Store as baseline values
    MYZpar$g0      <- MYZpar$g
    MYZpar$sigma0  <- MYZpar$sigma
    MYZpar$f0      <- MYZpar$f
    MYZpar$q0      <- MYZpar$q
    MYZpar$nu0     <- MYZpar$nu

    # The EIP model and the eip
    MYZpar$EIPmod <- EIPmod
    MYZpar$eip <- EIP(0, EIPmod)

    MYZpar$calK <- calK

    MYZpar$Omega <- make_Omega(g, sigma, calK, nPatches)
    MYZpar$Upsilon <- with(MYZpar, expm::expm(-Omega*eip))

    return(MYZpar)
})}

#' @title Setup initial values for the RM model
#' @description Implements [setup_MYZinits] for the RM model
#' @inheritParams setup_MYZinits
#' @return a [list]
#' @export
setup_MYZinits.RM = function(pars, s, MYZopts=list()){
  pars$MYZinits[[s]] = with(pars$MYZpar[[s]], make_MYZinits_RM_dde(nPatches, Upsilon, MYZopts))
  return(pars)
}


#' @title Make inits for RM adult mosquito model
#' @param nPatches the number of patches in the model
#' @param Upsilon a matrix describing survival and dispersal through the EIP
#' @param MYZopts a [list] of values that overwrites the defaults
#' @param M0 total mosquito density at each patch
#' @param P0 total parous mosquito density at each patch
#' @param Y0 infected mosquito density at each patch
#' @param Z0 infectious mosquito density at each patch
#' @return a [list]
#' @export
make_MYZinits_RM_dde = function(nPatches, Upsilon, MYZopts = list(),
                            M0=5, P0=1, Y0=1, Z0=1){
  with(MYZopts,{
    M = checkIt(M0, nPatches)
    P = checkIt(P0, nPatches)
    Y = checkIt(Y0, nPatches)
    Z = checkIt(Z0, nPatches)
    Upsilon = checkIt(as.vector(Upsilon), nPatches^2)
    dd = rep(0, nPatches)
    return(list(M=M, P=P, Y=Y, Z=Z, Upsilon=as.vector(Upsilon), d1=dd, d2=dd, d3=dd))
  })
}

#' @title Make inits for RM adult mosquito model
#' @param nPatches the number of patches in the model
#' @param MYZopts a [list] of values that overwrites the defaults
#' @param M0 total mosquito density at each patch
#' @param P0 total parous mosquito density at each patch
#' @param Y0 infected mosquito density at each patch
#' @param Z0 infectious mosquito density at each patch
#' @return a [list]
#' @export
make_MYZinits_RM_ode = function(nPatches, MYZopts = list(),
                                M0=5, P0=1, Y0=1, Z0=1){
  with(MYZopts,{
    M = checkIt(M0, nPatches)
    P = checkIt(P0, nPatches)
    Y = checkIt(Y0, nPatches)
    Z = checkIt(Z0, nPatches)
    return(list(M=M, P=P, Y=Y, Z=Z))
  })
}


#' @title Add indices for adult mosquitoes to parameter list
#' @description Implements [make_indices_MYZ] for the RM model.
#' @inheritParams make_indices_MYZ
#' @return a [list]
#' @importFrom utils tail
#' @export
make_indices_MYZ.RM_dde <- function(pars, s) {with(pars,{

  M_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(M_ix, 1)

  P_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(P_ix, 1)

  Y_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(Y_ix, 1)

  Z_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(Z_ix, 1)

  Upsilon_ix <- seq(from = max_ix+1, length.out = nPatches^2)
  max_ix <- tail(Upsilon_ix, 1)

  fqkappa_ix <- seq(from = max_ix+1, length.out = nPatches)
  max_ix <- tail(fqkappa_ix, 1)

  g_ix <- seq(from = max_ix+1, length.out = nPatches)
  max_ix <- tail(g_ix, 1)

  sigma_ix <- seq(from = max_ix+1, length.out = nPatches)
  max_ix <- tail(sigma_ix, 1)

  pars$max_ix = max_ix
  pars$ix$MYZ[[s]] = list(M_ix=M_ix, P_ix=P_ix, Y_ix=Y_ix, Z_ix=Z_ix,
              Upsilon_ix = Upsilon_ix, fqkappa_ix=fqkappa_ix,
              g_ix=g_ix, sigma_ix=sigma_ix)
  return(pars)
})}


#' @title Add indices for adult mosquitoes to parameter list
#' @description Implements [make_indices_MYZ] for the RM model.
#' @inheritParams make_indices_MYZ
#' @return a [list]
#' @importFrom utils tail
#' @export
make_indices_MYZ.RM_ode <- function(pars, s) {with(pars,{

  M_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(M_ix, 1)

  P_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(P_ix, 1)

  Y_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(Y_ix, 1)

  Z_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(Z_ix, 1)

  pars$max_ix = max_ix
  pars$ix$MYZ[[s]] = list(M_ix=M_ix, P_ix=P_ix, Y_ix=Y_ix, Z_ix=Z_ix)
  return(pars)
})}

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
#' @return a [list]
#' @export
make_parameters_MYZ_RM <- function(pars, g, sigma, f, q, nu, eggsPerBatch, eip, calK, solve_as = "dde") {
  stopifnot(is.numeric(g), is.numeric(sigma), is.numeric(f),
            is.numeric(q), is.numeric(nu), is.numeric(eggsPerBatch))

  MYZpar <- list()
  MYZpar$xde = solve_as
  class(MYZpar$xde) <- solve_as

  if(solve_as == 'dde') class(MYZpar) <- c('RM', 'RM_dde')
  if(solve_as == 'ode') class(MYZpar) <- c('RM', 'RM_ode')

  MYZpar$g      <- checkIt(g, pars$nPatches)
  MYZpar$sigma  <- checkIt(sigma, pars$nPatches)
  MYZpar$f      <- checkIt(f, pars$nPatches)
  MYZpar$q      <- checkIt(q, pars$nPatches)
  MYZpar$nu     <- checkIt(nu, pars$nPatches)
  MYZpar$eggsPerBatch <- eggsPerBatch

  # Store as baseline values
  MYZpar$g0      <- MYZpar$g
  MYZpar$sigma0  <- MYZpar$sigma
  MYZpar$f0      <- MYZpar$f
  MYZpar$q0      <- MYZpar$q
  MYZpar$nu0     <- MYZpar$nu

  Omega   <- make_Omega(g, sigma, calK, pars$nPatches)
  MYZpar$Omega <- Omega
  MYZpar$Upsilon <- expm(-Omega*eip)
  MYZpar$EIPmod <- setup_eip_static(eip=eip)
  MYZpar$eip <- eip
  MYZpar$calK <- calK
  MYZpar$nPatches <- pars$nPatches

  pars$MYZpar = list()
  pars$MYZpar[[1]] = MYZpar

  return(pars)
}

#' @title Make inits for RM adult mosquito model
#' @param pars a [list]
#' @param M0 total mosquito density at each patch
#' @param P0 total parous mosquito density at each patch
#' @param Y0 infected mosquito density at each patch
#' @param Z0 infectious mosquito density at each patch
#' @return a [list]
#' @export
make_inits_MYZ_RM_ode <- function(pars, M0, P0, Y0, Z0) {
  pars$MYZinits = list()
  pars$MYZinits[[1]] = list(M=M0, P=P0, Y=Y0, Z=Z0)
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
  pars$MYZinits = list()
  dmy = rep(0, pars$nPatches)
  pars$MYZinits[[1]] = list(M=M0, P=P0, Y=Y0, Z=Z0, Upsilon=Upsilon0, d1=dmy, d2=dmy, d3=dmy)
  return(pars)
}

#' @title Parse the output of deSolve and return variables for the RM model
#' @description Implements [parse_deout_MYZ] for the RM model
#' @inheritParams parse_deout_MYZ
#' @return a [list]
#' @export
parse_deout_MYZ.RM <- function(deout, pars, s) {with(pars$ix$MYZ[[s]],{
  time = deout[,1]
  M = deout[,M_ix+1]
  P = deout[,P_ix+1]
  Y = deout[,Y_ix+1]
  Z = deout[,Z_ix+1]
  y = Y/M
  z = Z/M
  parous = P/M
  return(list(time=time, M=M, P=P, Y=Y, Z=Z, y=y, z=z, parous))
})}

#' @title Return initial values as a vector
#' @description Implements [get_inits_MYZ] for the RM model.
#' @inheritParams get_inits_MYZ
#' @return [numeric]
#' @export
get_inits_MYZ.RM_dde <- function(pars, s) {with(pars$MYZinits[[s]],{
  c(M, P, Y, Z, as.vector(Upsilon), d1, d2, d3)
})}

#' @title Return initial values as a vector
#' @description Implements [get_inits_MYZ] for the RM model.
#' @inheritParams get_inits_MYZ
#' @return [numeric]
#' @export
get_inits_MYZ.RM_ode <- function(pars, s) {with(pars$MYZinits[[s]],{
  c(M, P, Y, Z)
})}

#' @title Update inits for RM adult mosquito model
#' @inheritParams update_inits_MYZ
#' @return a [list]
#' @export
update_inits_MYZ.RM_dde <- function(pars, y0, s) {with(pars$ix$MYZ[[s]],{
  M = y0[M_ix]
  P = y0[P_ix]
  Y = y0[Y_ix]
  Z = y0[Z_ix]
  Upsilon = y0[Upsilon_ix]
  pars$MYZinits[[s]]= make_MYZinits_RM_dde(pars$nPatches, Upsilon, list(), M0=M, P0=P, Y0=Y, Z0=Z)
  return(pars)
})}

#' @title Make inits for RM adult mosquito model
#' @inheritParams update_inits_MYZ
#' @return a [list]
#' @export
update_inits_MYZ.RM_ode <- function(pars, y0, s) {with(pars$ix$MYZ[[s]],{
  M = y0[M_ix]
  P = y0[P_ix]
  Y = y0[Y_ix]
  Z = y0[Z_ix]
  pars$MYZinits[[s]] = make_MYZinits_RM_ode(pars$nPatches, list(), M0=M, P0=P, Y0=Y, Z0=Z)
  return(pars)
})}




