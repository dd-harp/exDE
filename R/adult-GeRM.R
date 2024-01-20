# specialized methods for the adult mosquito GeRM model

#' @title Reset bloodfeeding and mortality rates to baseline
#' @description Implements [MBionomics] for the GeRM model
#' @inheritParams MBionomics
#' @return the model as a [list]
#' @export
MBionomics.GeRM <- function(t, y, pars, s) {
  with(pars,{
    pars$MYZpar[[s]]$f     <- F_f(t, MYZpar[[s]])
    pars$MYZpar[[s]]$q     <- F_q(t, MYZpar[[s]])
    pars$MYZpar[[s]]$g     <- F_g(t, MYZpar[[s]])
    pars$MYZpar[[s]]$sigma <- F_sigma(t, MYZpar[[s]])
    pars$MYZpar[[s]]$nu    <- F_nu(t, MYZpar[[s]])
    pars$MYZpar[[s]]$eip   <- EIP(t, MYZpar[[s]]$EIPmod)
    return(pars)
  })
}

#' @title Blood feeding rate of the infective mosquito population
#' @description Implements [F_fqZ] for the GeRM model.
#' @inheritParams F_fqZ
#' @return a [numeric] vector of length `nPatches`
#' @export
F_fqZ.GeRM <- function(t, y, pars, s) {
  with(pars$MYZpar[[s]], f*q)*y[pars$ix$MYZ[[s]]$Z_ix]
}

#' @title Blood feeding rate of the infective mosquito population
#' @description Implements [F_fqM] for the GeRM model.
#' @inheritParams F_fqM
#' @return a [numeric] vector of length `nPatches`
#' @export
F_fqM.GeRM <- function(t, y, pars, s) {
  with(pars$MYZpar[[s]], f*q)*y[pars$ix$MYZ[[s]]$M_ix]
}


#' @title Number of eggs laid by adult mosquitoes
#' @description Implements [F_eggs] for the GeRM model.
#' @inheritParams F_eggs
#' @return a [numeric] vector of length `nPatches`
#' @export
F_eggs.GeRM <- function(t, y, pars, s) {
  M <- y[pars$ix$MYZ[[s]]$M_ix]
  with(pars$MYZpar[[s]], {
    return(M*nu*eggsPerBatch)
  })
}

#' @title Derivatives for adult mosquitoes
#' @description Implements [dMYZdt] for the GeRM ODE model.
#' @inheritParams dMYZdt
#' @return a [numeric] vector
#' @export
dMYZdt.GeRM_ode <- function(t, y, pars, s) {

  Lambda = pars$Lambda[[s]]
  kappa = pars$kappa[[s]]

  with(pars$ix$MYZ[[s]],{
    M <- y[M_ix]
    G <- y[G_ix]
    Y <- y[Y_ix]
    Z <- y[Z_ix]

    with(pars$MYZpar[[s]],{

      Omega <- make_Omega(g, sigma, calK, nPatches)
      Upsilon <- expm::expm(-Omega*eip)

      dMdt <- Lambda - (Omega %*% M)
      dGdt <-  f*(M - G) - nu*G  - (Omega %*% G)
      dYdt <- f*q*kappa*(M - Y) - (Omega %*% Y)
      dZdt <- Upsilon %*% diag(f*q*kappa, nPatches) %*% (M - Y) - (Omega %*% Z)

      return(c(dMdt, dGdt, dYdt, dZdt))
    })
  })
}

#' @title Derivatives for adult mosquitoes
#' @description Implements [dMYZdt] for the GeRM DDE model.
#' @inheritParams dMYZdt
#' @return a [numeric] vector
#' @importFrom deSolve lagvalue
#' @importFrom deSolve lagderiv
#' @export
dMYZdt.GeRM_dde <- function(t, y, pars, s){

  Lambda = pars$Lambda[[s]]
  kappa = pars$kappa[[s]]

  with(pars$ix$MYZ[[s]],{
    M <- y[M_ix]
    G <- y[G_ix]
    Y <- y[Y_ix]
    Z <- y[Z_ix]
    Upsilon <- matrix(data = y[Upsilon_ix], nrow = nPatches, ncol = nPatches)

    with(pars$MYZpar[[s]],{

      if (t < eip) {
        M_eip <- MYZinits[[s]]$M0
        Y_eip <- MYZinits[[s]]$Y0
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
      dGdt <-  f*(M - G) - nu*G  - (Omega %*% G)
      dYdt <- f*q*kappa*(M - Y) - (Omega %*% Y)
      dZdt <- Upsilon %*% (fqkappa_eip * (M_eip - Y_eip)) - (Omega %*% Z)
      dUdt <- as.vector(((1-dEIPdt(t,EIPmod))*Omega_eip - Omega) %*% Upsilon)

      return(c(dMdt, dGdt, dYdt, dZdt, dUdt, f*q*kappa, g, sigma))
    })
  })
}

#' @title Setup MYZpar for the GeRM model
#' @description Implements [setup_MYZpar] for the GeRM model
#' @inheritParams setup_MYZpar
#' @return a [list] vector
#' @export
setup_MYZpar.GeRM = function(MYZname, pars, s, MYZopts=list(), EIPmod, calK){
  pars$MYZpar[[s]] = make_MYZpar_GeRM(pars$nPatches, MYZopts, EIPmod, calK)
  return(pars)
}

#' @title Make parameters for GeRM ODE adult mosquito model
#' @param nPatches is the number of patches, an integer
#' @param MYZopts a [list] of values that overwrites the defaults
#' @param EIPmod a [list] that defines the EIP model
#' @param calK a mosquito dispersal matrix of dimensions `nPatches` by `nPatches`
#' @param g mosquito mortality rate
#' @param g_par parameters to configure F_g
#' @param sigma emigration rate
#' @param sigma_par parameters to configure F_sigma
#' @param f feeding rate
#' @param f_par parameters to configure F_f
#' @param q human blood fraction
#' @param q_par parameters to configure F_q
#' @param nu oviposition rate, per mosquito
#' @param nu_par parameters to configure F_nu
#' @param eggsPerBatch eggs laid per oviposition
#' @return a [list]
#' @export
make_MYZpar_GeRM = function(nPatches, MYZopts=list(), EIPmod, calK,
                            g=1/12, g_par = list(),
                            sigma=1/8, sigma_par = list(),
                            f=0.3, f_par = list(),
                            q=0.95, q_par = list(),
                            nu=1, nu_par = list(),
                            eggsPerBatch=60){

  stopifnot(is.matrix(calK))
  stopifnot(dim(calK) == c(nPatches, nPatches))

  with(MYZopts,{
    MYZpar <- list()

    MYZpar$xde <- solve_as
    class(MYZpar$xde) <- solve_as
    if(solve_as == 'dde') class(MYZpar) <- c('GeRM', 'GeRM_dde')
    if(solve_as == 'ode') class(MYZpar) <- c('GeRM', 'GeRM_ode')

    MYZpar$nPatches <- nPatches
    with(MYZopts,{
      MYZpar <- list()

      MYZpar$xde <- solve_as
      class(MYZpar$xde) <- solve_as
      if(solve_as == 'dde') class(MYZpar) <- c('GeRM', 'GeRM_dde')
      if(solve_as == 'ode') class(MYZpar) <- c('GeRM', 'GeRM_ode')

      MYZpar$g <- checkIt(g, pars$nPatches)
      MYZpar$g0 <- MYZpar$g
      if(length(g_par) == 0){
        MYZpar$g_par <- list()
        class(MYZpar$g_par) <- "static"
      } else MYZpar$g_par <- g_par

      MYZpar$sigma <- checkIt(sigma, pars$nPatches)
      MYZpar$sigma0 <- MYZpar$sigma
      if(length(sigma_par) == 0){
        MYZpar$sigma_par <- list()
        class(MYZpar$sigma_par) <- "static"
      } else MYZpar$sigma_par <- sigma_par

      MYZpar$f <- checkIt(f, pars$nPatches)
      MYZpar$f0 <- MYZpar$f
      if(length(f_par) == 0){
        MYZpar$f_par <- list()
        class(MYZpar$f_par) <- "static"
      } else MYZpar$f_par = f_par

      MYZpar$q <- checkIt(q, pars$nPatches)
      MYZpar$q0 <- MYZpar$q
      if(length(q_par) == 0){
        MYZpar$q_par <- list()
        class(MYZpar$q_par) <- "static"
      } else MYZpar$q_par <- q_par

      MYZpar$nu <- checkIt(nu, pars$nPatches)
      MYZpar$nu0 <- MYZpar$nu
      if(length(nu_par) == 0){
        MYZpar$nu_par <- list()
        class(MYZpar$nu_par) <- "static"
      } else MYZpar$nu_par=nu_par


    # The EIP model and the eip
    MYZpar$EIPmod <- EIPmod
    MYZpar$eip <- EIPmod(0, EIPmod)

    MYZpar$calK <- calK

    MYZpar$Omega <- make_Omega(g, sigma, calK, nPatches)
    MYZpar$Upsilon <- with(MYZpar, expm::expm(-Omega*eip))

    return(MYZpar)
    })
  })
}

#' @title Setup the GeRM model
#' @description Implements [setup_MYZinits] for the GeRM model
#' @inheritParams setup_MYZinits
#' @return a [list] vector
#' @export
setup_MYZinits.GeRM = function(pars, s, MYZopts=list()){
  pars$MYZinits[[s]] = with(pars$MYZpar[[s]], make_MYZinits_GeRM_dde(nPatches, Upsilon, MYZopts))
}


#' @title Make inits for GeRM adult mosquito model
#' @param nPatches the number of patches in the model
#' @param Upsilon a matrix describing survival and dispersal through the EIP
#' @param MYZopts a [list] of values that overwrites the defaults
#' @param M0 total mosquito density at each patch
#' @param G0 total parous mosquito density at each patch
#' @param Y0 infected mosquito density at each patch
#' @param Z0 infectious mosquito density at each patch
#' @return a [list]
#' @export
make_MYZinits_GeRM_dde = function(nPatches, Upsilon, MYZopts = list(),
                                M0=5, G0=1, Y0=1, Z0=1){
  with(MYZopts,{
    M = checkIt(M0, nPatches)
    G = checkIt(G0, nPatches)
    Y = checkIt(Y0, nPatches)
    Z = checkIt(Z0, nPatches)
    dd = rep(0, nPatches)
    return(list(M=M, G=G, Y=Y, Z=Z, Upsilon=as.vector(Upsilon), d1=dd, d2=dd, d3=dd))
  })
}

#' @title Make inits for GeRM adult mosquito model
#' @param nPatches the number of patches in the model
#' @param MYZopts a [list] of values that overwrites the defaults
#' @param M0 total mosquito density at each patch
#' @param G0 total parous mosquito density at each patch
#' @param Y0 infected mosquito density at each patch
#' @param Z0 infectious mosquito density at each patch
#' @return a [list]
#' @export
make_MYZinits_GeRM_ode = function(nPatches, MYZopts = list(),
                                M0=5, G0=1, Y0=1, Z0=1){
  with(MYZopts,{
    M = checkIt(M0, nPatches)
    G = checkIt(G0, nPatches)
    Y = checkIt(Y0, nPatches)
    Z = checkIt(Z0, nPatches)
    return(list(M=M, G=G, Y=Y, Z=Z))
  })
}

#' @title Add indices for adult mosquitoes to parameter list
#' @description Implements [make_indices_MYZ] for the GeRM model.
#' @inheritParams make_indices_MYZ
#' @return none
#' @importFrom utils tail
#' @export
make_indices_MYZ.GeRM_dde <- function(pars, s) {with(pars,{
  M_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(M_ix, 1)

  G_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(G_ix, 1)

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
  pars$ix$MYZ[[s]] = list(M_ix=M_ix, G_ix=G_ix, Y_ix=Y_ix, Z_ix=Z_ix,
              Upsilon_ix = Upsilon_ix, fqkappa_ix=fqkappa_ix,
              g_ix=g_ix, sigma_ix=sigma_ix, max_ix=max_ix)
  return(pars)
})}


#' @title Add indices for adult mosquitoes to parameter list
#' @description Implements [make_indices_MYZ] for the GeRM model.
#' @inheritParams make_indices_MYZ
#' @return none
#' @importFrom utils tail
#' @export
make_indices_MYZ.GeRM_ode <- function(pars, s) {with(pars,{

  M_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(M_ix, 1)

  G_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(G_ix, 1)

  Y_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(Y_ix, 1)

  Z_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(Z_ix, 1)

  pars$max_ix = max_ix
  pars$ix$MYZ[[s]] = list(M_ix=M_ix, G_ix=G_ix, Y_ix=Y_ix, Z_ix=Z_ix)
  return(pars)
})}

#' @title Make parameters for GeRM ODE adult mosquito model
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
make_parameters_MYZ_GeRM <- function(pars, g, sigma, f, q, nu, eggsPerBatch, eip, calK, solve_as = 'dde') {
  stopifnot(is.numeric(g), is.numeric(sigma), is.numeric(f), is.numeric(q), is.numeric(nu), is.numeric(eggsPerBatch))

  MYZpar <- list()
  MYZpar$xde = solve_as
  class(MYZpar$xde) <- solve_as

  if(solve_as == 'dde') class(MYZpar) <- c('GeRM', 'GeRM_dde')
  if(solve_as == 'ode') class(MYZpar) <- c('GeRM', 'GeRM_ode')

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

  EIPmod <- setup_eip_static(eip=eip)
  MYZpar$EIPmod <- EIPmod
  MYZpar$eip <- eip
  MYZpar$calK <- calK

  MYZpar <- MYZpar
  pars$MYZpar = list()
  pars$MYZpar[[1]] = MYZpar

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
  pars$MYZinits = list()
  MYZinits = list(M=M0, G=G0, Y=Y0, Z=Z0)
  pars$MYZinits[[1]] = MYZinits
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
  pars$MYZinits = list()
  dmy = rep(0, pars$nPatches)
  MYZinits = list(M=M0, G=G0, Y=Y0, Z=Z0, Upsilon=Upsilon0, d1=dmy, d2=dmy, d3=dmy)
  pars$MYZinits[[1]] = MYZinits
  return(pars)
}

#' @title Parse the output of deSolve and return variables for the GeRM model
#' @description Implements [parse_deout_MYZ] for the GeRM model
#' @inheritParams parse_deout_MYZ
#' @return a [list]
#' @export
parse_deout_MYZ.GeRM <- function(deout, pars, s) {with(pars$ix$MYZ[[s]],{
  time = deout[,1]
  M = deout[,M_ix+1]
  G = deout[,G_ix+1]
  Y = deout[,Y_ix+1]
  Z = deout[,Z_ix+1]
  y = Y/M
  z = Z/M
  return(list(time=time, M=M, G=G, Y=Y, Z=Z, y=y, z=z))
})}

#' @title Return initial values as a vector
#' @description Implements [get_inits_MYZ] for the GeRM model.
#' @inheritParams get_inits_MYZ
#' @return none
#' @export
get_inits_MYZ.GeRM_ode <- function(pars, s) {with(pars$MYZinits[[s]],{
  c(M, G, Y, Z)
})}

#' @title Return initial values as a vector
#' @description Implements [get_inits_MYZ] for the GeRM model.
#' @inheritParams get_inits_MYZ
#' @return [numeric]
#' @export
get_inits_MYZ.GeRM_dde <- function(pars, s) {with(pars$MYZinits[[s]],{
  c(M, G, Y, Z, as.vector(Upsilon), d1, d2, d3)
})}


#' @title Update the initial values for GeRM adult mosquito model
#' @inheritParams update_inits_MYZ
#' @return none
#' @export
update_inits_MYZ.GeRM_dde <- function(pars, y0, s) {with(pars$ix$MYZ[[s]],{
  M = y0[M_ix]
  G = y0[G_ix]
  Y = y0[Y_ix]
  Z = y0[Z_ix]
  Upsilon = y0[Upsilon_ix]
  pars$MYZinits[[s]] =  make_MYZinits_GeRM_dde(pars, Upsilon, list(), M0=M, G0=G, Y0=Y, Z0=Z )
  return(pars)
})}

#' @title Make inits for GeRM adult mosquito model
#' @inheritParams update_inits_MYZ
#' @return none
#' @export
update_inits_MYZ.GeRM_ode <- function(pars, y0, s) {with(pars$ix$MYZ[[s]],{
  M = y0[M_ix]
  G = y0[G_ix]
  Y = y0[Y_ix]
  Z = y0[Z_ix]
  pars$MYZinits[[s]] = make_MYZinits_GeRM_ode(pars, list(), M0=M, G0=G, Y0=Y, Z0=Z)
  return(pars)
})}


