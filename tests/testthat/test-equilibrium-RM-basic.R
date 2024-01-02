library(expm)
library(MASS)
library(deSolve)


test_that("test equilibrium with RM adults (ODE), basic competition", {

  numeric_tol <- 1e-5

  # set number of patches and strata
  nPatches <- 2
  nStrata <- 2
  nHabitats <- 2

  # parameters
  b <- 0.55
  c <- 0.15
  r <- 1/200
  wf <- rep(1, nStrata)

  f <- 0.3
  q <- 0.9
  g <- 1/10
  sigma <- 1/100
  nu <- 1/2
  eggsPerBatch <- 30

  eip <- 11

  # mosquito movement calK
  calK <- matrix(0, nPatches, nPatches)
  calK[upper.tri(calK)] <- rexp(sum(1:(nPatches-1)))
  calK[lower.tri(calK)] <- rexp(sum(1:(nPatches-1)))
  calK <- calK/rowSums(calK)
  calK <- t(calK)

  # omega matrix
  Omega <- make_Omega(g, sigma, calK, nPatches)
  Omega_inv <- solve(Omega)
  OmegaEIP <- expm::expm(-Omega * eip)
  OmegaEIP_inv <- expm::expm(Omega * eip)

  # human PfPR and H
  pfpr <- rep(0.3, times = nStrata)
  H <- rpois(n = nStrata, lambda = 1000)
  I <- rbinom(n = nStrata, size = H, prob = pfpr)

  # TaR
  Psi <- matrix(rexp(n = nStrata*nPatches), nStrata, nPatches)
  Psi <- Psi/rowSums(Psi)
  Psi <- t(Psi)

  # derived EIR to sustain equilibrium pfpr
  EIR <- diag(1/b, nPatches, nPatches) %*% ((r*I) / (H - I))

  # ambient pop
  W <- Psi %*% H

  # biting distribution matrix
  beta <- diag(wf) %*% t(Psi) %*% diag(1/as.vector(W), nPatches , nPatches)

  # kappa
  kappa <- t(beta) %*% (I*c)

  # equilibrium solutions for adults
  Z <- diag(1/(f*q), nPatches, nPatches) %*% ginv(beta) %*% EIR
  MY <- diag(1/as.vector(f*q*kappa), nPatches, nPatches) %*% OmegaEIP_inv %*% Omega %*% Z
  Y <- Omega_inv %*% (diag(as.vector(f*q*kappa), nPatches, nPatches) %*% MY)
  M <- MY + Y
  P <- solve(diag(f, nPatches) + Omega) %*% diag(f, nPatches) %*% M
  Lambda <- Omega %*% M

  # equilibrium solutions for aquatic
  calN <- matrix(0, nPatches, nHabitats)
  diag(calN) <- 1

  calU <- matrix(0, nHabitats, nPatches)
  diag(calU) <- 1

  alpha <- as.vector(solve(calN) %*% Lambda)

  psi <- 1/10
  phi <- 1/12
  eta <- as.vector(calU %*% M * nu * eggsPerBatch)

  L <- alpha/psi
  theta <- (eta - psi*L - phi*L)/(L^2)

  # parameters for exDE
  params <- make_parameters_xde()
  params$nStrata <- nStrata
  params$nPatches <- nPatches
  params$nHabitats <- nHabitats
  params$egg_laying[[1]] = list()
  params$egg_laying[[1]]$calU <- calU
  params$calN <- calN

  params = make_parameters_MYZ_RM(pars = params, g = g, sigma = sigma, calK = calK, eip = eip, f = f, q = q, nu = nu, eggsPerBatch = eggsPerBatch, solve_as = "ode")
  params = make_inits_MYZ_RM_ode(pars = params, M0 = as.vector(M), P0 = as.vector(P), Y0 = as.vector(Y), Z0 = as.vector(Z))
  params = make_parameters_L_basic(pars = params, psi = psi, phi = phi, theta = theta)
  params = make_inits_L_basic(pars = params,  L0 = L)

  params = make_indices(params)
  params$kappa = as.vector(kappa)

  # set initial conditions
  y0 <- get_inits(params)


  # run simulation
  out <- deSolve::ode(y = y0, times = c(0,190), func = xDE_diffeqn_mosy, parms = params, method = "lsoda")

  expect_equal(as.vector(out[2, params$ix$L[[1]]$L_ix+1]), as.vector(L), tolerance = numeric_tol)
  expect_equal(as.vector(out[2, params$ix$MYZ[[1]]$M_ix+1]), as.vector(M), tolerance = numeric_tol)
  expect_equal(as.vector(out[2, params$ix$MYZ[[1]]$P_ix+1]), as.vector(P), tolerance = numeric_tol)
  expect_equal(as.vector(out[2, params$ix$MYZ[[1]]$Y_ix+1]), as.vector(Y), tolerance = numeric_tol)
  expect_equal(as.vector(out[2, params$ix$MYZ[[1]]$Z_ix+1]), as.vector(Z), tolerance = numeric_tol)
})

test_that("test equilibrium with RM adults (DDE), basic competition", {

  numeric_tol <- 1e-5

  # set number of patches and strata
  nPatches <- 2
  nStrata <- 2
  nHabitats <- 2

  # parameters
  b <- 0.55
  c <- 0.15
  r <- 1/200
  wf <- rep(1, nStrata)

  f <- rep(0.3, nPatches)
  q <- rep(0.9, nPatches)
  g <- rep(1/10, nPatches)
  sigma <- rep(1/100, nPatches)
  nu <- rep(1/2, nPatches)
  eggsPerBatch <- 30

  eip <- 11

  # mosquito movement calK
  calK <- matrix(0, nPatches, nPatches)
  calK[upper.tri(calK)] <- rexp(sum(1:(nPatches-1)))
  calK[lower.tri(calK)] <- rexp(sum(1:(nPatches-1)))
  calK <- calK/rowSums(calK)
  calK <- t(calK)

  # omega matrix
  Omega <- make_Omega(g, sigma, calK, nPatches)
  Omega_inv <- solve(Omega)
  OmegaEIP <- expm::expm(-Omega * eip)
  OmegaEIP_inv <- expm::expm(Omega * eip)

  # human PfPR and H
  pfpr <- rep(0.3, times = nStrata)
  H <- rpois(n = nStrata, lambda = 1000)
  I <- rbinom(n = nStrata, size = H, prob = pfpr)

  # TaR
  Psi <- matrix(rexp(n = nStrata*nPatches), nStrata, nPatches)
  Psi <- Psi/rowSums(Psi)
  Psi <- t(Psi)

  # derived EIR to sustain equilibrium pfpr
  EIR <- diag(1/b, nPatches, nPatches) %*% ((r*I) / (H - I))

  # ambient pop
  W <- Psi %*% H

  # biting distribution matrix
  beta <- diag(wf) %*% t(Psi) %*% diag(1/as.vector(W), nPatches , nPatches)

  # kappa
  kappa <- t(beta) %*% (I*c)

  # equilibrium solutions for adults
  Z <- diag(1/(f*q), nPatches, nPatches) %*% ginv(beta) %*% EIR
  MY <- diag(1/as.vector(f*q*kappa), nPatches, nPatches) %*% OmegaEIP_inv %*% Omega %*% Z
  Y <- Omega_inv %*% (diag(as.vector(f*q*kappa), nPatches, nPatches) %*% MY)
  M <- MY + Y
  P <- solve(diag(f, nPatches) + Omega) %*% diag(f, nPatches) %*% M
  Lambda <- Omega %*% M

  # equilibrium solutions for aquatic
  calN <- matrix(0, nPatches, nHabitats)
  diag(calN) <- 1

  calU <- matrix(0, nHabitats, nPatches)
  diag(calU) <- 1

  alpha <- as.vector(solve(calN) %*% Lambda)

  psi <- 1/10
  phi <- 1/12
  eta <- as.vector(calU %*% M * nu * eggsPerBatch)

  L <- alpha/psi
  theta <- (eta - psi*L - phi*L)/(L^2)

  params <- make_parameters_xde()
  params$nStrata <- nStrata
  params$nPatches <- nPatches
  params$nHabitats <- nHabitats
  params$egg_laying[[1]] = list()
  params$egg_laying[[1]]$calU <- calU
  params$calN <- calN

  # parameters for exDE
  params = make_parameters_MYZ_RM(pars = params, g = g, sigma = sigma, calK = calK, eip = eip, f = f, q = q, nu = nu, eggsPerBatch = eggsPerBatch)
  params = make_inits_MYZ_RM_dde(pars = params, M0 = as.vector(M), P0 = as.vector(P), Y0 = as.vector(Y), Z0 = as.vector(Z), Upsilon0=as.vector(OmegaEIP))
  params = make_parameters_L_basic(pars = params, psi = psi, phi = phi, theta = theta)
  params = make_inits_L_basic(pars = params,  L0 = L)

  params = make_indices(params)

  # set initial conditions
  y0 <- get_inits(params)
  params$kappa = as.vector(kappa)


  # run simulation
  out <- deSolve::dede(y = y0, times = c(0,90), func = xDE_diffeqn_mosy, parms = params, method = "lsoda")

  expect_equal(as.vector(out[2, params$ix$L[[1]]$L_ix+1]), as.vector(L), tolerance = numeric_tol)
  expect_equal(as.vector(out[2, params$ix$MYZ[[1]]$M_ix+1]), as.vector(M), tolerance = numeric_tol)
  expect_equal(as.vector(out[2, params$ix$MYZ[[1]]$P_ix+1]), as.vector(P), tolerance = numeric_tol)
  expect_equal(as.vector(out[2, params$ix$MYZ[[1]]$Y_ix+1]), as.vector(Y), tolerance = numeric_tol)
  expect_equal(as.vector(out[2, params$ix$MYZ[[1]]$Z_ix+1]), as.vector(Z), tolerance = numeric_tol)
})
