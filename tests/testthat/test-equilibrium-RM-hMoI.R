library(expm)
library(MASS)
library(deSolve)

numeric_tol <- 1e-5

test_that("test equilibrium with RM adults (ODE), hMoI humans, trace", {

  # set number of patches and strata
  nPatches <- 2
  nStrata <- 2
  nHabitats <- 2

  # parameters
  b <- 0.55
  c <- 0.15
  c1 <- 0.05
  c2 <- 0.25
  r1 <- 1/250
  r2 <- 1/50
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

  # MoI at equilibrium
  H <- c(100, 120)
  residence = c(1,2)
  searchWtsH = c(1,1)
  m20 <- 1.5
  h <- r2*m20
  m10 <- h/r1

  # TaR
  TaR <- matrix(rexp(n = nStrata*nPatches), nStrata, nPatches)
  TaR <- TaR/rowSums(TaR)
  TaR <- t(TaR)

  # derived EIR to sustain equilibrium pfpr
  EIR <- rep(h/b, nStrata)

  # ambient pop
  W <- TaR %*% H

  # biting distribution matrix
  beta <- diag(wf) %*% t(TaR) %*% diag(1/as.vector(W), nPatches , nPatches)

  # kappa
  # kappa <- t(beta) %*% (X*c)
  x1 <- pexp(q = rep(m10, nStrata))
  x2 <- pexp(q = rep(m20, nStrata))
  x <- ((c2 * x2) + (c1 * (x1 - x2))) * H
  kappa <- t(beta) %*% x

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

  # parameters for exDE
  params <- make_parameters_xde()
  params$nStrata <- nStrata
  params$nPatches <- nPatches
  params$nHabitats <- nHabitats
  params$egg_laying[[1]] = list()
  params$egg_laying[[1]]$calU <- calU
  params$calN <- calN

  params = make_parameters_MYZ_RM(pars = params, g = g, sigma = sigma, calK = calK, eip = eip, f = f, q = q, nu = nu, eggsPerBatch = eggsPerBatch, solve_as="ode")
  params = make_inits_MYZ_RM_ode(pars = params, M0 = as.vector(M), P0 = as.vector(P), Y0 = as.vector(Y), Z0 = as.vector(Z))
  params = make_parameters_demography_null(pars = params, H=H, residence=residence,
                                           searchWts=searchWtsH, TaR=TaR)
  params = make_parameters_X_hMoI(pars = params, b = b, c1 = c1, c2 = c2, r1 = r1, r2 = r2)
  params = make_inits_X_hMoI(pars = params, m10 = rep(m10,2), m20 = rep(m20,2))
  params = make_parameters_L_trace(pars = params, Lambda = as.vector(Lambda))

  params = make_indices(params)


  # set initial conditions
  y0 <- get_inits(params)

  # run simulation
  out <- deSolve::ode(y = y0, times = c(0,50), func = xDE_diffeqn, parms = params, method = "lsoda")

  expect_equal(as.vector(out[2, params$ix$MYZ[[1]]$M_ix+1]), as.vector(M), tolerance = numeric_tol)
  expect_equal(as.vector(out[2, params$ix$MYZ[[1]]$P_ix+1]), as.vector(P), tolerance = numeric_tol)
  expect_equal(as.vector(out[2, params$ix$MYZ[[1]]$Y_ix+1]), as.vector(Y), tolerance = numeric_tol)
  expect_equal(as.vector(out[2, params$ix$MYZ[[1]]$Z_ix+1]), as.vector(Z), tolerance = numeric_tol)
  expect_equal(as.vector(out[2L, params$ix$X[[1]]$m1_ix+1]), rep(m10, nStrata), tolerance = numeric_tol)
  expect_equal(as.vector(out[2L, params$ix$X[[1]]$m2_ix+1]), rep(m20, nStrata), tolerance = numeric_tol)
})

test_that("test equilibrium with RM adults (DDE), hMoI humans, trace", {

  # set number of patches and strata
  nPatches <- 2
  nStrata <- 2
  nHabitats <- 2

  # parameters
  b <- 0.55
  c <- 0.15
  c1 <- 0.05
  c2 <- 0.25
  r1 <- 1/250
  r2 <- 1/50
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

  # MoI at equilibrium
  H <- c(100, 120)
  residence = c(1,2)
  searchWtsH = c(1,1)
  m20 <- 1.5
  h <- r2*m20
  m10 <- h/r1

  # TaR
  TaR <- matrix(rexp(n = nStrata*nPatches), nStrata, nPatches)
  TaR <- TaR/rowSums(TaR)
  TaR <- t(TaR)

  # derived EIR to sustain equilibrium pfpr
  EIR <- rep(h/b, nStrata)

  # ambient pop
  W <- TaR %*% H

  # biting distribution matrix
  beta <- diag(wf) %*% t(TaR) %*% diag(1/as.vector(W), nPatches , nPatches)

  # kappa
  # kappa <- t(beta) %*% (X*c)
  x1 <- pexp(q = rep(m10, nStrata))
  x2 <- pexp(q = rep(m20, nStrata))
  x <- ((c2 * x2) + (c1 * (x1 - x2))) * H
  kappa <- t(beta) %*% x

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

  # parameters for exDE
  params <- make_parameters_xde()
  params$nStrata <- nStrata
  params$nPatches <- nPatches
  params$nHabitats <- nHabitats
  params$egg_laying[[1]] = list()
  params$egg_laying[[1]]$calU <- calU
  params$calN <- calN

  params = make_parameters_MYZ_RM(pars = params, g = g, sigma = sigma, calK = calK, eip = eip, f = f, q = q, nu = nu, eggsPerBatch = eggsPerBatch, solve_as="ode")
  params = make_inits_MYZ_RM_dde(pars = params, M0 = as.vector(M), P0 = as.vector(P), Y0 = as.vector(Y), Z0 = as.vector(Z), Upsilon0=OmegaEIP)
  params = make_parameters_demography_null(pars = params, H=H, residence=residence,
                                           searchWts=searchWtsH, TaR=TaR)
  params = make_parameters_X_hMoI(pars = params, b = b, c1 = c1, c2 = c2, r1 = r1, r2 = r2)
  params = make_inits_X_hMoI(pars = params, m10 = rep(m10,2), m20 = rep(m20,2))
  params = make_parameters_L_trace(pars = params, Lambda = as.vector(Lambda))

  params = make_indices(params)


  # set initial conditions
  y0 <- get_inits(params)

  # run simulation
  out <- deSolve::dede(y = y0, times = c(0,50), func = xDE_diffeqn, parms = params, method = "lsoda")

  expect_equal(as.vector(out[2, params$ix$MYZ[[1]]$M_ix+1]), as.vector(M), tolerance = numeric_tol)
  expect_equal(as.vector(out[2, params$ix$MYZ[[1]]$P_ix+1]), as.vector(P), tolerance = numeric_tol)
  expect_equal(as.vector(out[2, params$ix$MYZ[[1]]$Y_ix+1]), as.vector(Y), tolerance = numeric_tol)
  expect_equal(as.vector(out[2, params$ix$MYZ[[1]]$Z_ix+1]), as.vector(Z), tolerance = numeric_tol)
  expect_equal(as.vector(out[2L, params$ix$X[[1]]$m1_ix+1]), rep(m10, nStrata), tolerance = numeric_tol)
  expect_equal(as.vector(out[2L, params$ix$X[[1]]$m2_ix+1]), rep(m20, nStrata), tolerance = numeric_tol)
})
