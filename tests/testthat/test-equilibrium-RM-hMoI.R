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

  tau <- 11

  # mosquito movement calK
  calK <- matrix(0, nPatches, nPatches)
  calK[upper.tri(calK)] <- rexp(sum(1:(nPatches-1)))
  calK[lower.tri(calK)] <- rexp(sum(1:(nPatches-1)))
  calK <- calK/rowSums(calK)
  calK <- t(calK)

  # omega matrix
  Omega <- make_Omega(g, sigma, calK, nPatches)
  Omega_inv <- solve(Omega)
  OmegaEIP <- expm::expm(-Omega * tau)
  OmegaEIP_inv <- expm::expm(Omega * tau)

  # MoI at equilibrium
  H <- c(100, 120)
  m20 <- 1.5
  h <- r2*m20
  m10 <- h/r1

  # TaR
  Psi <- matrix(rexp(n = nStrata*nPatches), nStrata, nPatches)
  Psi <- Psi/rowSums(Psi)
  Psi <- t(Psi)

  # derived EIR to sustain equilibrium pfpr
  EIR <- rep(h/b, nStrata)

  # ambient pop
  W <- Psi %*% H

  # biting distribution matrix
  beta <- diag(wf) %*% t(Psi) %*% diag(1/as.vector(W), nPatches , nPatches)

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
  G <- solve(diag(nu+f, nPatches) + Omega) %*% diag(f, nPatches) %*% M
  Lambda <- Omega %*% M

  # equilibrium solutions for aquatic
  calN <- matrix(0, nPatches, nHabitats)
  diag(calN) <- 1

  calU <- matrix(0, nHabitats, nPatches)
  diag(calU) <- 1

  # parameters for exDE
  params <- new.env()
  params$nStrata <- nStrata
  params$nPatches <- nPatches
  params$nHabitats <- nHabitats
  params$calU <- calU
  params$calN <- calN

  make_parameters_MYZ_GeRM_ode(pars = params, g = g, sigma = sigma, calK = calK, tau = tau, f = f, q = q, nu = nu, eggsPerBatch = eggsPerBatch, M0 = as.vector(M), G0 = as.vector(G), Y0 = as.vector(Y), Z0 = as.vector(Z))
  make_parameters_X_hMoI(pars = params, b = b, c1 = c1, c2 = c2, r1 = r1, r2 = r2, Psi = Psi, m10 = m10, m20 = m20)
  make_parameters_demography_null(pars = params, H = H)
  make_parameters_L_trace(pars = params, Lambda = as.vector(Lambda))
  make_parameters_exogenous_null(pars = params)
  make_parameters_vc_null(pars = params)

  make_indices(params)

  # set initial conditions
  y <- rep(NaN, params$max_ix)
  y[params$M_ix] <- as.vector(M)
  y[params$G_ix] <- as.vector(G)
  y[params$Y_ix] <- as.vector(Y)
  y[params$Z_ix] <- as.vector(Z)
  y[params$Upsilon_ix] <- as.vector(OmegaEIP)
  y[params$m1_ix] <- m10
  y[params$m2_ix] <- m20

  # run simulation
  out <- deSolve::ode(y = y, times = c(0,50), func = xDE_diffeqn, parms = params, method = "lsoda")

  expect_equal(as.vector(out[2, params$M_ix+1]), as.vector(M), tolerance = numeric_tol)
  expect_equal(as.vector(out[2, params$G_ix+1]), as.vector(G), tolerance = numeric_tol)
  expect_equal(as.vector(out[2, params$Y_ix+1]), as.vector(Y), tolerance = numeric_tol)
  expect_equal(as.vector(out[2, params$Z_ix+1]), as.vector(Z), tolerance = numeric_tol)
  expect_equal(as.vector(out[2L, params$m1_ix+1]), rep(m10, nStrata), tolerance = numeric_tol)
  expect_equal(as.vector(out[2L, params$m2_ix+1]), rep(m20, nStrata), tolerance = numeric_tol)
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

  f <- 0.3
  q <- 0.9
  g <- 1/10
  sigma <- 1/100
  nu <- 1/2
  eggsPerBatch <- 30

  tau <- 11

  # mosquito movement calK
  calK <- matrix(0, nPatches, nPatches)
  calK[upper.tri(calK)] <- rexp(sum(1:(nPatches-1)))
  calK[lower.tri(calK)] <- rexp(sum(1:(nPatches-1)))
  calK <- calK/rowSums(calK)
  calK <- t(calK)

  # omega matrix
  Omega <- make_Omega(g, sigma, calK, nPatches)
  Omega_inv <- solve(Omega)
  OmegaEIP <- expm::expm(-Omega * tau)
  OmegaEIP_inv <- expm::expm(Omega * tau)

  # MoI at equilibrium
  H <- c(100, 120)
  m20 <- 1.5
  h <- r2*m20
  m10 <- h/r1

  # TaR
  Psi <- matrix(rexp(n = nStrata*nPatches), nStrata, nPatches)
  Psi <- Psi/rowSums(Psi)
  Psi <- t(Psi)

  # derived EIR to sustain equilibrium pfpr
  EIR <- rep(h/b, nStrata)

  # ambient pop
  W <- Psi %*% H

  # biting distribution matrix
  beta <- diag(wf) %*% t(Psi) %*% diag(1/as.vector(W), nPatches , nPatches)

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
  G <- solve(diag(nu+f, nPatches) + Omega) %*% diag(f, nPatches) %*% M
  Lambda <- Omega %*% M

  # equilibrium solutions for aquatic
  calN <- matrix(0, nPatches, nHabitats)
  diag(calN) <- 1

  calU <- matrix(0, nHabitats, nPatches)
  diag(calU) <- 1

  # parameters for exDE
  params <- new.env()
  params$nStrata <- nStrata
  params$nPatches <- nPatches
  params$nHabitats <- nHabitats
  params$calU <- calU
  params$calN <- calN

  make_parameters_MYZ_GeRM_dde(pars = params, g = g, sigma = sigma, calK = calK, tau = tau, f = f, q = q, nu = nu, eggsPerBatch = eggsPerBatch, M0 = as.vector(M), G0 = as.vector(G), Y0 = as.vector(Y), Z0 = as.vector(Z))
  make_parameters_X_hMoI(pars = params, b = b, c1 = c1, c2 = c2, r1 = r1, r2 = r2, Psi = Psi, m10 = m10, m20 = m20)
  make_parameters_demography_null(pars = params, H = H)
  make_parameters_L_trace(pars = params, Lambda = as.vector(Lambda))
  make_parameters_exogenous_null(pars = params)
  make_parameters_vc_null(pars = params)

  make_indices(params)

  # set initial conditions
  y <- rep(NaN, params$max_ix)
  y[params$M_ix] <- as.vector(M)
  y[params$G_ix] <- as.vector(G)
  y[params$Y_ix] <- as.vector(Y)
  y[params$Z_ix] <- as.vector(Z)
  y[params$Upsilon_ix] <- as.vector(OmegaEIP)
  y[params$m1_ix] <- m10
  y[params$m2_ix] <- m20

  # run simulation
  out <- deSolve::dede(y = y, times = c(0,50), func = xDE_diffeqn, parms = params, method = "lsoda")

  expect_equal(as.vector(out[2, params$M_ix+1]), as.vector(M), tolerance = numeric_tol)
  expect_equal(as.vector(out[2, params$G_ix+1]), as.vector(G), tolerance = numeric_tol)
  expect_equal(as.vector(out[2, params$Y_ix+1]), as.vector(Y), tolerance = numeric_tol)
  expect_equal(as.vector(out[2, params$Z_ix+1]), as.vector(Z), tolerance = numeric_tol)
  expect_equal(as.vector(out[2L, params$m1_ix+1]), rep(m10, nStrata), tolerance = numeric_tol)
  expect_equal(as.vector(out[2L, params$m2_ix+1]), rep(m20, nStrata), tolerance = numeric_tol)
})
