library(MASS)
library(expm)
library(deSolve)

numeric_tol <- 1e-5

test_that("Le Menach VC model with 0 coverage stays roughly at equilibrium", {

  pars <- new.env()
  pars$nPatches <- 3
  pars$nStrata <- 3
  pars$nHabitats <- 3

  # parameters
  b <- 0.55
  c <- 0.15
  r <- 1/200
  wf <- rep(1, pars$nStrata)

  f <- 0.3
  q <- 0.9
  g <- 1/10
  sigma <- 1/100
  nu <- 1/2
  eggsPerBatch <- 30

  tau <- 11

  # mosquito movement calK
  calK <- matrix(0, pars$nPatches, pars$nPatches)
  calK[upper.tri(calK)] <- 1/(pars$nPatches-1)
  calK[lower.tri(calK)] <- 1/(pars$nPatches-1)
  calK <- calK/rowSums(calK)
  calK <- t(calK)

  Omega <- make_Omega(g = g, sigma = sigma, K = calK, nPatches = pars$nPatches)
  Omega_inv <- solve(Omega)
  Upsilon <- expm::expm(-Omega * tau)
  Upsilon_inv <- expm::expm(Omega * tau)

  # human PfPR and H
  pfpr <- runif(n = pars$nStrata, min = 0.25, max = 0.35)
  H <- rpois(n = pars$nStrata, lambda = 1000)
  X <- rbinom(n = pars$nStrata, size = H, prob = pfpr)

  Psi <- matrix(
    data = c(
      0.9, 0.05, 0.05,
      0.05, 0.9, 0.05,
      0.05, 0.05, 0.9
    ), nrow = pars$nStrata, ncol = pars$nPatches, byrow = T
  )
  Psi <- t(Psi)

  # derived EIR to sustain equilibrium pfpr
  EIR <- diag(1/b, pars$nStrata) %*% ((r*X) / (H - X))

  # ambient pop
  W <- Psi %*% H

  # biting distribution matrix
  beta <- diag(wf) %*% t(Psi) %*% diag(1/as.vector(W), pars$nPatches)

  # kappa
  kappa <- t(beta) %*% (X*c)

  # equilibrium solutions
  Z <- diag(1/(f*q), pars$nPatches) %*% ginv(beta) %*% EIR
  MY <- diag(1/as.vector(f*q*kappa), pars$nPatches) %*% Upsilon_inv %*% Omega %*% Z
  Y <- Omega_inv %*% (diag(as.vector(f*q*kappa), pars$nPatches) %*% MY)
  M <- MY + Y
  G <- solve(diag(nu+f, pars$nPatches) + Omega) %*% diag(f, pars$nPatches) %*% M
  Lambda <- Omega %*% M

  # set parameters
  make_parameters_exogenous_null(pars)
  make_parameters_vc_lemenach(pars)
  make_parameters_MYZ_GeRM_ode(pars = pars, g = g, sigma = sigma, calK = calK, f = f, q = q, nu = nu, eggsPerBatch = eggsPerBatch, tau = tau, M0 = M, G0 = G, Y0 = Y, Z0 = Z)
  make_parameters_L_trace(pars = pars, Lambda = as.vector(Lambda))
  make_parameters_X_SIS(pars = pars, b = b, c = c, r = r, Psi = Psi, wf = wf, X0 = X, H = H)
  make_indices(pars = pars)

  pars$calU <- diag(pars$nPatches)
  pars$calN <- diag(pars$nHabitats)

  # ICs
  y <- rep(NaN, pars$max_ix)
  y[pars$X_ix] <- pars$Xpar$X0
  y[pars$M_ix] <- pars$MYZpar$M0
  y[pars$G_ix] <- pars$MYZpar$G0
  y[pars$Y_ix] <- pars$MYZpar$Y0
  y[pars$Z_ix] <- pars$MYZpar$Z0
  y[pars$Upsilon_ix] <- as.vector(Upsilon)

  # solve the model
  out = ode(y = y, times = c(0, 365), func = xDE_diffeqn, parms = pars)

  # check it stays at equilibrium with phi = 0 for all time
  expect_equal(out[1, -1], out[2, -1], tolerance = numeric_tol)
})
