library(expm)
library(MASS)
library(deSolve)

test_that("test equilibrium with RM adults (ODE), basic competition", {

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

  tau <- 11

  # mosquito movement calK
  K <- matrix(0, nPatches, nPatches)
  K[upper.tri(K)] <- rexp(sum(1:(nPatches-1)))
  K[lower.tri(K)] <- rexp(sum(1:(nPatches-1)))
  K <- K/rowSums(K)
  K <- t(K)

  # omega matrix
  Omega <- diag(g, nPatches) + (diag(sigma, nPatches) %*% (diag(nPatches) - K))
  Omega_inv <- solve(Omega)
  OmegaEIP <- expm::expm(-Omega * tau)
  OmegaEIP_inv <- expm::expm(Omega * tau)

  # human PfPR and H
  pfpr <- rep(0.3, times = nStrata)
  H <- rpois(n = nStrata, lambda = 1000)
  X <- rbinom(n = nStrata, size = H, prob = pfpr)

  # TaR
  Psi <- matrix(rexp(n = nStrata*nPatches), nStrata, nPatches)
  Psi <- Psi/rowSums(Psi)
  Psi <- t(Psi)

  # derived EIR to sustain equilibrium pfpr
  EIR <- diag(1/b, nPatches, nPatches) %*% ((r*X) / (H - X))

  # ambient pop
  W <- Psi %*% H

  # biting distribution matrix
  beta <- diag(wf) %*% t(Psi) %*% diag(1/as.vector(W), nPatches , nPatches)

  # kappa
  kappa <- t(beta) %*% (X*c)

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

  alpha <- as.vector(solve(calN) %*% Lambda)

  psi <- 1/10
  phi <- 1/12
  eta <- as.vector(calU %*% G * nu * eggsPerBatch)

  L <- alpha/psi
  theta <- (eta - psi*L - phi*L)/(L^2)

  # parameters for exDE
  params <- new.env()
  params$nStrata <- nStrata
  params$nPatches <- nPatches
  params$nHabitats <- nHabitats
  params$calU <- calU
  params$calN <- calN
  params$beta <- beta
  params$betaT <- t(beta)

  params$MYZpar <- make_parameters_MYZ_RM_ode(Omega = Omega, OmegaEIP = OmegaEIP, f = f, q = q, nu = nu, eggsPerBatch = eggsPerBatch, M0 = as.vector(M), G0 = as.vector(G), Y0 = as.vector(Y), Z0 = as.vector(Z))
  params$Lpar <- Lpar <- make_parameters_L_basic(psi = psi, phi = phi, theta = theta, L0 = L)

  params <- make_indices(params)

  # set initial conditions
  y <- rep(NaN, max(params$Z_ix))
  y[params$L_ix] <- as.vector(L)
  y[params$M_ix] <- as.vector(M)
  y[params$G_ix] <- as.vector(G)
  y[params$Y_ix] <- as.vector(Y)
  y[params$Z_ix] <- as.vector(Z)

  # run simulation
  out <- deSolve::ode(y = y, times = c(0,50), func = xDE_diffeqn_mosy, parms = params, method = "lsoda", kappa = as.vector(kappa))

  expect_equal(as.vector(out[2, params$L_ix+1]), as.vector(L), tolerance = 1e-4)
  expect_equal(as.vector(out[2, params$M_ix+1]), as.vector(M))
  expect_equal(as.vector(out[2, params$G_ix+1]), as.vector(G))
  expect_equal(as.vector(out[2, params$Y_ix+1]), as.vector(Y))
  expect_equal(as.vector(out[2, params$Z_ix+1]), as.vector(Z))
})

test_that("test equilibrium with RM adults (ODE), basic competition", {

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

  tau <- 11

  # mosquito movement calK
  K <- matrix(0, nPatches, nPatches)
  K[upper.tri(K)] <- rexp(sum(1:(nPatches-1)))
  K[lower.tri(K)] <- rexp(sum(1:(nPatches-1)))
  K <- K/rowSums(K)
  K <- t(K)

  # omega matrix
  Omega <- diag(g, nPatches) + (diag(sigma, nPatches) %*% (diag(nPatches) - K))
  Omega_inv <- solve(Omega)
  OmegaEIP <- expm::expm(-Omega * tau)
  OmegaEIP_inv <- expm::expm(Omega * tau)

  # human PfPR and H
  pfpr <- rep(0.3, times = nStrata)
  H <- rpois(n = nStrata, lambda = 1000)
  X <- rbinom(n = nStrata, size = H, prob = pfpr)

  # TaR
  Psi <- matrix(rexp(n = nStrata*nPatches), nStrata, nPatches)
  Psi <- Psi/rowSums(Psi)
  Psi <- t(Psi)

  # derived EIR to sustain equilibrium pfpr
  EIR <- diag(1/b, nPatches, nPatches) %*% ((r*X) / (H - X))

  # ambient pop
  W <- Psi %*% H

  # biting distribution matrix
  beta <- diag(wf) %*% t(Psi) %*% diag(1/as.vector(W), nPatches , nPatches)

  # kappa
  kappa <- t(beta) %*% (X*c)

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

  alpha <- as.vector(solve(calN) %*% Lambda)

  psi <- 1/10
  phi <- 1/12
  eta <- as.vector(calU %*% G * nu * eggsPerBatch)

  L <- alpha/psi
  theta <- (eta - psi*L - phi*L)/(L^2)

  # parameters for exDE
  params <- new.env()
  params$nStrata <- nStrata
  params$nPatches <- nPatches
  params$nHabitats <- nHabitats
  params$calU <- calU
  params$calN <- calN
  params$beta <- beta
  params$betaT <- t(beta)

  params$MYZpar <- make_parameters_MYZ_RM_dde(Omega = Omega, OmegaEIP = OmegaEIP, f = f, q = q, nu = nu, eggsPerBatch = eggsPerBatch, tau = tau, M0 = as.vector(M), G0 = as.vector(G), Y0 = as.vector(Y), Z0 = as.vector(Z))
  params$Lpar <- Lpar <- make_parameters_L_basic(psi = psi, phi = phi, theta = theta, L0 = L)

  params <- make_indices(params)

  # set initial conditions
  y <- rep(NaN, max(params$Z_ix))
  y[params$L_ix] <- as.vector(L)
  y[params$M_ix] <- as.vector(M)
  y[params$G_ix] <- as.vector(G)
  y[params$Y_ix] <- as.vector(Y)
  y[params$Z_ix] <- as.vector(Z)

  # run simulation
  out <- deSolve::dede(y = y, times = c(0,50), func = xDE_diffeqn_mosy, parms = params, method = "lsoda", kappa = t(cbind(kappa,kappa)))

  expect_equal(as.vector(out[2, params$L_ix+1]), as.vector(L), tolerance = 1e-4)
  expect_equal(as.vector(out[2, params$M_ix+1]), as.vector(M))
  expect_equal(as.vector(out[2, params$G_ix+1]), as.vector(G))
  expect_equal(as.vector(out[2, params$Y_ix+1]), as.vector(Y))
  expect_equal(as.vector(out[2, params$Z_ix+1]), as.vector(Z))
})

