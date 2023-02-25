library(expm)
library(deSolve)
library(MASS)

numeric_tol <- 1e-5

test_that("forced emergence works with equilibrium", {
  nPatches <- 3
  nHabitats <- 4
  f <- 0.3
  q <- 0.9
  g <- 1/20
  sigma <- 1/10
  tau <- 11
  nu <- 1/2
  eggsPerBatch <- 30

  calN <- matrix(0, nPatches, nHabitats)
  calN[1,1] <- 1
  calN[2,2] <- 1
  calN[3,3] <- 1
  calN[3,4] <- 1

  calU <- matrix(0, nHabitats, nPatches)
  calU[1,1] <- 1
  calU[2,2] <- 1
  calU[3:4,3] <- 0.5

  calK <- matrix(0, nPatches, nPatches)
  calK[1, 2:3] <- c(0.2, 0.8)
  calK[2, c(1,3)] <- c(0.5, 0.5)
  calK[3, 1:2] <- c(0.7, 0.3)
  calK <- t(calK)

  Omega <- make_Omega(g, sigma, calK, nPatches)
  OmegaEIP <- expm::expm(-Omega * tau)

  kappa <- c(0.1, 0.075, 0.025)
  Lambda <- c(5, 10, 8)

  # equilibrium solutions (forward)
  Omega_inv <- solve(Omega)
  OmegaEIP_inv <- expm::expm(Omega * tau)

  M_eq <- as.vector(Omega_inv %*% Lambda)
  G_eq <- as.vector(solve(diag(nu+f, nPatches) + Omega) %*% diag(f, nPatches) %*% M_eq)
  Y_eq <- as.vector(solve(diag(f*q*kappa) + Omega) %*% diag(f*q*kappa) %*% M_eq)
  Z_eq <- as.vector(Omega_inv %*% OmegaEIP %*% diag(f*q*kappa) %*% (M_eq - Y_eq))

  # the "Lambda" for the dLdt model
  alpha <- as.vector(ginv(calN) %*% Lambda)

  params <- make_parameters_xde()
  params$nPatches = nPatches
  params$nHabitats = nHabitats
  params$calU = calU
  params$calN = calN


  # ODE
  params = make_parameters_MYZ_GeRM(pars = params, g = g, sigma = sigma, calK = calK, tau = tau, f = f, q = q, nu = nu, eggsPerBatch = eggsPerBatch, solve_as = "ode")
  params = make_inits_MYZ_GeRM(pars = params, M0 = rep(0, nPatches), G0 = rep(0, nPatches), Y0 = rep(0, nPatches), Z0 = rep(0, nPatches), Upsilon0 = OmegaEIP)
  params = make_parameters_L_trace(pars = params, Lambda = alpha)
  params = make_parameters_vc_null(pars = params)
  params = make_parameters_exogenous_null(pars = params)

  params = make_indices(params)

  y0 <- get_inits(params)

  # mimic MosyBehavior
  MosyBehavior <- list()
  MosyBehavior$f <- rep(params$MYZpar$f, 2)
  attr(MosyBehavior$f, 'time') <- c(0, 0 - params$MYZpar$tau)
  MosyBehavior$q <- rep(params$MYZpar$q, 2)
  MosyBehavior$g <- rep(params$MYZpar$g, 2)

  out <- deSolve::ode(y = y0, times = c(0, 365), func = xDE_diffeqn_mosy, parms = params, method = 'lsoda', kappa = kappa, MosyBehavior = MosyBehavior)

  M_sim <- as.vector(out[2, params$M_ix+1])
  G_sim <- as.vector(out[2, params$G_ix+1])
  Y_sim <- as.vector(out[2, params$Y_ix+1])
  Z_sim <- as.vector(out[2, params$Z_ix+1])

  expect_equal(M_eq, M_sim, tolerance = numeric_tol)
  expect_equal(G_eq, G_sim, tolerance = numeric_tol)
  expect_equal(Y_eq, Y_sim, tolerance = numeric_tol)
  expect_equal(Z_eq, Z_sim, tolerance = numeric_tol)

})
