# test that passing trace is the same as setting EIR=0 and solving dynamically
library(deSolve)
library(expm)

test_that("SIS model; demography with trace is the same as EIR=0 dynamic demography", {

  n <- 4
  b <- c(2e-4, 0.0005, 0.0005/2, 1e-6)
  d <- c(1e-4, 0.0001, 0.0001/2, 1e-6)
  m <- c(1/16, 1/30, 1/14)

  # solve with trace based demography
  calD <- make_calD(d,m,b)
  calD_eig <- eigen(calD)

  params <- new.env()
  params$nStrata <- n
  params$X_ix <- 1:n

  H0 <- rep(1e3,n)

  Y <- calD_eig$vectors
  Q <- diag(calD_eig$values)
  Yinv <- solve(Y)

  Ht <- function(t) {
    as.vector(Y %*% (expm::expm(Q * t) %*% (Yinv %*% H0)))
  }

  make_parameters_X_SIS(pars = params, b = 0.55, c = 0.15, r = 1/200, Psi = diag(n), wf = rep(1,n), X0 = rep(0,n))
  make_parameters_demography_trace(pars = params, H = Ht)

  y <- rep(0,n)

  out_trace <- deSolve::ode(y = y, times = 0:1000, func = function(t, y, pars, EIR) {
    list(dXdt(t = t, y = y, pars = pars, EIR = EIR))
  }, parms = params, EIR = rep(0, n))

  # solve with dynamic demography
  params <- new.env()
  params$nStrata <- n
  params$X_ix <- 1:n
  params$H_ix <- (1:n)+n

  H0 <- rep(1e3,n)

  make_parameters_X_SIS(pars = params, b = 0.55, c = 0.15, r = 1/200, Psi = diag(n), wf = rep(1,n), X0 = rep(0,n))
  make_parameters_demography_dynamic(pars = params, b = list(b, b), d = list(d, d), m = list(m, m), H0 = H0)

  y <- rep(0,n*2)
  y[params$H_ix] <- H0

  out_dyn <- deSolve::ode(y = y, times = 0:1000, func = function(t, y, pars, EIR) {
    list(dXdt(t = t, y = y, pars = pars, EIR = EIR))
  }, parms = params, EIR = rep(0, n))

  expect_equal(out_dyn[nrow(out_dyn), params$X_ix+1], out_trace[nrow(out_trace), params$X_ix+1])
  expect_equal(as.numeric(out_dyn[nrow(out_dyn), params$X_ix+1]), rep(0, n))
  expect_equal(as.numeric(out_dyn[nrow(out_dyn), params$H_ix+1]), Ht(1e3))

})


test_that("SIP model; demography with trace is the same as EIR=0 dynamic demography", {

  n <- 4
  b <- c(2e-4, 0.0005, 0.0005/2, 1e-6)
  d <- c(1e-4, 0.0001, 0.0001/2, 1e-6)
  m <- c(1/16, 1/30, 1/14)

  # solve with trace based demography
  calD <- make_calD(d,m,b)
  calD_eig <- eigen(calD)

  params <- new.env()
  params$nStrata <- n
  params$X_ix <- 1:n
  params$P_ix <- (1:n)+n

  H0 <- rep(1e3,n)

  Y <- calD_eig$vectors
  Q <- diag(calD_eig$values)
  Yinv <- solve(Y)

  Ht <- function(t) {
    as.vector(Y %*% (expm::expm(Q * t) %*% (Yinv %*% H0)))
  }

  make_parameters_X_SIP(pars = params, b = 0.55, c = 0.15, r = 1/200, rho = 0.05, eta = 1/32, Psi = diag(n), wf = rep(1,n), X0 = rep(0,n), P0 = rep(0,n))
  make_parameters_demography_trace(pars = params, H = Ht)

  y <- rep(0,n*2)

  out_trace <- deSolve::ode(y = y, times = 0:1000, func = function(t, y, pars, EIR) {
    list(dXdt(t = t, y = y, pars = pars, EIR = EIR))
  }, parms = params, EIR = rep(0, n))

  # solve with dynamic demography
  params <- new.env()
  params$nStrata <- n
  params$X_ix <- 1:n
  params$P_ix <- (1:n)+n
  params$H_ix <- (1:n)+(n*2)

  H0 <- rep(1e3,n)

  make_parameters_X_SIP(pars = params, b = 0.55, c = 0.15, r = 1/200, rho = 0.05, eta = 1/32, Psi = diag(n), wf = rep(1,n), X0 = rep(0,n), P0 = rep(0,n))
  make_parameters_demography_dynamic(pars = params, b = list(b, b, b), d = list(d, d, d), m = list(m, m, m), H0 = H0)

  y <- rep(0,n*3)
  y[params$H_ix] <- H0

  out_dyn <- deSolve::ode(y = y, times = 0:1000, func = function(t, y, pars, EIR) {
    list(dXdt(t = t, y = y, pars = pars, EIR = EIR))
  }, parms = params, EIR = rep(0, n))

  expect_equal(out_dyn[nrow(out_dyn), params$X_ix+1], out_trace[nrow(out_trace), params$X_ix+1])
  expect_equal(as.numeric(out_dyn[nrow(out_dyn), params$X_ix+1]), rep(0, n))
  expect_equal(as.numeric(out_dyn[nrow(out_dyn), params$H_ix+1]), Ht(1e3))

})

# to-do hMoI model.
