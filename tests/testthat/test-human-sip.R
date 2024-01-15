library(deSolve)

numeric_tol <- 1e-5

test_that("human SIP model remains at equilibrium", {
  nStrata <- 3
  H <- c(100, 500, 250)
  residence = 1:nStrata
  searchWtsH = rep(1, nStrata)
  I <- c(20, 120, 80)
  b <- 0.55
  c <- 0.15
  r <- 1/200
  eta <- c(1/30, 1/40, 1/35)
  rho <- c(0.05, 0.1, 0.15)
  xi <- rep(0, 3)
  TaR <- matrix(data = 1,nrow = 1, ncol = nStrata)

  P <- diag(1/eta) %*% diag(rho/(1-rho)) %*% (r*I)
  EIR <- diag(1/b, nStrata) %*% diag(1/(1-rho)) %*% ((r*I)/(H-I-P))
  foi <- b*EIR

  params <- make_parameters_xde()
  params$nStrata <- nStrata
  params$nHosts <- 1
  params$nPatches <- 1

  params = make_parameters_demography_null(pars = params, H=H)
  params = setup_BloodFeeding(params, 1, 1, residence=residence, searchWts=searchWtsH)
  params$BFpar$TaR[[1]][[1]]=TaR
  params = make_parameters_X_SIP(pars = params, b = b, c = c, r = r, eta=eta, rho=rho, xi=xi)
  params = make_inits_X_SIP(pars = params, H-I-P, I, P)

  params = make_indices(params)
  params$FoI[[1]] <- foi


  # set initial conditions
  y0 <- get_inits(params)


  out <- deSolve::ode(y = y0, times = c(0, 365), func = function(t, y, pars, s) {
    list(dXdt(t, y, pars, s))
  }, parms = params, method = 'lsoda', s=1)

  expect_equal(as.vector(out[2L, params$ix$X[[1]]$I_ix+1]), I, tolerance = numeric_tol)
  expect_equal(as.vector(out[2L, params$ix$X[[1]]$P_ix+1]), as.vector(P), tolerance = numeric_tol)
})
