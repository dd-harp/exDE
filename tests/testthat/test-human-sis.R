library(deSolve)

numeric_tol <- 1e-5

test_that("human SIS model remains at equilibrium", {
  nStrata <- 3
  H <- c(100, 500, 250)
  residence = 1:nStrata
  searchWtsH = rep(1, nStrata)
  I <- c(20, 120, 80)
  b <- 0.55
  c <- 0.15
  r <- 1/200
  TaR <- matrix(data = 1,nrow = 1, ncol = nStrata)

  EIR <- diag(1/b, nStrata) %*% ((r*I)/(H-I))
  foi <- b*EIR

  params <- make_parameters_xde()
  params$nStrata <- nStrata
  params$nHosts <- 1
  params$nPatches <- 1

  params = make_parameters_demography_null(pars = params, H=H)
  params = setup_BloodFeeding(params, 1, 1, residence=residence, searchWts=searchWtsH)
  params$BFpar$TaR[[1]][[1]]=TaR
  params = make_parameters_X_SIS(pars = params, b = b, c = c, r = r)
  params = make_inits_X_SIS(pars = params, H-I, I)

  params = make_indices(params)
  params$FoI[[1]] <- foi

  # set initial conditions
  y0 <- get_inits(params)

  out <- deSolve::ode(y = y0, times = c(0, 365), func = function(t, y, pars, s) {
    list(dXdt(t, y, pars, s))
  }, parms = params, method = 'lsoda', s=1)

  expect_equal(as.vector(out[2L, params$ix$X[[1]]$I_ix+1]), I, tolerance = numeric_tol)
})
