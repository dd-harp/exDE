library(deSolve)

numeric_tol <- 1e-5

test_that("human SIS model remains at equilibrium", {
  nStrata <- 3
  H <- c(100, 500, 250)
  membershipH = 1:nStrata
  searchWtsH = rep(1, nStrata)
  X <- c(20, 120, 80)
  b <- 0.55
  c <- 0.15
  r <- 1/200
  Psi <- matrix(data = 1,nrow = 1, ncol = nStrata)

  EIR <- diag(1/b, nStrata) %*% ((r*X)/(H-X))

  params <- make_parameters_xde()
  params$nStrata <- nStrata
  params$nPatches <- 1

  params = make_parameters_demography_static(pars = params, H=H, membershipH=membershipH,
                                           searchWtsH=searchWtsH, TimeSpent=Psi,
                                           birthF = "constant", birthrate = 0,
                                           Hmatrix = diag(1, length(H)),
                                           birthsXstrata = rep(1,length(H)))
  params = make_parameters_X_SIS(pars = params, b = b, c = c, r = r)
  params = make_inits_X_SIS(pars = params, X)

  params = make_indices(params)


  # set initial conditions
  y0 <- get_inits(params)

  y0 <- rep(0, 3)
  y0[params$X_ix] <- X

  out <- deSolve::ode(y = y0, times = c(0, 365), func = function(t, y, pars, EIR) {
    list(dXdt(t, y, pars, EIR))
  }, parms = params, method = 'lsoda', EIR = as.vector(EIR))

  expect_equal(as.vector(out[2L, params$X_ix+1]), X, tolerance = numeric_tol)
})
