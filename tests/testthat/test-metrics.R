numeric_tol <- sqrt(.Machine$double.eps)

test_that("metrics calculations work in 1 patch/strata", {
  f <- 0.3
  q <- 0.9
  g <- 1/10
  tau <- 11
  p <- 1
  Omega <- diag(g, p)

  M <- 500
  W <- 800

  vc <- (f*q/g) * exp(-g*tau) * (f*q*M)/W
  calV <- metric_calV(f = f, q = q, Omega = Omega, tau = tau, M = M, W = W)

  expect_equal(as.vector(calV), vc, tolerance = numeric_tol)
  Psi <- matrix(1,1,1)
  beta <- t(Psi) %*% diag_inverse(W)
  b <- 0.55
  c <- 0.15
  r <- 1/200
  D <- c/r

  calD <- metric_calD(W = W, beta = beta, b = b, D = D, H = W)
  expect_equal(b*c/r, as.vector(calD), tolerance = numeric_tol)

  calR <- metric_calR(b = b, beta = beta, calV = calV, W = W, D = D, H = W)
  r0 <- (f*q/g) * exp(-g*tau) * (f*q*M)/W * (b*c/r)
  expect_equal(r0, as.vector(calR), tolerance = numeric_tol)

  calZ <- metric_calZ(Omega = Omega, tau = tau, f = f, q = q, M = M, W = W, calD = calD)
  expect_equal(r0, as.vector(calZ), tolerance = numeric_tol)
})
