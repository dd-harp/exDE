# specialized methods for the Le Menach model of ITN based vector control
# https://malariajournal.biomedcentral.com/articles/10.1186/1475-2875-6-10

#' @title Modify baseline values due to vector control
#' @description Implements [VectorControl] for the Le Menach ITN model of vector control
#' @inheritParams VectorControl
#' @return a named [list]
#' @export
VectorControl.lemenach <- function(t, y, pars, MosyBehavior) {

  # evaluate at one or multiple times?
  n <- length(MosyBehavior$f)
  MosyBehaviorControl <- MosyBehavior
  for (i in seq_len(n)) {

    t <- attr(MosyBehavior$f, 'time')[i]

    phi <- pars$VCpar$phi(t)
    tau0 <- F_tau(t, y, pars) # asks the mosy model if it can calculate tau1 and tau2
    if (is.null(tau0)) {
      tau0 <- (1/MosyBehavior$f[i]) * pars$VCpar$tau0_frac
    }

    p0 <- pexp(q = MosyBehavior$g[i] * tau0, lower.tail = FALSE)
    Q0 <- MosyBehavior$q[i]
    W <- (1-Q0) + Q0*(1-phi) + Q0*phi*pars$VCpar$s
    Z <- Q0*phi*pars$VCpar$r

    tau_phi <- tau0
    tau_phi[1] <- tau0[1]/(1-Z)

    f_phi <- 1 / sum(tau_phi) # feeding rate under control

    p_phi <- p0
    p_phi[1] <- (p0[1] * W) / (1 - Z*p0[1])

    g_phi <- -f_phi*log(prod(p_phi)) # mortality under control
    q_phi <- (Q0*(1-phi) + Q0*phi*pars$VCpar$s) / W # human feeding fraction under control

    MosyBehaviorControl$f[i] <- f_phi
    MosyBehaviorControl$q[i] <- q_phi
    MosyBehaviorControl$g[i] <- g_phi
  }

  return(MosyBehaviorControl)
}

#' @title Make parameters for Le Menach ITN model of vector control
#' @param pars an [environment]
#' @param tau0_frac a [numeric] vector giving the proportion of time spent
#' in host seeking/bloodfeeding and resting/oviposition
#' @param r probability of mosquito being repelled upon contact with ITN
#' @param s probability of mosquito successfully feeding upon contact with ITN
#' @return none
#' @export
make_parameters_vc_lemenach <- function(pars, tau0_frac = c(0.68/3, 2.32/3), r = 0.56, s = 0.03) {
  stopifnot(is.environment(pars))
  stopifnot(sum(tau0_frac) == 1)
  VCpar <- list()
  class(VCpar) <- 'lemenach'

  VCpar$tau0_frac <- tau0_frac
  VCpar$r <- r
  VCpar$s <- s

  pars$VCpar <- VCpar
}
