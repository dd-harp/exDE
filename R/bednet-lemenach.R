# specialized methods for the Le Menach model of ITN based vector control
# https://malariajournal.biomedcentral.com/articles/10.1186/1475-2875-6-10

#' @title Modify baseline values due to vector control
#' @description Implements [BedNets] for the Le Menach ITN model of vector control
#' @inheritParams BedNets
#' @return a named [list]
#' @importFrom stats pexp
#' @export
BedNets.lemenach <- function(t, pars) {

  n <- length(pars$MYZpar$f)

  for (i in seq_len(n)) {

    phi <- pars$ITNpar$cover(t)

    tau0 <- (1/pars$MYZpar$f[i]) * pars$ITNpar$tau0_frac

    p0 <- pexp(q = pars$MYZpar$g[i] * tau0, lower.tail = FALSE)
    Q0 <- pars$MYZpar$q[i]
    W <- (1-Q0) + Q0*(1-phi) + Q0*phi*pars$ITNpar$s
    Z <- Q0*phi*pars$ITNpar$r

    tau_phi <- tau0
    tau_phi[1] <- tau0[1]/(1-Z)

    f_phi <- 1 / sum(tau_phi) # feeding rate under control

    p_phi <- p0
    p_phi[1] <- (p0[1] * W) / (1 - Z*p0[1])

    g_phi <- -f_phi*log(prod(p_phi)) # mortality under control
    q_phi <- (Q0*(1-phi) + Q0*phi*pars$ITNpar$s) / W # human feeding fraction under control

    pars$MYZpar$f[i] <- f_phi
    pars$MYZpar$q[i] <- q_phi
    pars$MYZpar$g[i] <- g_phi
  }

  return(pars)
}

#' @title Make parameters for Le Menach ITN model of vector control
#' @description This model of ITN based vector control was originally described in \url{https://malariajournal.biomedcentral.com/articles/10.1186/1475-2875-6-10}.
#' @param pars a [list]
#' @param tau0_frac a [numeric] vector giving the proportion of time spent
#' in host seeking/bloodfeeding and resting/oviposition
#' @param r probability of mosquito being repelled upon contact with ITN
#' @param s probability of mosquito successfully feeding upon contact with ITN
#' @param phi a [function] that takes a single argument `t` and returns the level of ITN coverage at that time
#' @return none
#' @export
setup_itn_lemenach <- function(pars, tau0_frac = c(0.68/3, 2.32/3), r = 0.56, s = 0.03, phi = function(t) {.8} ) {
  stopifnot(sum(tau0_frac) == 1)
  stopifnot(phi(0) >= 0)
  stopifnot(phi(0) <= 1)
  ITNpar <- list()
  class(ITNpar) <- 'lemenach'

  ITNpar$tau0_frac <- tau0_frac
  ITNpar$r <- r
  ITNpar$s <- s
  ITNpar$cover <- phi

  pars$ITNpar <- ITNpar
  return(pars)
}
