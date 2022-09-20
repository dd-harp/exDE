# spatial metrics of transmission

#' @title Parasite dispersal by mosquitoes
#' @description Compute the `p` by `p` matrix \eqn{\mathcal{V}} whose columns describe
#' how infective bites arising from all the mosquitoes biting a single human on a
#' single day are dispersed to other patches, accounting for movement and mortality.
#' \deqn{\mathcal{V} = fq\Omega^{-1} \cdot e^{-\Omega\tau} \cdot \mbox{diag}\left(\frac{fqM}{W}\right)}
#' @param f the feeding rate
#' @param q fraction of bloodmeals taken on humans
#' @param Omega the mosquito demography matrix
#' @param tau duration of the extrinsic incubation period
#' @param M size of mosquito population in each patch
#' @param W ambient human population at each patch
#' @importFrom MASS ginv
#' @importFrom expm expm
#' @export
metric_calV <- function(f, q, Omega, tau, M, W) {
  stopifnot(inherits(Omega, 'matrix'))
  if (length(f) == 1L) {
    f <- rep(f, nrow(Omega))
  }
  if (length(q) == 1L) {
    q <- rep(q, nrow(Omega))
  }
  OmegaEIP <- expm(x = -Omega*tau)
  Omega_inv <- ginv(X = Omega)
  fq <- diag(f*q, nrow = nrow(Omega))
  fqMW <- diag(f*q*as.vector(M)/as.vector(W), nrow = nrow(Omega))
  return(fq %*% Omega_inv %*% OmegaEIP %*% fqMW)
}

#' @title Parasite dispersal by humans
#' @description Compute the `p` by `p` matrix \eqn{\mathcal{D}} whose columns describe
#' how potentially infectious person time from persons in that patch are dispersed
#' across other patches.
#' \deqn{\mathcal{D} = \mbox{diag}\left(W \right) \cdot \beta^T \cdot \mbox{diag}\left(bDH\right) \cdot \beta}
#' @param W ambient human population at each patch
#' @param beta the biting distribution matrix
#' @param b transmission efficiency from mosquitoes to humans
#' @param D human transmitting capacity
#' @param H human population size of each strata
#' @export
metric_calD <- function(W, beta, b, D, H) {
  stopifnot(inherits(beta, 'matrix'))
  D <- as.vector(D)
  H <- as.vector(H)
  stopifnot(length(D) == length(H))
  if (inherits(W, 'matrix')) {
    W <- diag(as.vector(W), nrow = length(as.vector(W)))
  } else {
    W <- diag(W, nrow = length(W))
  }
  bDH <- diag(b*D*H, nrow = length(H))
  return(W %*% t(beta) %*% bDH %*% beta)
}

#' @title Parasite Dispersal through one Parasite Generation (Humans)
#' @description Computes a `n` by `n` matrix describing parasite dispersal from infecteds (columns)
#' to infectees (rows).
#' \deqn{\mathcal{R} =   b \beta \cdot {\cal V}  \cdot \mbox{diag}\left(W \right) \cdot \beta^T  \cdot \mbox{diag}\left(DH\right)}
#' @param b transmission efficiency from mosquitoes to humans
#' @param beta the biting distribution matrix
#' @param calV parasite dispersal by mosquitoes matrix (see [xDE::metric_calV])
#' @param W ambient human population at each patch
#' @param D human transmitting capacity
#' @param H human population size of each strata
#' @export
metric_calR <- function(b, beta, calV, W, D, H) {
  stopifnot(inherits(beta, 'matrix'))
  D <- as.vector(D)
  H <- as.vector(H)
  stopifnot(length(D) == length(H))
  DH <- diag(D*H, nrow = length(H))
  if (inherits(W, 'matrix')) {
    W <- diag(as.vector(W), nrow = length(as.vector(W)))
  } else {
    W <- diag(W, nrow = length(W))
  }
  return((b*beta) %*% calV %*% W %*% t(beta) %*% DH)
}

#' @title Parasite Dispersal through one Parasite Generation (Mosquitoes)
#' @description Computes a `p` by `p` matrix describing parasite dispersal from infecteds (columns)
#' to infectees (rows).
#' \deqn{\mathcal{Z} = e^{-\Omega\tau} \cdot \mbox{diag}\left( \frac{fq M}{W} \right) \cdot {\cal D} \cdot  fq\Omega^{-1}}
#' @param Omega the mosquito demography matrix
#' @param tau duration of the extrinsic incubation period
#' @param f the feeding rate
#' @param q fraction of bloodmeals taken on humans
#' @param M size of mosquito population in each patch
#' @param W ambient human population at each patch
#' @param calD parasite dispersal by humans matrix (see [xDE::metric_calD])
#' @importFrom MASS ginv
#' @importFrom expm expm
#' @export
metric_calZ <- function(Omega, tau, f, q, M, W, calD) {
  stopifnot(inherits(Omega, 'matrix'))
  stopifnot(inherits(calD, 'matrix'))
  if (length(f) == 1L) {
    f <- rep(f, nrow(Omega))
  }
  if (length(q) == 1L) {
    q <- rep(q, nrow(Omega))
  }
  OmegaEIP <- expm(x = -Omega*tau)
  Omega_inv <- ginv(X = Omega)
  fq <- diag(f*q, nrow = length(f))
  # calc seperate; there might be patches w/out people, set those to 0
  fqMW_diag <- f*q*as.vector(M)/as.vector(W)
  fqMW_diag[!is.finite(fqMW_diag)] <- 0
  fqMW <- diag(fqMW_diag, nrow = nrow(Omega))
  return(OmegaEIP %*% fqMW %*% calD %*% fq %*% Omega_inv)
}
