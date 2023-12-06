
#' @title Make base parameters
#' @param solve_as, either "ode" or "dde"
#' @return a [list]
#' @export
make_parameters_xde = function(solve_as='ode'){
  pars = list()

  xde <- solve_as
  class(xde) <- xde
  pars$xde = xde

  pars <- setup_abiotic_null(pars)
  pars <- setup_shock_null(pars)
  pars <- setup_control_null(pars)
  pars <- setup_vc_null(pars)
  pars <- setup_behavior_null(pars)
  pars <- setup_visitors_null(pars)
  pars <- setup_resources_null(pars)
  pars <- setup_eip_static(pars)
  pars <- setup_travel_null(pars)
  pars <- setup_exposure_pois(pars)

  return(pars)
}

#' @title Set indices for generalized spatial model
#' @param pars a [list]
#' @return none
#' @export
make_indices <- function(pars) {
  pars$max_ix <- 0
  if ('Lpar' %in% names(pars)) {
    pars = make_indices_L(pars)
  }
  if ('MYZpar' %in% names(pars)) {
    pars = make_indices_MYZ(pars)
  }
  if ('Xpar' %in% names(pars)) {
    pars = make_indices_X(pars)
  }
  if ('Hpar' %in% names(pars)) {
    pars = make_indices_H(pars)
  }
  return(pars)
}

#' @title Get the initial values as a vector
#' @param pars a [list]
#' @return y a [numeric] vector assigned the class "dynamic"
#' @export
get_inits <- function(pars){
  if ('Lpar' %in% names(pars)) {
    Li = get_inits_L(pars)
  } else {Li = numeric(0)}
  if ('MYZpar' %in% names(pars)) {
    MYZi = get_inits_MYZ(pars)
  } else {MYZi = numeric(0)}
  if ('Xpar' %in% names(pars)) {
    Xi = get_inits_X(pars)
  } else {Xi = numeric(0)}
  if ('Hpar' %in% names(pars)) {
    Hi = get_inits_H(pars)
  } else {Hi = numeric(0)}
  y = c(L=Li, MYZ=MYZi, X=Xi, H=Hi)
  return(y)
}

#' @title Parse the output of an object returned by deSolve
#' @param pars a [list]
#' @param deout a [matrix] of orbits returned by deSolve
#' @return varslist a [list]
#' @export
parse_deout <- function(deout, pars){
  varslist = list()
  varslist$time = deout[,1]
  if ('Lpar' %in% names(pars)) {
    varslist$L = parse_deout_L(deout, pars)
  }
  if ('MYZpar' %in% names(pars)) {
    varslist$MYZ = parse_deout_MYZ(deout, pars)
  }
  if ('Xpar' %in% names(pars)) {
    varslist$XH = parse_deout_X(deout, pars)
  }
  varslist$terms = compute_terms(deout, pars)
  varslist$deout = deout
  return(varslist)
}

#' @title Parse the output of an object returned by deSolve
#' @param vec a [vector] with the variables, as returned by rootsolve
#' @param pars a [list]
#' @return varslist a [list]
#' @export
parse_deout_vec <- function(vec, pars){
  varslist = list()
  deout = rbind(c(0,vec), c(0,vec))
  varslist$vec = deout
  if ('Lpar' %in% names(pars)) {
    varslist$L = parse_deout_L(deout, pars)
  }
  if ('MYZpar' %in% names(pars)) {
    varslist$MYZ = parse_deout_MYZ(deout, pars)
  }
  if ('Xpar' %in% names(pars)) {
    varslist$X = parse_deout_X(deout, pars)
  }
  for(i in 1:length(varslist)){
    varslist[[i]] = tail(varslist[[i]],1)
  }
  varslist$terms = compute_terms_steady(vec, pars)
  return(varslist)
}


#' @title Invert a diagonal matrix
#' @description Invert a diagonal matrix which is passed as a vector. If any
#' elements are zero, set them to one.
#' @param x a [numeric] vector
#' @return a diagonal [matrix]
#' @export
diag_inverse <- function(x) {
  x <- as.vector(x)
  ix <- which(x == 0)
  if (length(ix) > 0) {
    x[ix] <- 1
  }
  return(diag(x = 1/x, nrow = length(x), ncol = length(x), names = FALSE))
}

#' @title Check if two numeric values are approximately equal
#' @param a a [numeric] object
#' @param b a [numeric] object
#' @param tol the numeric tolerance
#' @return a [logical] value
#' @export
approx_equal <- function(a, b, tol = sqrt(.Machine$double.eps)) {
  abs(a - b) < tol
}

#' @title Check the length of an input value
#' @param x a [numeric] object
#' @param lng a [numeric] object
#' @param type a [character] string specifying required typeof
#' @param fixit a [logical] value, if TRUE force length to lng
#' @return a [numeric] object
#' @export
checkIt = function(x, lng, type = "numeric", fixit=TRUE){
  stopifnot(is.numeric(x))
  if(type == "integer") x = as.integer(x)
  if(length(x)==1 & fixit) x=rep(x, lng)
  stopifnot(length(x)==lng)
  x
}

#' @title Set the initial values to the last values of the last simulation
#' @param pars a [list]
#' @return y a [numeric] vector
#' @export
last_to_inits <- function(pars){
  y0 <- tail(pars$orbits$deout, 1)[-1]
  pars <- update_inits(y0, pars)
  return(pars)
}

#' @title Set the initial values to the last values of the last simulation
#' @param y0 a [vector] of initial values
#' @param pars a [list]
#' @return y a [numeric] vector
#' @export
update_inits <- function(y0, pars){
  if ('Lpar' %in% names(pars)) {
    pars = update_inits_L(pars, y0)
  } else {Li = numeric(0)}
  if ('MYZpar' %in% names(pars)) {
    pars = update_inits_MYZ(pars, y0)
  } else {MYZi = numeric(0)}
  if ('Xpar' %in% names(pars)) {
    pars = update_inits_X(pars, y0)
  } else {Xi = numeric(0)}
  if ('Hpar' %in% names(pars)) {
    pars = update_inits_H(pars, y0)
  } else {Hi = numeric(0)}
  return(pars)
}
