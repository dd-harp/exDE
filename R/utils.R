
#' @title Make base parameters, assuming nVectors = nHosts = 1
#' @param solve_as, either "ode" or "dde"
#' @return a [list]
#' @export
make_parameters_xde = function(solve_as='dde'){
  pars = list()

  xde <- solve_as
  class(xde) <- xde
  pars$xde = xde

  pars$MYZpar = list()
  pars$Lpar = list()
  pars$Xpar = list()
  pars$Hpar = list()
  pars$vars = list()

  pars$Lambda = list()
  pars <- setup_EGGpar_static(pars)
  pars <- setup_BFpar_static(pars)

  pars$Linits = list()
  pars$MYZinits = list()
  pars$Xinits = list()

  pars$ix = list()
  pars$ix$X = list()
  pars$ix$MYZ = list()
  pars$ix$L = list()


  pars$outputs = list()
  pars$compute = list()

  pars$HostAvailability = list()

  pars <- setup_abiotic_null(pars)
  pars <- setup_shock_null(pars)
  pars <- setup_control_null(pars)
  pars <- setup_vc_null(pars)
  pars <- setup_behavior_null(pars)
  pars <- setup_habitat_dynamics_static(pars)
  pars <- setup_bionomics_static(pars)
  pars <- setup_visitors_static(pars)
  pars <- setup_resources_null(pars)
  pars <- setup_travel_static(pars)
  pars <- setup_exposure_pois(pars)

  return(pars)
}

#' @title Set indices for generalized spatial model
#' @param pars a [list]
#' @return none
#' @export
make_indices <- function(pars) {
  pars$max_ix <- 0

  s = length(pars$Linits)
  if(s>0)
    for(ix in 1:s)
      pars = make_indices_L(pars, ix)

  s = length(pars$MYZinits)
  if(s>0)
    for(ix in 1:s)
      pars = make_indices_MYZ(pars, ix)

  i = length(pars$Xinits)
  if(i>0)
    for(ix in 1:i)
      pars = make_indices_X(pars, ix)

  return(pars)
}

#' @title Get the initial values as a vector
#' @param pars a [list]
#' @return y a [numeric] vector assigned the class "dynamic"
#' @export
get_inits <- function(pars){

  Li = c()
  s = length(pars$Lpar)
  if(s>0)
    for(ix in 1:s)
      Li = c(Li, get_inits_L(pars, ix))

  MYZi = c()
  s = length(pars$MYZpar)
  if(s>0)
    for(ix in 1:s)
      MYZi = c(MYZi, get_inits_MYZ(pars, ix))

  Xi = c()
  i = length(pars$Xpar)
  if(i>0)
    for(ix in 1:i)
      Xi = c(Xi, get_inits_X(pars, ix))

  y = c(L=Li, MYZ=MYZi, X=Xi)
  return(y)
}

#' @title Parse the output of an object returned by deSolve
#' @param pars a [list]
#' @param deout a [matrix] of orbits returned by deSolve
#' @return varslist a [list]
#' @export
parse_deout <- function(deout, pars){
  varslist = list()

  s = length(pars$Lpar)
  if(s>0)
    for(ix in 1:s)
      varslist$L[[ix]]= parse_deout_L(deout, pars, ix)

  s = length(pars$MYZpar)
  if(s>0)
    for(ix in 1:s)
      varslist$MYZ[[ix]]= parse_deout_MYZ(deout, pars, ix)

  s = length(pars$Xpar)
  if(s>0)
    for(ix in 1:s)
      varslist$XH[[ix]]= parse_deout_X(deout, pars, ix)

  varslist$terms = compute_terms(varslist, deout, pars, 1, 1)
  varslist$deout = deout
  return(varslist)
}

#' @title Parse the output of an object returned by deSolve
#' @param vec a [vector] with the variables, as returned by rootsolve
#' @param pars a [list]
#' @return varslist a [list]
#' @export
parse_deout_vec <- function(vec, pars){
  deout = rbind(c(0,vec), c(0,vec))
  varslist = parse_deout(deout, pars)

  for(i in 1:length(varslist$XH))
    varslist$XH[[i]] = tail(varslist$XH[[i]],1)
  for(i in 1:length(varslist$MYZ))
    varslist$MYZ[[i]] = tail(varslist$MYZ[[i]],1)
  for(i in 1:length(varslist$L))
    varslist$L[[i]] = tail(varslist$L[[i]],1)

  varslist$terms = compute_terms_steady(varslist, vec, pars)
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

#' @title Check the shape and dimensions of an object
#' @param obj a [numeric] object
#' @param d1 an [integer]
#' @param d2 an [integer]
#' @return [matrix]
#' @export
shapeIt = function(obj, d1, d2){
  Obj = as.matrix(obj)
  dd = dim(Obj)
  stopifnot(d1 %in% dd)
  stopifnot(d2 %in% dd)
  if(dd[1]!=d1) obj = t(obj)
  return(obj)
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
  s = length(pars$Lpar)
  if(s>0)
    for(ix in 1:s)
      pars = update_inits_L(pars, y0)

  s = length(pars$MYZpar)
  if(s>0)
    for(ix in 1:s)
      pars = update_inits_MYZ(pars, y0)

  s = length(pars$Xpar)
  if(s>0)
    for(ix in 1:s)
      pars = update_inits_X(pars, y0)

  return(pars)
}
