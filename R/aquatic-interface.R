# generic methods for aquatic component


#' @title Set aquatic bionomic parameter rates relative to baseline
#' @description This method dispatches on the type of `pars$Lpar`. It should
#' compute the values of parameters as a function of exogenous variables
#' or reset the values of the bionomic parameters to baseline values.
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @return a [list]
#' @export
LBionomics <- function(t, y, pars) {
  UseMethod("LBionomics", pars$Lpar)
}

#' @title Number of newly emerging adults from each larval habitat
#' @description This method dispatches on the type of `pars$Lpar`.
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @return a [numeric] vector of length `nHabitats`
#' @export
F_alpha <- function(t, y, pars) {
  UseMethod("F_alpha", pars$Lpar)
}

#' @title Derivatives for aquatic stage mosquitoes
#' @description This method dispatches on the type of `pars$Lpar`.
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @param eta vector giving number of eggs being laid in each larval habitat
#' @return a [numeric] vector of length `pars$L_ix`
#' @export
dLdt <- function(t, y, pars, eta) {
  UseMethod("dLdt", pars$Lpar)
}

#' @title A function to set up Lpar
#' @description This method dispatches on `Lname`.
#' @param pars a [list]
#' @param Lname a [character] string
#' @param membership a [vector] of length `nHabitats`
#' @param searchQ a [vector] of length `nHabitats`
#' @param Lopts a [list]
#' @return none
#' @export
setup_L = function(pars, Lname,
                      membership=1, searchQ=1,
                      Lopts=list()){
  class(Lname) <- Lname
  UseMethod("setup_L", Lname)
}

#' @title Add indices for aquatic stage mosquitoes to parameter list
#' @description This method dispatches on the type of `pars$Lpar`. Adds field `L_ix`
#' to parameter list.
#' @param pars a [list]
#' @return none
#' @export
make_indices_L <- function(pars) {
  UseMethod("make_indices_L", pars$Lpar)
}

#' @title Parse the output of deSolve and return the variables by name in a list
#' @description This method dispatches on the type of `pars$Lpar`. Attaches the
#' state variables for the aquatic ecology model to varslist and returns it
#' @param varslist a [list] the object to be returned
#' @param deout a [matrix] of outputs from deSolve
#' @param pars a [list] that defines a model
#' @export
parse_deout_L <- function(varslist, deout, pars) {
  UseMethod("parse_deout_L", pars$Lpar)
}

#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Lpar`.
#' @param pars a [list]
#' @return none
#' @export
get_inits_L <- function(pars) {
  UseMethod("get_inits_L", pars$Lpar)
}

#' @title Make the habitat membership matrix, calN
#' @param nPatches is the number of patches
#' @param membership is a vector describing the patch where each habitat is found
#' @return a [matrix] of dimensions `nPatches` by `nHabitats`
#' @export
make_calN = function(nPatches, membership){
  nHabitats = length(membership)
  calN = matrix(0, nPatches, nHabitats)
  calN[cbind(membership, 1:nHabitats)]=1
  return(calN)
}

#' @title Make the egg distribution matrix, calU
#' @param calN the habitat membership matrix
#' @param searchQ the habitat search weights
#' @return a [matrix] of dimensions `nHabitats` by `nPatches`
#' @export
make_calU = function(calN, searchQ=1){
  calU = searchQ*t(calN)
  colNorms = colSums(calU)
  ix = which(colNorms == 0)
  if(length(ix)>0) colNorms[ix]=1
  calU = calU %*% diag(1/colNorms)
  return(calU)
}
