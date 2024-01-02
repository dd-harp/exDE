# generic methods for aquatic component


#' @title Set aquatic bionomic parameter rates relative to baseline
#' @description This method dispatches on the type of `pars$Lpar`. It should
#' compute the values of parameters as a function of exogenous variables
#' or reset the values of the bionomic parameters to baseline values.
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @param s the species index
#' @return a [list]
#' @export
LBionomics <- function(t, y, pars, s) {
  UseMethod("LBionomics", pars$Lpar[[s]])
}

#' @title Number of newly emerging adults from each larval habitat
#' @description This method dispatches on the type of `pars$Lpar`.
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @param s the species index
#' @return a [numeric] vector of length `nHabitats`
#' @export
F_alpha <- function(t, y, pars, s) {
  UseMethod("F_alpha", pars$Lpar[[s]])
}

#' @title Derivatives for aquatic stage mosquitoes
#' @description This method dispatches on the type of `pars$Lpar`.
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @param s the species index
#' @param eta vector giving number of eggs being laid in each larval habitat
#' @return a [numeric] vector of length `pars$L_ix`
#' @export
dLdt <- function(t, y, pars, eta, s) {
  UseMethod("dLdt", pars$Lpar[[s]])
}


#' @title A function to set up adult mosquito models
#' @description This method dispatches on `Lname`.
#' @param Lname the class name of the aquatic model
#' @param pars a [list]
#' @param s the species index
#' @param Lopts a [list]
#' @return [list]
#' @export
setup_Lpar = function(Lname, pars, s, Lopts=list()){
  class(Lopts) <- Lname
  UseMethod("setup_Lpar", Lopts)
}

#' @title A function to set up adult mosquito models
#' @description This method dispatches on `Lname`.
#' @param pars a [list]
#' @param s the species index
#' @param Lopts a [list]
#' @return [list]
#' @export
setup_Linits = function(pars, s, Lopts=list()){
  UseMethod("setup_Linits", pars$Lpar[[s]])
}

#' @title Add indices for aquatic stage mosquitoes to parameter list
#' @description This method dispatches on the type of `pars$Lpar`. Adds field `L_ix`
#' to parameter list.
#' @param pars a [list]
#' @param s the species index
#' @return none
#' @export
make_indices_L <- function(pars, s) {
  UseMethod("make_indices_L", pars$Lpar[[s]])
}

#' @title Parse the output of deSolve and return the variables by name in a list
#' @description This method dispatches on the type of `pars$Lpar`. Attaches the
#' state variables for the aquatic ecology model to a list and returns it
#' @param deout a [matrix] of outputs from deSolve
#' @param pars a [list] that defines the model
#' @param s the species index
#' @export
parse_deout_L <- function(deout, pars, s) {
  UseMethod("parse_deout_L", pars$Lpar[[s]])
}

#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Lpar`.
#' @param pars a [list]
#' @param s the species index
#' @return none
#' @export
get_inits_L <- function(pars, s) {
  UseMethod("get_inits_L", pars$Lpar[[s]])
}

#' @title Set the initial values from a vector of model states
#' @description This method dispatches on the type of `pars$Lpar`.
#' @param pars a [list]
#' @param y0 a vector of variable values from a simulation
#' @param s the species index
#' @return none
#' @export
update_inits_L <- function(pars, y0, s) {
  UseMethod("update_inits_L", pars$Lpar[[s]])
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
