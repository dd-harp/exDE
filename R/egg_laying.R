#' @title Setup egg laying for most models
#' @description Sets up the egg-deposition matrix calU for the s^th species
#' @param pars the model object
#' @param s the vector species index
#' @param searchQ the membership matrix
#' @param Lopts a [list] of options to override defaults
#' @return a [list] vector
#' @export
setup_EggLaying = function(pars, s, searchQ=1, Lopts=list()){with(Lopts,{

  gl = list()
  gl$searchQ = checkIt(searchQ, pars$nHabitats, F)
  gl$calU = make_calU(pars$calN, searchQ)
  pars$egg_laying[[s]] = gl
  return(pars)
})}

#' @title Compute eggs laid
#' @description Computes eggs laid for the s^th species
#' @param t the time
#' @param y the state variables
#' @param pars the model object
#' @param s the vector species index
#' @return a [numeric] vector
#' @export
LayEggs = function(t, y, pars, s){
  with(pars$egg_laying[[s]],{
  return(calU %*% F_eggs(t,y,pars,s))
})}
