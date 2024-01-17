
#' @title Make a time spent matrix, called TimeSpent
#' @param pars an [list]
#' @param i the host species index
#' @param TimeSpent a matrix; or a string
#' @param opts a list of options to configure TimeSpent
#' @return a [list]
#' @export
make_TimeSpent = function(pars, i, TimeSpent, opts = list()){
  if(is.matrix(TimeSpent)) class(opts) <- "as_matrix"
  if(is.character(TimeSpent)) class(opts) <- TimeSpent
  UseMethod("make_TimeSpent", opts)
}

#' @title Make a mosquito dispersal matrix, called TimeSpent with a here / away
#' @description Implements [make_TimeSpent] for as_matrix
#' @inheritParams make_TimeSpent
#' @return a [list]
#' @export
make_TimeSpent.athome = function(pars, i, TimeSpent = "athome", opts = list()){

  residence = pars$BFpar$residence[[i]]
  TiSp = make_TimeSpent_athome(pars$nPatches, residence, opts)
  pars$BFpar$TimeSpent[[i]] = TiSp
  for(s in 1:pars$nVectors) pars = make_TaR(t, pars, i, s)

  return(pars)
}

#' @title Make a mosquito dispersal matrix, called TimeSpent
#' @param nPatches is the number of patches
#' @param residence is the home patch for each stratum
#' @param atHome is the fraction of time spent at home
#' @param opts is a set of options that overwrites the defaults
#' @param travel is the fraction of time spent traveling
#' @return a [matrix]
#' @export
make_TimeSpent_athome = function(nPatches, residence, opts=list(), atHome=1, travel=0) {with(opts,{
  nStrata = length(residence)
  away = ifelse(nPatches == 1, 0, (1-atHome-travel)/(nPatches-1))
  atHome = ifelse(nPatches == 1, 1-travel, atHome)
  TiSp <- matrix(away, nPatches, length(residence))
  TiSp[cbind(residence, c(1:nStrata))] <- atHome
  return(TiSp)
})}

#' @title Pass a pre-configured TimeSpent
#' @description Implements [make_TimeSpent] for as_matrix
#' @inheritParams make_TimeSpent
#' @return a [list]
#' @export
make_TimeSpent.as_matrix = function(pars, i, TimeSpent, opts=list()){
  pars$BFpar$TimeSpent[[i]] = TimeSpent
  for(s in 1:pars$nVectors) pars = make_TaR(t, pars, i, s)
  return(pars)
}

#' @title Develop a mosquito dispersal matrix from a kernel and xy-coordinates
#' @description Implements [make_TimeSpent] for kernels
#' @inheritParams make_TimeSpent
#' @return a [list]
#' @export
make_TimeSpent.xy = function(pars, i, TimeSpent = "xy", opts=list()) {
  residence = pars$BFpar[[i]]$residence
  TiSp = with(opts, make_TimeSpent_xy(xy, residence, kern, stay, travel))
  pars$BFpar$TimeSpent[[i]] = TiSp
  for(s in 1:pars$nVectors) pars = make_TaR(t, pars, i, s)
  return(pars)
}

#' @title Make a mosquito dispersal matrix, called TimeSpent
#' @param xy is the xy-locations of the patches
#' @param residence is the home patch for each stratum
#' @param kern is a function that gives weight by distance
#' @param stay is the fraction of time spent at home
#' @param travel is the fraction of time spent traveling
#' @return a [matrix]
#' @export
make_TimeSpent_xy = function(xy, residence, kern, stay, travel) {
  nPatches = dim(xy)[1]
  nStrata = length(residence)
  stopifnot(length(stay)==nStrata)
  stopifnot(length(travel)==nStrata)
  TimeSpent = matrix(0, nPatches, nStrata)
  for(i in 1:nStrata){
    j = residence[i]
    dd = sqrt((xy[j,1] - xy[,1])^2 + (xy[j,2] - xy[,2])^2)
    wts = kern(dd)
    wts[j] = 0
    wts = (1-stay[i]-travel[i])*wts/sum(wts[-j])
    wts[j] = stay[i]
    TimeSpent[,i] = wts
  }
  return(TimeSpent)
}
