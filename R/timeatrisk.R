
#' @title Make a mosquito dispersal matrix, called TaR
#' @param nPatches is the number of patches
#' @param residence is the residence vector
#' @param TaR a matrix; or a string
#' @param opts a list of options to configure TaR
#' @return a [matrix]
#' @export
make_TaR = function(nPatches, residence, TaR, opts = list()){
  if(is.matrix(TaR)) class(opts) <- "as_matrix"
  if(is.character(TaR)) class(opts)<- TaR
  UseMethod("make_TaR", opts)
}

#' @title Make a mosquito dispersal matrix, called TaR with a here / away
#' @description Implements [make_TaR] for as_matrix
#' @inheritParams make_TaR
#' @return a [matrix]
#' @export
make_TaR.athome = function(nPatches, residence, TaR = "athome", opts = list()){
  with(opts,make_TaR_athome(nPatches, residence, opts))
}

#' @title Make a mosquito dispersal matrix, called TaR
#' @param nPatches is the number of patches
#' @param residence is the home patch for each stratum
#' @param atHome is the fraction of time spent at home
#' @param opts is a set of options that overwrites the defaults
#' @param travel is the fraction of time spent traveling
#' @return a [matrix]
#' @export
make_TaR_athome = function(nPatches, residence, opts, atHome=0.95, travel=0.01) {with(opts,{
  nStrata = length(residence)
  away = ifelse(nPatches == 1, 0, (1-atHome-travel)/(nPatches-1))
  atHome = ifelse(nPatches == 1, 1-travel, atHome)
  TiSp <- matrix(away, nPatches, length(residence))
  TiSp[cbind(residence, c(1:nStrata))] <- atHome
  return(TiSp)
})}

#' @title Pass a pre-configured TaR
#' @description Implements [make_TaR] for as_matrix
#' @inheritParams make_TaR
#' @return a [matrix]
#' @export
make_TaR.as_matrix = function(nPatches, residence, TaR, opts=list()){
  return(TaR)
}

#' @title Develop a mosquito dispersal matrix from a kernel and xy-coordinates
#' @description Implements [make_TaR] for kernels
#' @inheritParams make_TaR
#' @return a [matrix]
#' @export
make_TaR.xy = function(nPatches, residence, TaR = "xy", opts=list()) {
  with(opts,make_TaR_xy(xy, residence, kern, stay, travel))
}

#' @title Make a mosquito dispersal matrix, called TaR
#' @param xy is the xy-locations of the patches
#' @param residence is the home patch for each stratum
#' @param kern is a function that gives weight by distance
#' @param stay is the fraction of time spent at home
#' @param travel is the fraction of time spent traveling
#' @return a [matrix]
#' @export
make_TaR_xy = function(xy, residence, kern, stay, travel) {
  nPatches = dim(xy)[1]
  nStrata = length(residence)
  stopifnot(length(stay)==nStrata)
  stopifnot(length(travel)==nStrata)
  TaR = matrix(0, nPatches, nStrata)
  for(i in 1:nStrata){
    j = residence[i]
    dd = sqrt((xy[j,1] - xy[,1])^2 + (xy[j,2] - xy[,2])^2)
    wts = kern(dd)
    wts[j] = 0
    wts = (1-stay[i]-travel[i])*wts/sum(wts[-j])
    wts[j] = stay[i]
    TaR[,i] = wts
  }
  return(TaR)
}
