
#' @title Make a mosquito dispersal matrix, called calK
#' @param nPatches is the number of patches
#' @param calK a matrix or a string
#' @param opts a list of options to configure calK
#' @return a [matrix]
#' @export
make_calK = function(nPatches, calK, opts = list()){
  if(is.matrix(calK)) class(opts) <- "as_matrix"
  if(is.character(calK)) class(opts)<- calK
  UseMethod("make_calK", opts)
}

#' @title Dispersal to every other patch, with equal probability
#' @description Implements [make_calK] for the herethere model
#' @inheritParams make_calK
#' @return a [matrix]
#' @export
make_calK.herethere = function(nPatches, calK = "herethere", opts = list()){
  calK <- matrix(1/(nPatches-1), nPatches, nPatches)
  diag(calK) <- 0
  return(calK)
}

#' @title Pass a pre-configured calK
#' @description Implements [make_calK] for as_matrix
#' @inheritParams make_calK
#' @return a [matrix]
#' @export
make_calK.as_matrix = function(nPatches, calK, opts=list()){
  return(calK)
}

#' @title Develop a mosquito dispersal matrix from a kernel and xy-coordinates
#' @description Implements [make_calK] for kernels
#' @inheritParams make_calK
#' @return a [matrix]
#' @export
make_calK.xy = function(nPatches, calK = "xy", opts=list()) {
  with(opts, make_calK_xy(xy, ker))
}

#' @title Develop a mosquito dispersal matrix from a kernel and xy-coordinates
#' @param xy is a vector of the xy-coordinates of patch locations
#' @param ker is a function that weights putative locations by distance
#' @export
make_calK_xy = function(xy, ker) {
  dmat <- as.matrix(stats::dist(xy), upper=T)
  calK <- ker(dmat)
  diag(calK) <- 0
  calK = calK %*% diag(1/rowSums(calK))
  return(calK)
}
