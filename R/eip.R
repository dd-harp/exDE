# generic methods to compute the extrinsic incubation period (EIP)

#' @title Compute the EIP
#' @description This method dispatches on the type of `EIPmod`.
#' @param t current simulation time
#' @param EIPmod a [list]
#' @return [numeric]
#' @export
EIP <- function(t, EIPmod) {
  UseMethod("EIP", EIPmod)
}

#' @title Compute the derivative of the EIP as a function of time
#' @description This method dispatches on the type of `EIPmod`.
#' @param t current simulation time
#' @param EIPmod a [list]
#' @return [numeric]
#' @export
dEIPdt <- function(t, EIPmod) {
  UseMethod("dEIPdt", EIPmod)
}

#' @title Modify parameters due to exogenous forcing by all kinds of control
#' @description Implements [EIP] for the static model (the EIP is constant)
#' @inheritParams EIP
#' @return [numeric]
#' @export
EIP.static <- function(t, EIPmod){EIPmod$eip}

#' @title This function computes the negative derivative of the EIP as a function of time
#' @description Implements [dEIPdt] for the static model (the dEIPdt=0)
#' @inheritParams dEIPdt
#' @return [numeric]
#' @export
dEIPdt.static <- function(t, EIPmod){0}

#' @title Set up the static model for control forcing (do nothing)
#' @param EIPname the class name of the function
#' @param EIPopts is a [list] that overwrites default options
#' @return [list]
#' @export
setup_EIP <- function(EIPname = 'static', EIPopts = list()) {
  class(EIPname) <- EIPname
  UseMethod("setup_EIP", EIPname)
}

#' @title Set up a static model for the EIP
#' @inheritParams setup_EIP
#' @return [list]
#' @export
setup_EIP.static<- function(EIPname, EIPopts=list()){
  setup_eip_static(EIPopts)
}

#' @title Set up a static model for the EIP
#' @param EIPopts a [list]
#' @param eip the extrinsic incubation period (in days)
#' @return [list]
#' @export
setup_eip_static = function(EIPopts=list(), eip=11){with(EIPopts,{
  EIPmod <- list()
  class(EIPmod) <- 'static'
  EIPmod$eip = eip
  return(EIPmod)
})}
