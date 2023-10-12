% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adult-GeRM.R
\name{parse_deout_MYZ.GeRM}
\alias{parse_deout_MYZ.GeRM}
\title{Parse the output of deSolve and return variables for the GeRM model}
\usage{
\method{parse_deout_MYZ}{GeRM}(varslist, deout, pars)
}
\arguments{
\item{varslist}{a \link{list} the object to be returned}

\item{deout}{a \link{matrix} of outputs from deSolve}

\item{pars}{a \link{list} that defines a model}
}
\value{
none
}
\description{
Implements \link{parse_deout_MYZ} for the GeRM model.
}