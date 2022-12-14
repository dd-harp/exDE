% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/demography-interface.R
\name{F_H}
\alias{F_H}
\title{Size of human population denominators}
\usage{
F_H(t, y, pars)
}
\arguments{
\item{t}{current simulation time}

\item{y}{state vector}

\item{pars}{an \link{environment}}
}
\value{
a \link{numeric} vector of length \code{nStrata}
}
\description{
This method dispatches on the type of \code{pars$Hpar}.
}
