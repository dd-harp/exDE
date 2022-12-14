% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adult-GeRM.R
\name{F_Z_lag.GeRM}
\alias{F_Z_lag.GeRM}
\title{Density of lagged infectious mosquitoes}
\usage{
\method{F_Z_lag}{GeRM}(t, y, pars, lag)
}
\arguments{
\item{t}{current simulation time}

\item{y}{state vector}

\item{pars}{an \link{environment}}

\item{lag}{duration of lag \code{t-lag}}
}
\value{
a \link{numeric} vector of length \code{nPatches}
}
\description{
Implements \link{F_Z_lag} for the generalized GeRM model.
}
