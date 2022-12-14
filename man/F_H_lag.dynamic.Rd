% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/demography-dynamic.R
\name{F_H_lag.dynamic}
\alias{F_H_lag.dynamic}
\title{Size of lagged human population denominators}
\usage{
\method{F_H_lag}{dynamic}(t, y, pars, lag)
}
\arguments{
\item{t}{current simulation time}

\item{y}{state vector}

\item{pars}{an \link{environment}}

\item{lag}{duration of lag \code{t-lag}}
}
\value{
a \link{numeric} vector of length \code{nStrata}
}
\description{
Implements \link{F_H_lag} for dynamic models.
}
