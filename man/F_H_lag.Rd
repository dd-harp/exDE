% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/demography-interface.R
\name{F_H_lag}
\alias{F_H_lag}
\title{Size of lagged human population denominators}
\usage{
F_H_lag(t, y, pars, lag)
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
This method dispatches on the type of \code{pars$Hpar}.
}
