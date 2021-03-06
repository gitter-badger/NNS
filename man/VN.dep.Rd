% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Dependence.R
\name{VN.dep}
\alias{VN.dep}
\title{VN Dependence}
\usage{
VN.dep(x, y, order = 2, degree = 0, print.map = TRUE)
}
\arguments{
\item{x}{Variable 1}

\item{y}{Variable 2}

\item{order}{Number of partial moment quadrants to be generated}

\item{degree}{Defaults to 0 for smaller number of observations}

\item{print.map}{Displays partition mapping onto plot.  Defaults to TRUE.}
}
\description{
Returns the dependence between two variables based on higher order partial moment correlations measured by frequency or area.
}
\examples{
set.seed(123)
x<-rnorm(100); y<-rnorm(100)
\dontrun{VN.dep(x,y)}
}
\author{
Fred Viole, OVVO Financial Systems
}
\references{
Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
\url{http://amzn.com/1490523995}
}
\keyword{dependence}

