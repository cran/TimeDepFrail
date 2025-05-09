% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check.R
\name{check.poslegend}
\alias{check.poslegend}
\title{Check Correctness of Legend Position}
\usage{
check.poslegend(pos_legend)
}
\arguments{
\item{pos_legend}{Either a numerical vector of length 2, with the x and y coordinates, or a string with the exact position.}
}
\value{
An error if any condition is not satisfied.
}
\description{
The function controls that the provided position of the legend is correct. It can be either
a vector of length 2, giving the x and y coordinates, or a string, giving the exact position among different possibilities.
}
\keyword{internal}
