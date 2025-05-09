% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check.R
\name{check.structure_post_frailty_var}
\alias{check.structure_post_frailty_var}
\title{Check Structure of Posterior Frailty Variances}
\usage{
check.structure_post_frailty_var(post_frailty_var, n_intervals, n_centres)
}
\arguments{
\item{post_frailty_var}{Posterior frailty variances S3 object of class 'PFV.AdPaik', composed of three elements:
\itemize{
\item 'alphaVar': posterior frailty variance for \eqn{\alpha_j, \forall j}. It is a vector of length equal to the number of centres.
\item 'epsVar': posterior frailty variance for \eqn{\epsilon_{jk}, \forall j,k}. It is a matrix of dimension (n_centres, n_intervals).
\item 'ZVar': posterior frailty variance for \eqn{Z_{jk} = \alpha_j + \epsilon_{jk}, \forall j,k}. It is a matrix of dimension (n_centres, n_intervals)
}}

\item{n_intervals}{Number of intervals of the time-domain}

\item{n_centres}{Number of centres/clusters.}
}
\value{
An error if any condition is not satisfied.
}
\description{
The function controls that the structure of the 'Posterior Frailty Variances' coincides with the theoretical one.
}
\keyword{internal}
