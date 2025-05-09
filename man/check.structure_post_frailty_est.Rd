% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check.R
\name{check.structure_post_frailty_est}
\alias{check.structure_post_frailty_est}
\title{Check Structure of Posterior Frailty Estimates}
\usage{
check.structure_post_frailty_est(post_frailty_est, n_intervals, n_centres)
}
\arguments{
\item{post_frailty_est}{Posterior frailty estimates S3 object of class 'PFE.AdPaik', composed of three elements:
\itemize{
\item 'alpha': posterior frailty estimates for \eqn{\alpha_j, \forall j}. It is a vector of length equal to the number of centres.
\item 'eps': posterior frailty estimates for \eqn{\epsilon_{jk}, \forall j,k}. It is a matrix of dimension (n_centres, n_intervals).
\item 'Z': posterior frailty estimates for \eqn{Z_{jk} = \alpha_j + \epsilon_{jk}, \forall j,k}. It is a matrix of dimension (n_centres, n_intervals)
}}

\item{n_intervals}{Number of intervals of the time-domain}

\item{n_centres}{Number of centres/clusters.}
}
\value{
An error if any condition is not satisfied.
}
\description{
The function controls that the structure of the 'Posterior Frailty Estimates' coincides with the theoretical one.
}
\keyword{internal}
