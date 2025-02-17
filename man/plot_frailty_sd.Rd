% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot_frailty_sd}
\alias{plot_frailty_sd}
\title{Plot for the Frailty Standard Deviation or Variance}
\usage{
plot_frailty_sd(
  result,
  flag_variance = FALSE,
  flag_sd_external = FALSE,
  frailty_sd = NULL,
  xlim = c(min(result$TimeDomain), max(result$TimeDomain)),
  ylim = NULL,
  xlab = "Time",
  ylab = "Values",
  main = NULL,
  pch = 21,
  color_bg = "blue",
  cex_points = 0.7
)
}
\arguments{
\item{result}{An S3 object of class 'AdPaik', returned by the main model call 'AdPaikModel(...)'.}

\item{flag_variance}{A boolean flag indicating whether to plot the frailty variance (\code{TRUE}) or the frailty standard deviation (\code{FALSE}). Default is \code{FALSE}.}

\item{flag_sd_external}{A logical flag indicating whether the user is providing an external frailty standard deviation vector.}

\item{frailty_sd}{A numerical vector representing the evaluated frailty standard deviation, with length equal to the number of time-domain intervals.
Its elements must be non-negative. Default is \code{NULL}.}

\item{xlim}{A numeric vector specifying the range for the x-axis (intervals). If NULL, default is set to the interval min-max of the time-domain.}

\item{ylim}{A numeric vector specifying the range for the y-axis (intervals). If NULL, default is 0 to the maximum value of the frailty variance/standard deviation.}

\item{xlab}{A string for the x-axis label. Default is \code{'Intervals'}.}

\item{ylab}{A string for the y-axis label. Default is \code{'Values'}.}

\item{main}{A string for the plot title. Default title is \code{'Frailty Standard Deviation'} or \code{'Frailty Variance'} according to the produced plot (flag_variance).}

\item{pch}{A numeric or character symbol used for plotting the frailty standard deviation values. Default is a dot (\code{21}).}

\item{color_bg}{A string specifying the color used for the symbols. Default is \code{'blue'}.}

\item{cex_points}{A numeric value indicating the size of the plotting symbols. Default is \code{0.7}.}
}
\value{
A plot displaying either the frailty standard deviation or variance across the specified intervals.
}
\description{
This function generates a plot of either the frailty standard deviation or the frailty variance for the intervals in the time-domain.
}
\details{
The plot represents the values of the frailty standard deviation or variance for each time interval (represented by its mid point).
It connects these points to illustrate the trend of the chosen metric.

This function supports two modes of operation:
\itemize{
\item Plotting the frailty standard deviation or variance retrieved from the main model (contained in the S3 object of class 'AdPaik').
\item Plotting a user-provided vector of frailty standard deviations, which can be computed using the method \code{frailty.sd}. This allows for flexibility in analysis without re-optimizing the log-likelihood function. For instance, users can compare frailty standard deviations computed with different model specifications (e.g., including only time-dependent terms).
}

The output will differentiate between these two cases, ensuring the correct values are plotted regardless of the source.
}
\examples{
# Import data
data(data_dropout)

# Define the variables needed for the model execution
eps_paik <- 1e-10
categories_range_min <- c(-8, -2, eps_paik, eps_paik, eps_paik)
categories_range_max <- c(-eps_paik, 0.4, 1 - eps_paik, 1, 10)
time_axis <- c(1.0, 1.4, 1.8, 2.3, 3.1, 3.8, 4.3, 5.0, 5.5, 5.8, 6.0)
formula <- time_to_event ~ Gender + CFUP + cluster(group)

# Call the main model function
\donttest{
result <- AdPaikModel(formula, data_dropout, time_axis, categories_range_min, categories_range_max)

plot_frailty_sd(result)
}
}
