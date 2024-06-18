% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{rose_stats}
\alias{rose_stats}
\title{Show Average Direction and Spread in Rose Diagram}
\usage{
rose_stats(
  x,
  weights = NULL,
  axial = TRUE,
  avg = c("mean", "median"),
  spread = c("CI", "sd", "IQR"),
  avg.col = "#85112AFF",
  avg.lty = 2,
  avg.lwd = 1.5,
  spread.col = ggplot2::alpha("#85112AFF", 0.2),
  spread.border = FALSE,
  spread.lty = NULL,
  spread.lwd = NULL,
  add = TRUE,
  ...
)
}
\arguments{
\item{x}{Data to be plotted. A numeric vector containing angles (in degrees).}

\item{weights}{Optional vector of numeric weights associated with x.}

\item{axial}{Logical. Whether data are uniaxial (\code{axial=FALSE})
or biaxial (\code{TRUE}, the default).}

\item{avg}{character. The average estimate for x. Either the circular mean
(\code{"mean"}, the default) or the circular Quasi Median (\code{"median"})}

\item{spread}{character. The measure of spread to be plotted as a fan.
Either the 95\% confidence interval (\code{"CI"}, the default), the circular
standard deviation (\code{"sd"}), or the Quasi interquartile range on the circle
(\code{"IQR"}). \code{NULL} if no fan should be drawn.}

\item{avg.col}{color for the average line}

\item{avg.lty}{line type of the average line}

\item{avg.lwd}{line width of the average line}

\item{spread.col}{color of the spread fan}

\item{spread.border}{logical. Whether to draw a border of the fan or not.}

\item{spread.lty}{line type of the spread fan's border}

\item{spread.lwd}{line width of the spread fan's border}

\item{add}{logical.}

\item{...}{optional arguments to \code{rose_baseplot()} if add is \code{FALSE}.}
}
\value{
No return value, called for side effects
}
\description{
Adds the average direction (and its spread) to an existing rose diagram.
}
\examples{
data("san_andreas")
rose(san_andreas$azi, weights = 1 / san_andreas$unc, muci = FALSE)
rose_stats(san_andreas$azi, weights = 1 / san_andreas$unc, avg = "median", spread = "IQR")
}
\seealso{
\code{\link[=rose]{rose()}} for plotting the rose diagram, and
\code{\link[=circular_mean]{circular_mean()}}, \code{\link[=circular_median]{circular_median()}}, \code{\link[=confidence_interval]{confidence_interval()}},
\code{\link[=circular_sd]{circular_sd()}}, \code{\link[=circular_IQR]{circular_IQR()}} for statistical parameters.
}