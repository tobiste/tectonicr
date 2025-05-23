% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_shmax.R
\name{data2PoR}
\alias{data2PoR}
\title{Transforms coordinates and azimuths into PoR coordinates system}
\usage{
data2PoR(x, PoR)
}
\arguments{
\item{x}{\code{sf} object or a \code{data.frame} containing the coordinates of the
point(s) (\code{lat}, \code{lon} columns). \code{x} must contain the direction of
\eqn{\sigma_{Hmax}}{SHmax} as column \code{azi}, its standard deviation
(column \code{unc}) is optional).}

\item{PoR}{\code{data.frame} or object of class \code{euler.pole} containing the
geographical coordinates of the Eule pole.}
}
\value{
\code{sf} object in PoR CRS with additional columns \code{lon.PoR},
\code{lat.PoR}, and \code{azi.PoR}
}
\description{
Convenience function to add PoR coordinates and PoR azimuths to data
}
\examples{
por <- subset(nuvel1, nuvel1$plate.rot == "na")
data2PoR(san_andreas, por)
}
