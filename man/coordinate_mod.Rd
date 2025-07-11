% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coordinates.R
\name{coordinate_mod}
\alias{coordinate_mod}
\alias{longitude_modulo}
\alias{latitude_modulo}
\title{Coordinate Correction}
\usage{
longitude_modulo(x)

latitude_modulo(x)
}
\arguments{
\item{x}{Longitude(s) or latitude(s) in degrees}
}
\value{
numeric
}
\description{
Corrects the longitudes or latitudes to value between -180.0 and
180.0 or -90 and 90 degree
}
\examples{
longitude_modulo(-361 + 5 * 360) # -1
latitude_modulo(-91 + 5 * 180) # 89
}
