% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coordinates.R
\name{por_transformation_quat}
\alias{por_transformation_quat}
\alias{geographical_to_PoR_quat}
\alias{PoR_to_geographical_quat}
\title{Conversion between PoR to geographical coordinate system using quaternions}
\usage{
geographical_to_PoR_quat(x, PoR)

PoR_to_geographical_quat(x, PoR)
}
\arguments{
\item{x, PoR}{two-column vectors containing the lat and lon coordinates}
}
\value{
two-element numeric vector
}
\description{
Helper function for the transformation from PoR to geographical coordinate
system or vice versa
}
\keyword{internal}
