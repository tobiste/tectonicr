% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coordinates.R
\name{coordinates}
\alias{coordinates}
\alias{cartesian_to_geographical}
\alias{geographical_to_cartesian}
\alias{geographical_to_spherical}
\title{Coordinate Transformations}
\usage{
cartesian_to_geographical(n)

geographical_to_cartesian(p)

geographical_to_spherical(p)
}
\arguments{
\item{n}{Cartesian coordinates (x, y, z) as vector}

\item{p}{Geographical coordinates (latitude, longitude) as vector}
}
\value{
Functions return a (2- or 3-dimensional) vector representing a
point in the requested coordinate system.
}
\description{
Converts vector between Cartesian and geographical coordinate
systems
}
\examples{
n <- c(1, -2, 3)
cartesian_to_geographical(n)
p <- c(50, 10)
geographical_to_cartesian(p)
}
\seealso{
\code{\link[=cartesian_to_spherical]{cartesian_to_spherical()}} and \code{\link[=spherical_to_cartesian]{spherical_to_cartesian()}} for
conversions to spherical coordinates
}
