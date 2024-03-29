% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coordinates.R
\name{PoR_crs}
\alias{PoR_crs}
\title{PoR coordinate reference system}
\usage{
PoR_crs(x)
}
\arguments{
\item{x}{\code{"data.frame"} or \code{"euler.pole"} object containing the
geographical coordinates of the Euler pole}
}
\value{
Object of class \code{crs}
}
\description{
Create the reference system transformed in Euler pole
coordinate
}
\details{
The PoR coordinate reference system is oblique transformation of the
geographical coordinate system with the Euler pole coordinates being the the
translation factors.
}
\examples{
data("nuvel1")
por <- subset(nuvel1, nuvel1$plate.rot == "na") # North America relative to Pacific plate
PoR_crs(por)
}
\seealso{
\code{\link[sf:st_crs]{sf::st_crs()}}
}
