% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{PoR_map}
\alias{PoR_map}
\title{Map of data in Pole of Rotation reference frame}
\usage{
PoR_map(
  x,
  PoR,
  pb = NULL,
  type = c("none", "in", "out", "right", "left"),
  show.deviation = FALSE,
  ...
)
}
\arguments{
\item{x, pb}{\code{sf} objects of the data points and the plate
boundary geometries in the geographical coordinate system}

\item{PoR}{Pole of Rotation. \code{"data.frame"} or object of class
\code{"euler.pole"}
containing the geographical coordinates of the Pole of Rotation}

\item{type}{Character. Type of plate boundary (optional). Can be
\code{"out"}, \code{"in"}, \code{"right"}, or
\code{"left"} for outward, inward, right-lateral, or left-lateral
moving plate boundaries, respectively. If \code{"none"} (the default), only
the PoR-equivalent azimuth is returned.}

\item{show.deviation}{logical.
Whether the data should be color-coded according to the deviation from the
prediction, or according to the stress regime? Is ignored if \code{type=='none'}.}

\item{...}{optional arguments passed to \code{\link[=tectonicr.colors]{tectonicr.colors()}}}
}
\value{
plot
}
\description{
Transforms the spatial data and azimuths into the PoR reference frame and
shows them in a map
}
\examples{
data("nuvel1")
na_pa <- subset(nuvel1, nuvel1$plate.rot == "na")

data("plates")
plate_boundary <- subset(plates, plates$pair == "na-pa")

data("san_andreas")
PoR_map(san_andreas, PoR = na_pa, pb = plate_boundary, type = "right", show.deviation = TRUE)
}
\seealso{
\code{\link[=PoR_shmax]{PoR_shmax()}}, \code{\link[=axes]{axes()}}, \code{\link[=tectonicr.colors]{tectonicr.colors()}}
}
