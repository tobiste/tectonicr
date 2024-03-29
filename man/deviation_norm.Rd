% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_shmax.R
\name{deviation_norm}
\alias{deviation_norm}
\title{Normalize Angle Between Two Directions}
\usage{
deviation_norm(x)
}
\arguments{
\item{x}{numeric vector containing angles in degrees}
}
\value{
numeric vector, acute angles between two directions, i.e. values
between 0 and 90°
}
\description{
Normalizes the angle between two directions to the acute angle
in between, i.e. angles between 0 and 90°
}
\examples{
deviation_norm(91)
deviation_norm(c(-91, NA, 23497349))
}
\author{
Tobias Stephan
}
