% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/various.R
\name{quantise_wsm_quality}
\alias{quantise_wsm_quality}
\title{Numerical values to World Stress Map Quality Ranking}
\usage{
quantise_wsm_quality(regime)
}
\arguments{
\item{regime}{Either a string or a character vector of WSM quality ranking}
}
\value{
\code{"integer"} or vector of type \code{"integer"}
}
\description{
Assigns numeric values of the precision of each measurement to the
categorical quality ranking of the World Stress Map (A, B, C, D).
}
\examples{
quantise_wsm_quality(c("A", "B", "C", "D", NA))
data("san_andreas")
quantise_wsm_quality(san_andreas$quality)
}
\references{
Heidbach, O., Barth, A., Müller, B., Reinecker, J.,
Stephansson, O., Tingay, M., Zang, A. (2016). WSM quality
ranking scheme, database description and analysis guidelines for stress
indicator. \emph{World Stress Map Technical Report} \strong{16-01}, GFZ German Research
Centre for Geosciences. \doi{10.2312/wsm.2016.001}
}
