% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/various.R
\name{parse_wsm}
\alias{parse_wsm}
\alias{parse_wsm_quality}
\alias{quantise_wsm_quality}
\title{Numerical values to World Stress Map Quality Ranking}
\usage{
parse_wsm_quality(x)

quantise_wsm_quality(x)
}
\arguments{
\item{x}{Either a string or a character/factor vector of WSM quality ranking}
}
\value{
\code{"numeric"}. the standard deviation of stress azimuth
}
\description{
Assigns numeric values of the precision (sd.) of each measurement to the
categorical quality ranking of the World Stress Map (A, B, C, D, E, X).
}
\examples{
parse_wsm_quality(c("A", "B", "C", "D", NA, "E", "X"))
data("san_andreas")
head(parse_wsm_quality(san_andreas$quality))
}
\references{
Heidbach, O., Barth, A., M<U+00FC>ller, B., Reinecker, J.,
Stephansson, O., Tingay, M., Zang, A. (2016). WSM quality
ranking scheme, database description and analysis guidelines for stress
indicator. \emph{World Stress Map Technical Report} \strong{16-01}, GFZ German Research
Centre for Geosciences. \doi{10.2312/wsm.2016.001}
}
