% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rotation.R
\name{abs_vel}
\alias{abs_vel}
\title{Absolute Plate Velocity}
\usage{
abs_vel(w, alpha, r = earth_radius())
}
\arguments{
\item{w}{Angular velocity or rate or angle of rotation}

\item{alpha}{Angular distance to Euler pole or small circle around Euler pole}

\item{r}{Radius. Default is WGS84 Earth's radius (6371.009 km)}
}
\value{
numeric (unit of velocity: km/Myr)
}
\description{
Calculates the absolute angular velocity of plate motion
}
\examples{
abs_vel(0.21, 0)
abs_vel(0.21, 45)
abs_vel(0.21, 90)
}
\seealso{
\code{\link[=earth_radius]{earth_radius()}}
}
