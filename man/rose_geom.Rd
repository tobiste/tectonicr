% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{rose_geom}
\alias{rose_geom}
\alias{rose_line}
\alias{rose_fan}
\title{Direction Lines and Fans in Circular Diagram}
\usage{
rose_line(x, radius = 1, axial = TRUE, add = TRUE, ...)

rose_fan(x, d, radius = 1, axial = TRUE, add = TRUE, ...)
}
\arguments{
\item{x}{angles in degrees}

\item{radius}{of the plotted circle}

\item{axial}{Logical. Whether \code{x} are uniaxial (\code{axial=FALSE})
or biaxial (\code{TRUE}, the default).}

\item{add}{logical. Add to existing plot?}

\item{...}{optional arguments passed to \code{\link[graphics:segments]{graphics::segments()}} or
\code{\link[graphics:polygon]{graphics::polygon()}}}

\item{d}{width of a fan (in degrees)}
}
\value{
No return value, called for side effects
}
\description{
Direction Lines and Fans in Circular Diagram
}
\examples{
angles <- c(0, 10, 45)
radius <- c(.7, 1, .2)
lwd <- c(2, 1, .75)
col <- c(1, 2, 3)
rose_line(angles, radius = radius, axial = FALSE, add = FALSE, lwd = lwd, col = col)
}