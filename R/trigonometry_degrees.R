#' @title Degrees to Radians
#'
#' @description Transforms between angles in degrees and radians.
#'
#' @param deg	 (array of) angles in degrees.
#' @param rad (array of) angles in radians.
#' @return The angle in degrees or radians.
#' @source \code{\link[pracma]{deg2rad}}, \code{\link[pracma]{rad2deg}} from
#' package "pracma".
#' @examples
#' deg2rad(c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90))
#' rad2deg(seq(-pi / 2, pi / 2, length = 19))
#' @name angle-conversion
NULL

#' @rdname angle-conversion
#' @export
rad2deg <- function(rad) {
  stopifnot(is.numeric(rad))
  (deg <- rad / (pi / 180))
}
#' @rdname angle-conversion
#' @export
deg2rad <- function(deg) {
  stopifnot(is.numeric(deg))
  (rad <- (pi / 180) * deg)
}

#' @title Trigonometric Functions in Degrees
#' @description Trigonometric functions expecting input in degrees, not radians.
#'
#' @param x,x1,x2 Numeric or complex vectors.
#' @return Returns a scalar or vector of numeric values.
#' @source \code{\link[pracma]{sind}}, \code{\link[pracma]{cosd}}, ...  from
#' package "pracma".
#' @examples
#' x <- seq(-3, 7, by = 1 / 8)
#' tx <- cbind(x, cos(pi * x), cospi(x), sin(pi * x), sinpi(x),
#'   tan(pi * x), tanpi(x),
#'   deparse.level = 2
#' )
#' op <- options(digits = 4, width = 90) # for nice formatting
#' head(tx)
#' tx[(x %% 1) %in% c(0, 0.5), ]
#' options(op)
#' @name trigon
NULL

#' @rdname trigon
#' @export
sind <- function(x) {
  sinpi(x / 180)
}
#' @rdname trigon
#' @export
cosd <- function(x) {
  cospi(x / 180)
}
#' @rdname trigon
#' @export
tand <- function(x) {
  sinpi(x / 180) / cospi(x / 180)
}
#'
#' @rdname trigon
#' @export
asind <- function(x) {
  asin(x) * 180 / pi
}
#' @rdname trigon
#' @export
acosd <- function(x) {
  acos(x) * 180 / pi
}
#' @rdname trigon
#' @export
atand <- function(x) {
  atan(x) * 180 / pi
}
#'
#' @rdname trigon
#' @export
atan2d <- function(x1, x2) {
  atan2(x1, x2) * 180 / pi
}
