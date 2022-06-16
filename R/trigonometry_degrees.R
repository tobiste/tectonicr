#' @title Degrees to Radians
#'
#' @description Helper functions to transform between angles in degrees and
#' radians.
#'
#' @param deg	 (array of) angles in degrees.
#' @param rad (array of) angles in radians.
#' @return The angle in degrees or radians.
#' @examples
#' deg2rad(seq(-90, 90, 15))
#' rad2deg(seq(-pi / 2, pi / 2, length = 13))
#' @name angle-conversion
NULL

#' @rdname angle-conversion
#' @export
rad2deg <- function(rad) {
  stopifnot(is.numeric(rad))
  rad * 180 / pi
}
#' @rdname angle-conversion
#' @export
deg2rad <- function(deg) {
  stopifnot(is.numeric(deg))
  deg * pi / 180
}

#' @title Trigonometric Functions in Degrees
#' @description Trigonometric functions expecting input in degrees, not radians.
#'
#' @param x,x1,x2 Numeric or complex vectors.
#' @return Returns a scalar or vector of numeric values.
#' @keywords internal
#' @name trigon
NULL

#' @rdname trigon
sind <- function(x) {
  sinpi(x / 180)
}
#' @rdname trigon
cosd <- function(x) {
  cospi(x / 180)
}
#' @rdname trigon
tand <- function(x) {
  sinpi(x / 180) / cospi(x / 180)
}
#'
#' @rdname trigon
asind <- function(x) {
  asin(x) * 180 / pi
}
#' @rdname trigon
acosd <- function(x) {
  acos(x) * 180 / pi
}
#' @rdname trigon
atand <- function(x) {
  atan(x) * 180 / pi
}
#'
#' @rdname trigon
atan2d <- function(x1, x2) {
  atan2(x1, x2) * 180 / pi
}

#' Quadrant-specific inverse of the tangent
#'
#' returns the quadrant specific inverse of the tangent
#'
#' @param x,y dividend and divisor that comprise the sum of sines and cosines,
#' respectively.
#'
#' @references Jammalamadaka, S. Rao, and Ambar Sengupta. Topics in circular
#' statistics. Vol. 5. world scientific, 2001.
#' @name spec_atan
NULL

#' @rdname spec_atan
#' @export
atan2_spec <- function(x, y) {
  if (y > 0 & x >= 0) {
    atan(x / y)
  } else if (y == 0 & x > 0) {
    pi / 2
  } else if (y < 0) {
    atan(x / y + pi)
  } else if (y > 0 & x < 0) {
    atan(x / y) + 2 * pi
  } else if (y == 0 & x == 0) {
    Inf
  }
}

#' @rdname spec_atan
#' @export
atan2d_spec <- function(x, y) {
  atan2_spec(x, y) * 180 / pi
}


#' @title Angle Between Two Vectors
#'
#' @description Calculates the angle between two vectors
#' @param x,y Vectors in Cartesian coordinates. Can be vectors of three numbers
#'  or a matrix of 3 columns (x, y, z)
#' @return Numeric; angle in degrees
#' @export
#' @examples
#' u <- c(1, -2, 3)
#' v <- c(-2, 1, 1)
#' angle_vectors(u, v)
angle_vectors <- function(x, y) {
  stopifnot(is.numeric(x) & is.numeric(y))
  if (length(x) == length(y)) {
    angle.d <- rad2deg(
      acos(sum(x * y) / (sqrt(sum(x * x)) * sqrt(sum(y * y))))
    )
    return(angle.d)
  }
}
