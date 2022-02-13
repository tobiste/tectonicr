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

#' @title Rotation Matrix
#'
#' @description Calculates the rotation matrix using the rotation axis and the
#' angle of rotation
#' @param n Rotation axis (in Cartesian coordinates). Can be a vector of three
#' numbers or a matrix of 3 columns (x, y, z)
#' @param alpha Rotation angle in degrees
#' @return \code{matrix}
#' @export
#' @examples
#' w <- c(0, 1, 0)
#' rotation_matrix(w, 90)
rotation_matrix <- function(n, alpha) {
  stopifnot(is.numeric(n) & is.numeric(alpha))

  n <- n / sqrt(sum(n^2)) # unit vector
  R <- matrix(nrow = 3, ncol = 3)
  R[1, 1] <- n[1]^2 * (1 - cosd(alpha)) + cosd(alpha)
  R[1, 2] <- n[1] * n[2] * (1 - cosd(alpha)) - n[3] * sind(alpha)
  R[1, 3] <- n[1] * n[3] * (1 - cosd(alpha)) + n[2] * sind(alpha)
  R[2, 1] <- n[2] * n[1] * (1 - cosd(alpha)) + n[3] * sind(alpha)
  R[2, 2] <- n[2]^2 * (1 - cosd(alpha)) + cosd(alpha)
  R[2, 3] <- n[2] * n[3] * (1 - cosd(alpha)) - n[1] * sind(alpha)
  R[3, 1] <- n[3] * n[1] * (1 - cosd(alpha)) - n[2] * sind(alpha)
  R[3, 2] <- n[3] * n[2] * (1 - cosd(alpha)) + n[1] * sind(alpha)
  R[3, 3] <- n[3]^2 * (1 - cosd(alpha)) + cosd(alpha)
  return(R)
}

#' @title Longitude Correction
#'
#' @description Corrects the longitude value to values between -180 and +180
#' degree
#' @param longitude Longitude(s) in degree
#' @return Number
#' @export
#' @examples
#' longitude_modulo(-361)
longitude_modulo <- function(longitude) {
  # longitude.mod <- (longitude %% 360 + 540) %% 360 - 180
  (longitude + 540) %% 360 - 180
}

#' @title Coordinate Transformations
#'
#' @description Converts vector between Cartesian and geographical coordinate
#' systems
#' @param n Cartesian coordinates (x, y, z) as vector
#' @param p Geographical coordinates (latitude, longitude) as vector
#' @return Functions return a (2- or 3-dimensional) vector representing a
#' point in the requested coordinate system.
#' @examples
#' n <- c(1, -2, 3)
#' cartesian_to_geographical(n)
#' p <- c(50, 10)
#' geographical_to_cartesian(p)
#' @name coordinates
NULL

#' @rdname coordinates
#' @export
cartesian_to_geographical <- function(n) {
  stopifnot(is.numeric(n) & length(n) == 3)
  r <- sqrt(n[1]^2 + n[2]^2 + n[3]^2)
  lat <- rad2deg(asin(n[3] / r))
  lon <- rad2deg(atan2(n[2], n[1]))
  if (lat <= -90) {
    lat <- 180 + lat
  }
  if (lat >= 90) {
    lat <- 180 - lat
  }
  lon <- longitude_modulo(lon)
  return(c(lat, lon))
}


#' @rdname coordinates
#' @export
geographical_to_cartesian <- function(p) {
  stopifnot(is.numeric(p) & length(p) == 2)
  x <- c()
  x[1] <- cosd(p[1]) * cosd(p[2])
  x[2] <- cosd(p[1]) * sind(p[2])
  x[3] <- sind(p[1])
  return(x)
}

#' @title Absolute Plate Velocity
#'
#' @description Calculates the absolute angular velocity of plate motion
#' @param w Angular velocity or rate or angle of rotation
#' @param alpha Angular distance to Euler pole or small circle around Euler pole
#' @param r Radius. Default is Earth's radius (6371.00887714 km)
#' @return Number
#' @export
#' @examples
#' abs_vel(0.21, 0)
#' abs_vel(0.21, 45)
#' abs_vel(0.21, 90)
abs_vel <- function(w, alpha, r = 6371.00887714) {
  stopifnot(is.numeric(w) & is.numeric(alpha))
  v <- deg2rad(w) * r * sind(alpha)
  return(v)
}
