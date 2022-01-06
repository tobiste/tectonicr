#' @title Angle between two vectors
#'
#' @description Calculates the angle between two vectors
#' @author Tobias Stephan
#' @param x Vector 1 in cartesian coordinates. Can be a vector of three numbers
#'  or a matrix of 3 columns (first one is x, second is y, thrid one is z)
#' @param y Vector 2. As above
#' @return numeric; angle in degrees
#' @export
#' @examples
#' u <- c(1, -2, 3)
#' v <- c(-2, 1, 1)
#' angle_vectors(u, v)
angle_vectors <- function(x, y) {
  angle.d <- pracma::rad2deg(
    acos(sum(x * y) / (sqrt(sum(x * x)) * sqrt(sum(y * y))))
  )
  return(angle.d)
}

#' @title Rotation matrix
#'
#' @description Calculates the rotation matrix using the rotation axis and the angle of rotation
#' @author Tobias Stephan
#' @param n Rotation axis (in cartesian coordinates). Can be a vector of three numbers
#'  or a matrix of 3 columns (first one is x, second is y, thrid one is z)
#' @param alpha Rotation angle in degrees
#' @return matrix
#' @export
#' @examples
#' w <- c(0, 1, 0)
#' rotation_matrix(w, 90)
rotation_matrix <- function(n, alpha) {
  n <- ppls::normalize.vector(n) # unit vector
  R <- matrix(nrow = 3, ncol = 3)
  R[1, 1] <- n[1]^2 * (1 - pracma::cosd(alpha)) + pracma::cosd(alpha)
  R[1, 2] <- n[1] * n[2] * (1 - pracma::cosd(alpha)) - n[3] * pracma::sind(alpha)
  R[1, 3] <- n[1] * n[3] * (1 - pracma::cosd(alpha)) + n[2] * pracma::sind(alpha)
  R[2, 1] <- n[2] * n[1] * (1 - pracma::cosd(alpha)) + n[3] * pracma::sind(alpha)
  R[2, 2] <- n[2]^2 * (1 - pracma::cosd(alpha)) + pracma::cosd(alpha)
  R[2, 3] <- n[2] * n[3] * (1 - pracma::cosd(alpha)) - n[1] * pracma::sind(alpha)
  R[3, 1] <- n[3] * n[1] * (1 - pracma::cosd(alpha)) - n[2] * pracma::sind(alpha)
  R[3, 2] <- n[3] * n[2] * (1 - pracma::cosd(alpha)) + n[1] * pracma::sind(alpha)
  R[3, 3] <- n[3]^2 * (1 - pracma::cosd(alpha)) + pracma::cosd(alpha)
  return(R)
}

#' @title Longitude correction
#'
#' @description Corrects the longitude value to values between -180 and +180 degree
#' @param longitude longitude(s) in degree
#' @return number
#' @export
#' @examples
#' longitude_modulo(-361)
longitude_modulo <- function(longitude) {
  longitude.mod <- (longitude %% 360 + 540) %% 360 - 180
  return(longitude.mod)
}

#' @title Cartesian to geographical coordinates
#'
#' @description Converts vector from cartesian into geographical coordinates
#' @param n Vector of three numbers (x, y, z)
#' @return vector of two numbers (latitude, longitude)
#' @export
#' @examples
#' u <- c(1, -2, 3)
#' cartesian_to_geographical(u)
cartesian_to_geographical <- function(n) {
  r <- sqrt(n[1]^2 + n[2]^2 + n[3]^2)
  lat <- pracma::rad2deg(asin(n[3] / r))
  lon <- pracma::rad2deg(atan2(n[2], n[1]))
  if (lat <= -90) {
    lat <- 180 + lat
  }
  if (lat >= 90) {
    lat <- 180 - lat
  }
  lon <- longitude_modulo(lon)
  return(c(lat, lon))
}

#' @title Geographical to cartesian coordinates
#'
#' @description Converts vector from geographical into cartesian  coordinates
#' @param p vector of two numbers (latitude, longitude)
#' @param r radius. default is 1.
#' @return Vector of three numbers (x, y, z)
#' @export
#' @examples
#' u <- c(50, 10)
#' geographical_to_cartesian(u)
geographical_to_cartesian <- function(p, r = 1) {
  if (missing(r)) {
    r <- 1
  }
  x <- c()
  x[1] <- r * pracma::cosd(p[1]) * pracma::cosd(p[2])
  x[2] <- r * pracma::cosd(p[1]) * pracma::sind(p[2])
  x[3] <- r * pracma::sind(p[1])
  return(x)
}

#' @title Absolute plate velocity
#'
#' @description Calculates the absolute angular velocity of plate motion
#' @param w Angular velocity or rate or angle of rotation
#' @param alpha Angular distance to Euler pole or smallcircle around Euler pole
#' @param r radius. Default is Earth's radius (r=6371.00887714)
#' @return number
#' @export
#' @examples
#' abs_vel(0.21, 0)
#' abs_vel(0.21, 45)
#' abs_vel(0.21, 90)
abs_vel <- function(w, alpha, r = 6371.00887714) {
  v <- pracma::deg2rad(w) * r * pracma::sind(alpha)
  return(v)
}
