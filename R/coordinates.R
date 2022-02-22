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
  c(lat, lon)
}


#' @rdname coordinates
#' @export
geographical_to_cartesian <- function(p) {
  stopifnot(is.numeric(p) & length(p) == 2)
  x <- c()
  x[1] <- cosd(p[1]) * cosd(p[2])
  x[2] <- cosd(p[1]) * sind(p[2])
  x[3] <- sind(p[1])
  x
}
