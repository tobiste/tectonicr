#' @title Azimuth Between two Points
#'
#' @description Calculate initial bearing (or forward azimuth/direction) to go
#' from point `a` to point `b` following great circle arc on a
#' sphere.
#'
#' @param lat_a,lat_b Numeric. Latitudes of a and b (in degrees).
#' @param lon_a,lon_b Numeric. Longitudes of a and b (in degrees).
#' @param gc character. How is the great circle angle computated. Either
#' `"orthodrome"` for the spherical law of cosines, "`"haversine"` for the
#' haversine formula (the default), or `"vincenty"` for the  Vincenty formula.
#' @details [get_azimuth_tan()] is based on the spherical law of tangents and is
#' computational faster. This formula is for the initial bearing (sometimes referred to as
#' forward azimuth) which if followed in a straight line along a great circle
#' arc will lead from the start point `a` to the end point `b`.
#' \deqn{\theta = \arctan2 (\sin \Delta\lambda
#' \cos\psi_2, \cos\psi_1 \sin\psi_1-\sin\psi_1 \cos\psi_2 \cos\Delta\lambda)}
#' where  \eqn{\psi_1, \lambda_1} is the start point, \eqn{\psi_2},
#' \eqn{\lambda_2} the end point (\eqn{\Delta\lambda} is the difference in
#' longitude)
#'
#' [get_azimuth_sin()] uses the Law of sines and the great circle distance.
#' \deqn{
#'  \theta = \asin( \frac{\sin\Delta\phi * \cos\lambda_2}{\sin \gamma} )
#' }
#' where  \eqn{\gamma} is the great-circle distance between the points.
#'
#' [get_azimuth_half()] uses the half-angle formula
#' \deqn{\theta = 2 \arctan \sqrt(\frac{\sin(s - (\pi/2 -\lambda_1)) \sin \gamma}{\sin s \sin(s - (\pi/2 - \lambda_2))})
#' }
#' where \eqn{s = \frac{1}{2} (\pi/2 -\lambda_1 + \pi/2 -\lambda_2 + \gamma)}.
#'
#' @seealso [orthodrome()], [haversine()], [vincenty()]
#' @references \url{http://www.movable-type.co.uk/scripts/latlong.html}
#' @return Azimuth in degrees
#' @name get_azimuth_alt
#' @examples
#' berlin <- c(52.517, 13.4) # Berlin
#' tokyo <- c(35.7, 139.767) # Tokyo
#' get_azimuth(berlin[1], berlin[2], tokyo[1], tokyo[2])
#' get_azimuth_tan(berlin[1], berlin[2], tokyo[1], tokyo[2])
#' get_azimuth_sin(berlin[1], berlin[2], tokyo[1], tokyo[2])
#' get_azimuth_half(berlin[1], berlin[2], tokyo[1], tokyo[2])
NULL

#' @export
#' @rdname get_azimuth_alt
get_azimuth_tan <- function(lat_a, lon_a, lat_b, lon_b) {
  # stopifnot(is.numeric(lat_a), is.numeric(lat_b), is.numeric(lon_a), is.numeric(lon_b))

  # convert deg into rad
  phi1 <- pi / 180 * lat_a
  phi2 <- pi / 180 * lat_b

  d.lambda <- (lon_b - lon_a) * (pi / 180)
  cos_phi_2 <- cos(phi2)

  y <- sin(d.lambda) * cos_phi_2
  x <- cos(phi1) * sin(phi2) -
    sin(phi1) * cos_phi_2 * cos(d.lambda)
  theta <- atan2d(y, x)

  # Normalize result to a compass bearing (0-360)
  (theta + 360) %% 360
}

#' @export
#' @rdname get_azimuth_alt
get_azimuth_sin <- function(lat_a, lon_a, lat_b, lon_b, gc = c("haversine", "orthodrome", "vincenty")) {
  # stopifnot(is.numeric(lat_a), is.numeric(lat_b), is.numeric(lon_a), is.numeric(lon_b))

  # convert deg into rad
  x <- c(lat_a, lon_a, lat_b, lon_b) * pi / 180

  # dlon <- lon2 - lon1
  # acos(sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(dlon))

  gc <- match.arg(gc)
  if (gc == "haversine") {
    gamma <- haversine(x[1], x[2], x[3], x[4])
  } else if (gc == "vincenty") {
    gamma <- vincenty(x[1], x[2], x[3], x[4])
  } else {
    gamma <- orthodrome(x[1], x[2], x[3], x[4])
  }

  if (cos(gamma) > sin(x[1]) * sin(x[3])) {
    y <- sin(x[4] - x[2]) * cos(x[3])
  } else {
    y <- sin(pi + x[2] - x[4]) * cos(x[3])
  }
  theta <- asind(y / sin(gamma))

  # Normalize result to a compass bearing (0-360)
  (theta + 360) %% 360
}

#' @export
#' @rdname get_azimuth_alt
get_azimuth_sin_fast <- function(lat_a, lon_a, lat_b, lon_b) {
  # convert deg into rad
  x <- c(lat_a, lon_a, lat_b, lon_b) * pi / 180

  dlon <- x[4] - x[2]
  cos_x3 <- cos(x[3])
  sin_x1__sin_x3 <- sin(x[1]) * sin(x[3])
  cos_gamma <- sin_x1__sin_x3 + cos_x3 * cos(x[1]) * cos(dlon)

  if (cos_gamma > sin_x1__sin_x3) {
    y <- sin(dlon) * cos_x3
  } else {
    y <- sin(pi - dlon) * cos_x3
  }
  theta <- asind(y / sin(acos(cos_gamma)))

  # Normalize result to a compass bearing (0-360)
  (theta + 360) %% 360
}


#' @export
#' @rdname get_azimuth_alt
get_azimuth_half <- function(lat_a, lon_a, lat_b, lon_b, gc = c("haversine", "orthodrome", "vincenty")) {
  # Half-angle formula
  # stopifnot(is.numeric(lat_a), is.numeric(lat_b), is.numeric(lon_a), is.numeric(lon_b))

  # convert deg into rad
  x <- c(lat_a, lon_a, lat_b, lon_b) * pi / 180

  gc <- match.arg(gc)
  if (gc == "haversine") {
    c <- haversine(x[1], x[2], x[3], x[4])
  } else if (gc == "vincenty") {
    c <- vincenty(x[1], x[2], x[3], x[4])
  } else {
    c <- orthodrome(x[1], x[2], x[3], x[4])
  }

  b <- pi / 2 - x[1]
  a <- pi / 2 - x[3]

  s <- (a + b + c) / 2

  x <- sin(s - b) * sin(s - c)
  y <- sin(s) * sin(s - a)

  theta <- 2 * atand(sqrt(x / y))

  # Normalize result to a compass bearing (0-360)
  (theta + 360) %% 360
}


# Rotation matrices ####

#' Basic rotation matrices
#'
#' Helper functions for constructing the matrices for rotations about the
#' x-, y-, and z angles with the angle x
#'
#' @param x angle (in degree)
#' @return matrix
#' @name transform_matrices
NULL

#' @rdname transform_matrices
rotmat_x <- function(x) {
  x <- deg2rad(as.numeric(x))
  rbind(
    c(1, 0, 0),
    c(0, cos(x), -sin(x)),
    c(0, sin(x), cos(x))
  )
}

#' @rdname transform_matrices
rotmat_y <- function(x) {
  x <- deg2rad(as.numeric(x))
  rbind(
    c(cos(x), 0, sin(x)),
    c(0, 1, 0),
    c(-sin(x), 0, cos(x))
  )
}

#' @rdname transform_matrices
rotmat_z <- function(x) {
  x <- deg2rad(as.numeric(x))
  rbind(
    c(cos(x), -sin(x), 0),
    c(sin(x), cos(x), 0),
    c(0, 0, 1)
  )
}


