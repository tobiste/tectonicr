#' @title Degrees to Radians
#'
#' @description Helper functions to transform between angles in degrees and
#' radians.
#'
#' @param deg	 (array of) angles in degrees.
#' @param rad (array of) angles in radians.
#'
#' @returns numeric. angle in degrees or radians.
#'
#' @examples
#' deg2rad(seq(-90, 90, 15))
#' rad2deg(seq(-pi / 2, pi / 2, length = 13))
#' @name angle-conversion
NULL

#' @rdname angle-conversion
#' @export
rad2deg <- function(rad) {
  # stopifnot(is.numeric(rad))
  rad * 180 / pi
}
#' @rdname angle-conversion
#' @export
deg2rad <- function(deg) {
  # stopifnot(is.numeric(deg))
  deg * pi / 180
}

#' @title Trigonometric Functions in Degrees
#'
#' @description Trigonometric functions expecting input in degrees.
#'
#' @param x,x1,x2 Numeric or complex vectors.
#'
#' @returns scalar or vector of numeric values.
#'
#' @keywords internal
#'
#' @name trigon
NULL

dir2ax <- function(x) {
  (x / 2) %% 180
}

ax2dir <- function(x) {
  (2 * x) %% 360
}


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

#' @rdname trigon
cot <- function(x) {
  1 / tan(x)
}

#' @rdname trigon
cotd <- function(x) {
  1 / tand(x)
}


#' Quadrant-specific inverse of the tangent
#'
#' Returns the quadrant specific inverse of the tangent
#'
#' @param x,y dividend and divisor that comprise the sum of sines and cosines,
#' respectively.
#'
#' @references Jammalamadaka, S. Rao, and Ambar Sengupta (2001). Topics in
#' circular statistics. Vol. 5. world scientific.
#'
#' @returns numeric.
#'
#' @name spec_atan
NULL

#' @rdname spec_atan
#' @export
atan2_spec <- function(x, y) {
  # if (y > 0 && x >= 0) {
  #   atan(x / y)
  # } else if (y == 0 && x > 0) {
  #   pi / 2
  # } else if (y < 0) {
  #   atan(x / y + pi)
  # } else if (y > 0 && x < 0) {
  #   atan(x / y) + 2 * pi
  # } else if (y == 0 && x == 0) {
  #   Inf
  # }
  # dplyr::case_when(
  #   y > 0 & x >= 0 ~ atan(x / y),
  #   y == 0 & x > 0 ~ pi / 2,
  #   y < 0 ~ atan(x / y + pi),
  #   y > 0 & x < 0 ~ atan(x / y) + 2 * pi,
  #   y == 0 & x == 0 ~ Inf
  # )
  angle <- atan2(x, y)
  ifelse(angle < 0, angle + 2 * pi, angle)
}

#' @rdname spec_atan
#' @export
atan2d_spec <- function(x, y) {
  atan2_spec(x, y) * 180 / pi
}



#' @title Angle Between Two Vectors
#'
#' @description Calculates the angle between two vectors
#'
#' @param x,y Vectors in Cartesian coordinates. Can be vectors of three numbers
#'  or a matrix of 3 columns (x, y, z)
#'
#' @returns numeric. angle in degrees
#'
#' @export
#'
#' @examples
#' u <- c(1, -2, 3)
#' v <- c(-2, 1, 1)
#' angle_vectors(u, v) # 96.26395
angle_vectors <- function(x, y) {
  stopifnot(length(x) == length(y))

  dot <- sum(x * y)
  norm_x <- sqrt(sum(x^2))
  norm_y <- sqrt(sum(y^2))

  cos_angle <- dot / (norm_x * norm_y)
  cos_angle <- pmin(pmax(cos_angle, -1), 1)

  # cos_angle <- dot / (sqrt(sum(x^2) * sum(y^2)))

  angle_deg <- acos(cos_angle) * 180 / pi

  return(angle_deg)
}


#' @keywords internal
hav <- function(x) {
  sin(x / 2)^2
}

#' @keywords internal
ahav <- function(x) {
  2 * asin(sqrt(pmin(pmax(x, 0), 1)))
}

#' Angle along great circle on spherical surface
#'
#' Smallest angle between two points on the surface of a sphere, measured along
#' the surface of the sphere
#'
#' @param lat1,lat2 numeric vector. latitudes of point 1 and 2 (in radians)
#' @param lon1,lon2 numeric vector. longitudes of point 1 and 2 (in radians)
#' @details \describe{
#' \item{\code{"orthodrome"}}{based on the spherical law of cosines}
#' \item{\code{"haversine"}}{uses haversine formula that is
#' optimized for 64-bit floating-point numbers}
#' \item{\code{"vincenty"}}{uses Vincenty formula for an ellipsoid
#' with equal major and minor axes}
#' }
#'
#' @returns numeric. Angle in radians
#'
#' @references
#' * Imboden, C. & Imboden, D. (1972). Formel fuer Orthodrome und Loxodrome bei
#' der Berechnung von Richtung und Distanz zwischen Beringungs- und
#' Wiederfundort.
#' *Die Vogelwarte* **26**, 336-346.
#' * Sinnott, Roger W. (1984). Virtues of the Haversine. *Sky and telescope*
#' **68**(2), 158.
#' Vincenty, T. (1975). Direct and inverse solutions of geodesics on the
#' ellipsoid with application of nested equations. *Survey Review*, **23**(176),
#' 88<U+2013>93. \doi{10.1179/sre.1975.23.176.88}.
#' * \url{http://www.movable-type.co.uk/scripts/latlong.html}
#' * \url{http://www.edwilliams.org/avform147.htm}
#'
#' @name spherical_angle
#'
#' @examples
#' berlin <- c(52.52, 13.41) |> deg2rad()
#' calgary <- c(51.04, -114.072) |> deg2rad()
#' orthodrome(berlin[1], berlin[2], calgary[1], calgary[2]) # 1.176406
#' haversine(berlin[1], berlin[2], calgary[1], calgary[2]) # 1.176406
#' vincenty(berlin[1], berlin[2], calgary[1], calgary[2]) # 1.176406
NULL

#' @rdname spherical_angle
#' @export
orthodrome <- function(lat1, lon1, lat2, lon2) {
  dlon <- lon2 - lon1

  sin_lat1 <- sin(lat1)
  sin_lat2 <- sin(lat2)
  cos_lat1 <- cos(lat1)
  cos_lat2 <- cos(lat2)
  cos_dlon <- cos(dlon)

  acos(pmin(pmax(sin_lat1 * sin_lat2 + cos_lat1 * cos_lat2 * cos_dlon, -1), 1))
}

#' @rdname spherical_angle
#' @export
haversine <- function(lat1, lon1, lat2, lon2) {
  dlon <- lon2 - lon1
  dlat <- lat2 - lat1

  hav_dlat <- hav(dlat)
  hav_dlon <- hav(dlon)
  hav_lat_sum <- hav(lat2 + lat1)

  arg <- hav_dlat + (1 - hav_dlat - hav_lat_sum) * hav_dlon

  ahav(arg)
}

#' @rdname spherical_angle
#' @export
vincenty <- function(lat1, lon1, lat2, lon2) {
  dlon <- lon2 - lon1

  sin_lat1 <- sin(lat1)
  cos_lat1 <- cos(lat1)
  sin_lat2 <- sin(lat2)
  cos_lat2 <- cos(lat2)
  sin_dlon <- sin(dlon)
  cos_dlon <- cos(dlon)

  y <- sqrt((cos_lat2 * sin_dlon)^2 +
              (cos_lat1 * sin_lat2 - sin_lat1 * cos_lat2 * cos_dlon)^2)
  x <- sin_lat1 * sin_lat2 + cos_lat1 * cos_lat2 * cos_dlon

  atan2(y, x)
}

#' @keywords internal
ddistance <- function(theta1, phi1, theta2, phi2, r = earth_radius()) {
  # Method after to Ziegler & Heidbach:  (2017)
  x1 <- r * cos(theta1) * cos(phi1)
  y1 <- r * cos(theta1) * sin(phi1)
  z1 <- r * sin(theta1)
  x2 <- r * cos(theta2) * cos(phi2)
  y2 <- r * cos(theta2) * sin(phi2)
  z2 <- r * sin(theta2)

  sqrt((x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2) # distance in km
}




#' Distance between points
#'
#' Returns the great circle distance between a location and all grid point in km
#'
#' @param lat1,lon1 numeric vector. coordinate of point(s) 1 (degrees).
#' @param lat2,lon2 numeric vector. coordinates of point(s) 2 (degrees).
#' @param r numeric. radius of the sphere (default = 6371.0087714 km, i.e. the
#' radius of the Earth)
#' @param method Character. Formula for calculating great circle distance,
#' one of:
#' \describe{
#' \item{\code{"haversine"}}{great circle distance based on the haversine
#' formula that is optimized for 64-bit floating-point numbers (the default)}
#' \item{\code{"orthodrome"}}{great circle distance based on the spherical law of cosines}
#' \item{\code{"vincenty"}}{distance based on the Vincenty formula for an
#' ellipsoid with equal major and minor axes}
#' \item{"euclidean"}{Euclidean distance (not great circle distance!)}
#' }
#'
#' @returns numeric vector with length equal to `length(lat1)` or `length(lat2)`
#'
#' @export
#'
#' @seealso [orthodrome()], [haversine()], [vincenty()]
#'
#' @examples
#' # Haversine: (4149.157, 2296.583) km
#' dist_greatcircle(lat1 = 20, lon1 = 12, lat2 = c(50, 30), lon2 = c(40, 32))
#'
#' # Orthodrome: (4149.157, 2296.583) km
#' dist_greatcircle(
#'   lat1 = 20, lon1 = 12, lat2 = c(50, 30), lon2 = c(40, 32),
#'   method = "orthodrome"
#' )
#'
#' # Vincenty: (4149.157, 2296.583) km
#' dist_greatcircle(
#'   lat1 = 20, lon1 = 12, lat2 = c(50, 30), lon2 = c(40, 32),
#'   method = "vincenty"
#' )
#'
#' # Euclidean (4076.220, 2284.169) km
#' dist_greatcircle(
#'   lat1 = 20, lon1 = 12, lat2 = c(50, 30), lon2 = c(40, 32),
#'   method = "euclidean"
#' )
dist_greatcircle <- function(lat1, lon1, lat2, lon2,
                             r = earth_radius(),
                             method = c("haversine", "orthodrome", "vincenty", "euclidean")) {
  method <- match.arg(method)
  #n <- length(lat1)
  if (!is.numeric(r) || length(r) != 1) stop("'r' must be a single numeric value.")
  # if (length(lon1) != n || length(lat2) != n || length(lon2) != n) {
  #   stop("lat1, lon1, lat2, and lon2 must all have the same length.")
  # }

  lat1_rad <- deg2rad(lat1)
  lon1_rad <- deg2rad(lon1)
  lat2_rad <- deg2rad(lat2)
  lon2_rad <- deg2rad(lon2)

  d <- switch(method,
              "haversine"  = haversine(lat1_rad, lon1_rad, lat2_rad, lon2_rad) * r,
              "orthodrome" = orthodrome(lat1_rad, lon1_rad, lat2_rad, lon2_rad) * r,
              "vincenty"   = vincenty(lat1_rad, lon1_rad, lat2_rad, lon2_rad) * r,
              "euclidean"  = ddistance(lat1_rad, lon1_rad, lat2_rad, lon2_rad, r)
  )
  return(d)
}

#' Shortest distance between pairs of geometries
#'
#' The shortest Great Circle distance between pairs of geometries
#'
#' @param x,line objects of class `sfg`, `sfc` or `sf`
#' @param ellipsoidal Logical. Whether the distance is calculated using
#' spherical distances ([sf::st_distance()]) or
#' ellipsoidal distances (`lwgeom::st_geod_distance()`).
#'
#' @returns numeric. Shortest distance in meters
#'
#' @export
#'
#' @examples
#' plate_boundary <- subset(plates, plates$pair == "na-pa")
#' shortest_distance_to_line(san_andreas, plate_boundary) |>
#'   head()
shortest_distance_to_line <- function(x, line, ellipsoidal = FALSE) {
  n <- sf::st_coordinates(x) |>
    nrow()

  suppressMessages(
    sf::sf_use_s2(!ellipsoidal)
  )

  line_pts <- sf::st_cast(line, "MULTILINESTRING", warn = FALSE) |>
    sf::st_cast("LINESTRING", warn = FALSE) |>
    sf::st_cast("POINT", warn = FALSE)

  dmat <- sf::st_distance(
    x,
    line_pts,
    which = "Great Circle",
    radius = 1000 * earth_radius()
  )
  suppressMessages(
    sf::sf_use_s2(ellipsoidal)
  )

  if (n == 1) {
    min(dmat)
  } else {
    sapply(1:n, function(i) min((dmat[i, ])))
  }
}


#' @title Azimuth Between two Points
#'
#' @description Calculate initial bearing (or forward azimuth/direction) to go
#' from point `a` to point `b` following great circle arc on a
#' sphere.
#'
#' @param lat_a,lat_b Numeric. Latitudes of a and b (in degrees).
#' @param lon_a,lon_b Numeric. Longitudes of a and b (in degrees).
#' @details [get_azimuth()] is based on the spherical law of tangents.
#' This formula is for the initial bearing (sometimes referred to as
#' forward azimuth) which if followed in a straight line along a great circle
#' arc will lead from the start point `a` to the end point `b`.
#' \deqn{\theta = \arctan2 (\sin \Delta\lambda
#' \cos\psi_2, \cos\psi_1 \sin\psi_1-\sin\psi_1 \cos\psi_2 \cos\Delta\lambda)}
#' where  \eqn{\psi_1, \lambda_1} is the start point, \eqn{\psi_2},
#' \eqn{\lambda_2} the end point (\eqn{\Delta\lambda} is the difference in
#' longitude).
#'
#' @returns numeric. Azimuth in degrees
#'
#' @references \url{http://www.movable-type.co.uk/scripts/latlong.html}
#'
#' @export
#'
#' @examples
#' berlin <- c(52.517, 13.4) # Berlin
#' tokyo <- c(35.7, 139.767) # Tokyo
#' get_azimuth(berlin[1], berlin[2], tokyo[1], tokyo[2]) # 41.57361
get_azimuth <- function(lat_a, lon_a, lat_b, lon_b) {
  dr <- pi / 180 # degrees to radians

  la <- lat_a * dr
  lb <- lat_b * dr

  dphi <- (lon_b - lon_a) * dr
  cos_lb <- cos(lb)

  y <- sin(dphi) * cos_lb
  x <- cos(la) * sin(lb) - sin(la) * cos_lb * cos(dphi)
  theta <- atan2(y, x)

  (theta / dr) %% 360 # Convert to degrees and wrap to [0, 360)
}
