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


#' Azimuth relative to PoR
#'
#' Transform azimuth from geographical to PoR system according to Wdowinski 1998
#'
#' @param x \code{"data.frame"}. coordinates of data point (lat, lon) and the azimuth
#' @param euler \code{"data.frame"}. coordinates of Euler pole (PoR)
#' @return \code{"numeric"}. Transformed azimuth relative to PoR
#' @references Wdowinski, S., 1998, A theory of intraplate
#'   tectonics: Journal of Geophysical Research: Solid Earth, v. 103, p.
#'   5037-5059, \doi{10.1029/97JB03390}.
#' @examples
#' \dontrun{
#' data("nuvel1")
#' # North America relative to Pacific plate:
#' euler <- subset(nuvel1, nuvel1$plate.rot == "na")
#'
#' data("san_andreas")
#' res <- PoR_azimuth(san_andreas, euler)
#' head(res)
#' }
PoR_azimuth <- function(x, euler) {
  # .Deprecated("PoR_shmax")
  # convert latitude to colatitude and to radians
  x <- data.frame(lat = x$lat, lon = x$lon, azi = x$azi) * pi / 180
  ep <- data.frame(lat = euler$lat, lon = euler$lon) * pi / 180

  S <-
    cos(ep$lat) * cos(ep$lon) * cos(x$lat) * cos(x$lon) +
    cos(ep$lat) * sin(ep$lon) * cos(x$lat) * sin(x$lon) +
    sin(ep$lat) * sin(x$lat)

  beta <- acos(
    (sin(ep$lat) - S * sin(x$lat)) /
      (sqrt(1 - S * S) * abs(cos(x$lat)))
  )

  azi.por <- x$azi - beta
  (azi.por * 180 / pi + 360) %% 180
}


#' Displacement vector and stress matrix in PoR
#'
#' @param x \code{sf} object of the data points in the geographical coordinate
#' system
#' @param euler \code{data.frame} of the geographical coordinates of the Euler pole
#' (\code{lat}, \code{lon})
#' @param positive Numeric. Equatorial relative plate boundary motion
#' (displacement).
#' @param positive Sign of the equatorial relative plate boundary motion
#' (displacement).
#' If \code{tangential == TRUE}, positive values indicate right-lateral plate
#' boundary displacement and negative values represent left-lateral displacement.
#' If \code{tangential == FALSE}, positive values indicate outward-moving plate
#' boundary offset and negative values represent inward-moving plate boundary
#' displacement.
#' @param v Numeric. Poisson's ratio. default is 0.25
#' @param E Numeric. Young's modulus. default is 50 (GPa)
#' @param tangential Logical Whether the plate boundary is a tangential
#' boundary (\code{TRUE}) or an inward and outward boundary (\code{FALSE}, the
#' default).
#' @returns \code{sf} object of the data points in the PoR coordinate system
#' and the components of the displacement vector or the stress matrix.
#' @references Wdowinski, S., 1998, A theory of intraplate
#'   tectonics. *Journal of Geophysical Research: Solid Earth*, **103**,
#'   5037-5059, doi: 10.1029/97JB03390.
#' @importFrom magrittr %>%
#' @importFrom sf st_coordinates st_as_sf
#' @name stressstrain
#' @examples
#' \dontrun{
#' data("nuvel1")
#' na_pa <- subset(nuvel1, nuvel1$plate.rot == "na")
#'
#' data("plates")
#' plate_boundary <- subset(plates, plates$pair == "na-pa")
#'
#' data("san_andreas")
#' res <- displacement_vector(x = san_andreas, ep = na_pa, tangential = TRUE, positive = FALSE)
#' head(res)
#'
#' res2 <- stress_matrix(x = san_andreas, euler = na_pa, tangential = TRUE, positive = FALSE)
#' head(res2)
#' }
NULL

#' @rdname stressstrain
displacement_vector <- function(x, euler, tangential = FALSE, positive = TRUE) {
  if (positive) {
    u <- abs(euler$angle)
  } else {
    u <- -abs(euler$angle)
  }
  x.por <- geographical_to_PoR(x, euler) %>%
    sf::st_coordinates()

  lon1 <- min(x.por[, 1])
  lon2 <- max(x.por[, 1])

  lat1 <- min(x.por[, 2])
  lat2 <- max(x.por[, 2])

  u_lat <- u_lon <- c()
  for (i in seq_along(x.por[, 1])) {
    if (!tangential) {
      u_lat[i] <- 0
      u_lon[i] <- u * sind(90 - x.por[i, 2])^2 * (lon1 - x.por[i, 1]) / (lon2 - lon1)
    }

    if (tangential) {
      u_lat[i] <- 0
      u_lon[i] <- u * sind(90 - x.por[i, 2])^2 * (x.por[i, 2] - lat1) / (lat1 - lat2)
    }
  }

  sf::st_as_sf(
    x = data.frame(x.por, d_sc = u_lon, d_gc = u_lat),
    coords = c(1, 2),
    crs = PoR_crs(euler)
  )
}

#' @rdname stressstrain
stress_matrix <- function(x, euler, tangential = FALSE, positive = FALSE, v = .25, E = 50) {
  if (positive) {
    u <- abs(euler$angle)
  } else {
    u <- -abs(euler$angle)
  }
  x.por <- geographical_to_PoR(x, euler) %>%
    sf::st_coordinates()

  lon1 <- min(x.por[, 1])
  lon2 <- max(x.por[, 1])

  lat1 <- min(x.por[, 2])
  lat2 <- max(x.por[, 2])

  s_xx <- s_xz <- s_zx <- s_zz <- c()
  if (!tangential) {
    d <- lat2 - lat1
    A <- matrix(data = NA, nrow = 2, ncol = 2)
    A[1, 1] <- v / (1 - v)
    A[1, 2] <- 0
    A[2, 2] <- 1
    A[2, 1] <- 0

    Q <- E / (1 + v) * A

    for (i in seq_along(x.por[, 1])) {
      S <- -(u * sind(90 - x.por[i, 2])^2 / d) * Q

      s_xx[i] <- S[1, 1]
      s_xz[i] <- S[1, 2]
      s_zx[i] <- S[2, 1]
      s_zz[i] <- S[2, 2]
    }
  }

  if (tangential) {
    d <- lon2 - lon1
    A <- matrix(data = NA, nrow = 2, ncol = 2)
    A[1, 1] <- 0
    A[1, 2] <- 1
    A[2, 2] <- 0
    A[2, 1] <- 1

    Q <- E / (1 + v) * A

    for (i in seq_along(x.por[, 1])) {
      S <- -(u * sind(90 - x.por[i, 2])^2 / 2 * d) * Q

      s_xx[i] <- S[1, 1]
      s_xz[i] <- S[1, 2]
      s_zx[i] <- S[2, 1]
      s_zz[i] <- S[2, 2]
    }
  }

  sf::st_as_sf(
    x = data.frame(x.por, s_xx, s_xz, s_zx, s_zz),
    coords = c(1, 2),
    crs = PoR_crs(euler)
  )
}

# Coordinate transformation using rotation matrices ####

#' Conversion between spherical PoR to geographical coordinate system
#'
#' Transformation from spherical PoR to geographical coordinate system and
#' vice versa
#'
#' @param x \code{sf} or \code{data.frame} containing (co)lat and lon coordinates
#' (\code{lat}, \code{lon}) of the points to be transformed
#' @param euler \code{data.frame} of the geographical coordinates of the Euler pole
#' (\code{(co)lat}, \code{lon})
#' @return \code{data.frame} with the PoR coordinates
#' (\code{colat.PoR}, \code{lon.PoR} or \code{lat.PoR}, \code{lon.PoR})
#' @param spherical logical. Whether x or the return are in spherical
#' coordinates, i.e. colat-lon (TRUE, the default), or in lat-lon (FALSE)
#' @name por_conversion_df
#' @references Wdowinski, S., 1998, A theory of intraplate
#'   tectonics: Journal of Geophysical Research: Solid Earth, v. 103, p.
#'   5037-5059, \doi{10.1029/97JB03390}.
#' @examples
#' data("nuvel1")
#' euler <- subset(nuvel1, nuvel1$plate.rot == "na") # North America relative to Pacific plate
#' data("san_andreas")
#' san_andreas.por <- geographical_to_PoR2(san_andreas, euler)
#' head(san_andreas.por)
#' head(PoR_to_geographical2(san_andreas.por, euler))
NULL

#' @rdname por_conversion_df
#' @export
PoR_to_geographical2 <- function(x, euler, spherical = TRUE) {
  ep.geo <- c(euler$lat, euler$lon)
  lat <- lon <- c()
  for (i in seq_along(x$lon)) {
    if (spherical) {
      x_geo.i <- PoR_to_geographical_vec(c(x$colat.PoR[i], x$lon.PoR[i]), euler = ep.geo, spherical = TRUE)
    } else {
      x_geo.i <- PoR_to_geographical_vec(c(x$lat.PoR[i], x$lon.PoR[i]), euler = ep.geo, spherical = FALSE)
    }
    lat[i] <- x_geo.i[1]
    lon[i] <- x_geo.i[2]
  }
  data.frame(lat = lat, lon = lon)
}

#' @rdname por_conversion_df
#' @export
geographical_to_PoR2 <- function(x, euler, spherical = TRUE) {
  ep.geo <- c(euler$lat, euler$lon)
  colat.PoR <- lon.PoR <- c()
  for (i in seq_along(x$lon)) {
    x_por.i <- geographical_to_PoR_vec(c(x$lat[i], x$lon[i]), euler = ep.geo, spherical)
    colat.PoR[i] <- x_por.i[1]
    lon.PoR[i] <- x_por.i[2]
  }
  if (spherical) {
    data.frame(colat.PoR = colat.PoR, lon.PoR = lon.PoR)
  } else {
    data.frame(lat.PoR = colat.PoR, lon.PoR = lon.PoR)
  }
}

#' Conversion between PoR to geographical coordinate system
#'
#' Helper function for the transformation from PoR to geographical coordinate
#' system or vice versa
#'
#' @param x,euler two-column vectors containing the (co)lat and lon coordinates
#' @param spherical logical. Whether x or the return are in spherical
#' coordinates
#' @references Wdowinski, S., 1998, A theory of intraplate
#'   tectonics: Journal of Geophysical Research: Solid Earth, v. 103, p.
#'   5037-5059, \doi{10.1029/97JB03390}.
#' @name por_conversion_vec
#' @examples
#' \dontrun{
#' ep.geo <- c(20, 33)
#' q.geo <- c(10, 45)
#' q.por <- geographical_to_PoR_vec(q.geo, ep.geo)
#' q.por
#' PoR_to_geographical_vec(q.por, ep.geo)
#' }
NULL

#' @rdname por_conversion_vec
PoR_to_geographical_vec <- function(x, euler, spherical = TRUE) {
  x[2] <- longitude_modulo(x[2] + 180)

  if (spherical) {
    x.por.cart <- spherical_to_cartesian(x)
  } else {
    x.por.cart <- geographical_to_cartesian(x)
  }
  x.cart <- t(rotmat_z(180 - euler[2])) %*% t(rotmat_y(90 - euler[1])) %*% x.por.cart
  cartesian_to_geographical(x.cart)
}

#' @rdname por_conversion_vec
geographical_to_PoR_vec <- function(x, euler, spherical = TRUE) {
  x.cart <- geographical_to_cartesian(x)
  x.cart.por <- rotmat_y(90 - euler[1]) %*% rotmat_z(180 - euler[2]) %*% x.cart

  if (spherical) {
    res <- cartesian_to_spherical(x.cart.por)
  } else {
    res <- cartesian_to_geographical(x.cart.por)
  }
  res[2] <- longitude_modulo(res[2] - 180)
  return(res)
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

# Quaternions ####

# Q4 <- function(q) {
#   stopifnot(is.Q4(q))
#   q$Sc + q$Vec
# }

rotmat_to_QScVec <- function(x, normalize = FALSE) {
  stopifnot(is.matrix(x), is.logical(normalize))
  q0 <- sqrt(abs((1 + x[1, 1] + x[2, 2] + x[3, 3]) / 4))
  q1 <- sqrt(abs((1 + x[1, 1] - x[2, 2] - x[3, 3]) / 4))
  q2 <- sqrt(abs((1 - x[1, 1] + x[2, 2] - x[3, 3]) / 4))
  q3 <- sqrt(abs((1 - x[1, 1] - x[2, 2] + x[3, 3]) / 4))

  maxq <- max(c(q0, q1, q2, q3))

  if (maxq == q0) {
    q1 <- (x[3, 2] - x[2, 3]) / 4 * q0
    q2 <- (x[1, 3] - x[3, 1]) / 4 * q0
    q3 <- (x[2, 1] - x[1, 2]) / 4 * q0
  } else if (maxq == q1) {
    q0 <- (x[3, 2] - x[2, 3]) / 4 * q1
    q2 <- (x[1, 2] + x[2, 1]) / 4 * q1
    q3 <- (x[1, 3] + x[3, 1]) / 4 * q1
  } else if (maxq == q2) {
    q0 <- (x[1, 3] - x[3, 1]) / 4 * q2
    q1 <- (x[1, 2] + x[2, 1]) / 4 * q2
    q3 <- (x[2, 3] + x[3, 2]) / 4 * q2
  } else if (maxq == q3) {
    q0 <- (x[2, 1] - x[1, 2]) / 4 * q3
    q1 <- (x[1, 3] + x[3, 1]) / 4 * q3
    q2 <- (x[2, 3] + x[3, 2]) / 4 * q3
  }
  Q4_to_QScVec(c(q0, q1, q2, q3), normalize)
}

QScVec_to_rotmat <- function(x, normalize = TRUE) {
  stopifnot(is.Q4(x), is.logical(normalize))
  q <- QScVec_to_Q4(x)

  m <- rbind(
    c(q[1]^2 + q[2]^2 - q[3]^2 - q[4]^2, 2 * q[2] * q[3] - 2 * q[1] * q[4], 2 * q[2] * q[4] + 2 * q[1] * q[3]),
    c(2 * q[2] * q[3] + 2 * q[1] * q[4], q[1]^2 - q[2]^2 + q[3]^2 - q[4]^2, 2 * q[3] * q[4] - 2 * q[1] * q[2]),
    c(2 * q[2] * q[4] - 2 * q[1] * q[3], 2 * q[3] * q[4] + 2 * q[1] * q[2], q[1]^2 - q[2]^2 - q[3]^2 + q[4]^2)
  )
  if (normalize) {
    m <- normalize_matrix(m)
  }
  return(m)
}

normalize_matrix <- function(x) {
  x / max(x)
}


# microbenchmark::microbenchmark(
#   geographical_to_PoR_quat(q.geo, ep.geo),
#   geographical_to_PoR_vec(q.geo, ep.geo, spherical = FALSE),
#   geographical_to_PoR_sf(data.frame(lat = q.geo[1], lon = q.geo[2]) %>% sf::st_as_sf(coords=c('lon', 'lat')), euler = data.frame(lat = ep.geo[1], lon = ep.geo[2]))
# )





#' @title Median and statistics on Pi-periodic Data
#'
#' @description Calculate the mean, median, quartile, interquartile range,
#' variance, deviation, and error of directional data.
#'
#' @param x Numeric vector in degrees.
#' @param na.rm logical. Should missing values (including `NaN`) be removed?
#'
#' @return Numeric vector
#'
#' @details Quasi median on the circle, quasi quartiles on a circle, quasi
#' interquartile range on a circle.
#'
#' @source [median()], [quantile()], and [IQR()] are the
#' equivalents for non-periodic data.
#'
#' @references
#' * Ratanaruamkarn, S., Niewiadomska-Bugaj, M., Wang, J.-C. (2009).
#' A New Estimator of a Circular Median. *Communications in Statistics -
#' Simulation and Computation*, **38**(6), 1269-1291.
#' \doi{10.1080/03610910902899950}.
#' * Reiter, K., Heidbach, O., Schmitt, D., Haug, K., Ziegler, M., & Moeck, I.
#' (2014). A revised crustal stress orientation database for Canada.
#' *Tectonophysics*, **636**, 111-124. \doi{10.1016/j.tecto.2014.08.006}
#'
#' @importFrom stats median
#' @examples
#' x <- c(0, 45, 55, 40 + 180, 50 + 180, NA)
#' circular_mean(x)
#' circular_quasi_median(x)
#' circular_quasi_quantile(x)
#' circular_quasi_IQR(x)
#' circular_var(x)
#' circular_mean_deviation(x)
#' circular_median_deviation(x)
#' circular_mean_error(x)
#'
#' data("san_andreas")
#' circular_quasi_median(san_andreas$azi)
#' @name circle_median
NULL

#' @rdname circle_median
#' @export
circular_quasi_median <- function(x, na.rm = TRUE) {
  stopifnot(any(is.numeric(x)), is.logical(na.rm))

  if (na.rm) {
    x <- as.numeric(na.omit(x))
  }
  x <- deg2rad(x)
  x <- sort(x[!is.na(x)])
  n <- length(x)

  if (n %% 2 != 0) { # if odd
    m <- (n - 1) / 2
    ss <- sin(x[m + 1])
    cs <- cos(x[m + 1])
  } else { # if even
    m <- n / 2
    ss <- sin(x[m]) + sin(x[m + 1])
    cs <- cos(x[m]) + cos(x[m + 1])
  }
  atan2d_spec(ss, cs) %% 180
}

#' @rdname circle_median
#' @export
circular_mean <- function(x, na.rm = TRUE) {
  circular_weighted_mean(x, w = NULL, na.rm)
}

#' @rdname circle_median
#' @export
circular_quasi_quantile <- function(x, na.rm = TRUE) {
  stopifnot(any(is.numeric(x)), is.logical(na.rm))

  if (na.rm) {
    x <- as.numeric(na.omit(x))
  }
  x <- x
  x <- sort(x[!is.na(x)])
  n <- length(x)

  if (n > 3) {
    med <- circular_quasi_median(x)

    if (n %% 4 == 0) {
      m <- n / 4
      lq <- atand(
        sind(x[m + 1]) / cosd(x[m + 1])
      )
      uq <- atand(
        sind(x[3 * m + 1]) / cosd(x[3 * m + 1])
      )
    } else if (n %% 4 == 1) {
      m <- (n - 1) / 4
      lq <- atand(
        (3 * sind(x[m]) + sind(x[m + 1])) /
          (3 * cosd(x[m]) + cosd(x[m + 1]))
      )
      uq <- atand(
        (3 * sind(x[3 * m]) + sind(x[3 * m + 1])) /
          (3 * cosd(x[3 * m]) + cosd(x[3 * m + 1]))
      )
    } else if (n %% 4 == 2) {
      m <- (n - 2) / 4
      lq <- atand((sind(x[m]) + sind(x[m + 1])) /
        (cosd(x[m]) + cosd(x[m + 1])))
      uq <- atand((sind(x[3 * m]) + sind(x[3 * m + 1])) /
        (cosd(x[3 * m]) + cosd(x[3 * m + 1])))
    } else { # if (n %% 4 == 3) {
      m <- (n - 2) / 4
      lq <- atand((sind(x[m]) + 3 * sind(x[m + 1])) /
        (cosd(x[m]) + 3 * cosd(x[m + 1])))
      uq <- atand((sind(x[3 * m]) +
        3 * sind(x[3 * m + 1])) /
        (cosd(x[3 * m]) +
          3 * cosd(x[3 * m + 1])))
    }

    quantiles <- c(x[1], lq, med, uq, x[length(x)])
    names(quantiles) <- c("0%", "25%", "50%", "75%", "100%")
    return(as.numeric(quantiles))
  } else {
    message("x needs more than 3 values")
    return(NULL)
  }
}

#' @rdname circle_median
#' @export
circular_quasi_IQR <- function(x, na.rm = TRUE) {
  quantiles <- circular_quasi_quantile(x, na.rm)
  deviation_norm(quantiles[4] - quantiles[2])
}

#' @rdname circle_median
#' @export
circular_var <- function(x, na.rm = TRUE, Batschelet = FALSE) {
  weighted_circular_variance(x, w = NULL, na.rm, Batschelet)
}

#' @rdname circle_median
#' @export
circular_standard_deviation <- function(x, na.rm = TRUE) {
  circular_weighted_sd(x, w = NULL, na.rm = TRUE)
}

#' @rdname circle_median
#' @export
circular_mean_deviation <- function(x, na.rm = TRUE) {
  stopifnot(any(is.numeric(x)), is.logical(na.rm))

  if (na.rm) {
    x <- as.numeric(na.omit(x))
  }
  x <- x %% 180
  n <- length(x)

  for (i in 1:n) {
    k <- abs(
      180 - abs(x[i] - circular_quasi_median(x))
    )
  }
  180 - (1 / n * sum(k))
}

#' @rdname circle_median
#' @export
circular_median_deviation <- function(x, na.rm = TRUE) {
  stopifnot(any(is.numeric(x)), is.logical(na.rm))

  if (na.rm) {
    x <- as.numeric(na.omit(x))
  }
  x <- x %% 180

  for (i in seq_along(x)) {
    k <- 180 - abs(180 - abs(x[i] - circular_quasi_median(x)))
  }
  stats::median(k)
}

#' @rdname circle_median
#' @export
circular_mean_error <- function(x, na.rm = TRUE) {
  stopifnot(any(is.numeric(x)), is.logical(na.rm))

  if (na.rm) {
    x <- as.numeric(na.omit(x))
  }

  x <- x %% 180
  n <- length(x)

  for (i in 1:n) {
    k <- abs(180 - abs(x[i] - circular_quasi_median(x)))
  }
  180 - (1 / n * sum(k))
}


#' Circular distance between two angles
#'
#' @param x,y vector of numeric values in degrees
#' @param na.rm logical value indicating whether \code{NA} values in \code{x}
#' should be stripped before the computation proceeds.
#' @param axial logical. Whether the data are axial, i.e. pi-periodical
#' (TRUE, the default) or circular, i.e. 2pi-periodical (FALSE).
#' @examples
#' \dontrun{
#' a <- c(0, 2, 359, 6, 354)
#' b <- a + 20
#' circular_distance(a, b, axial)
#' }
circular_distance <- function(x, y, axial = TRUE, na.rm = TRUE) {
  if (axial) {
    f <- 2
  } else {
    f <- 1
  }

  if (na.rm) {
    x <- na.omit(x)
  }

  cdist <- 1 - cosd(f * (x - y))
  cdist # / f
}

#' Circular dispersion
#'
#' Dispersion of angles about a given angle.
#'
#' @param x Data values. A vector of numeric values in degrees
#' @param from numeric. the angle (in degree) for which the dispersion should be measured.
#' If NULL, the circular median is used.
#' @param w (optional) Weights. A vector of positive numbers, of the same length as
#' \code{x}.
#' @param na.rm logical value indicating whether \code{NA} values in \code{x}
#' should be stripped before the computation proceeds.
#' @param axial logical. Whether the data are axial, i.e. pi-periodical
#' (TRUE, the default) or circular, i.e. 2pi-periodical (FALSE).
#' @importFrom tidyr drop_na
#' @examples
#' \dontrun{
#' data("san_andreas")
#' circular_dispersion(san_andreas$azi)
#'
#' data("nuvel1")
#' ep <- subset(nuvel1, nuvel1$plate.rot == "na")
#' sa.por <- PoR_shmax(san_andreas, ep, "right")
#' circular_dispersion(sa.por$azi.PoR, from = 135)
#' }
circular_dispersion <- function(x, from = NULL, w = NULL, axial = TRUE, na.rm = TRUE) {
  stopifnot(any(is.numeric(x)), is.logical(na.rm), is.logical(axial))

  if (is.null(from)) {
    from <- circular_mean(x, w, axial, na.rm)
  }

  if (axial) {
    f <- 2
  } else {
    f <- 1
  }

  if (is.null(w)) {
    w <- rep(1, times = length(x))
  } else {
    w <- as.numeric(w)
  }

  data <- data.frame(x, w)
  if (na.rm) {
    data <- tidyr::drop_na(data)
  }

  w <- data$w
  Z <- sum(w)
  cosf <- w * cosd(f * (data$x - from))
  sum(1 - cosf) / Z
}

sample_circular_dispersion <- function(x, axial = TRUE) {
  if (axial) {
    f <- 2
  } else {
    f <- 1
  }

  R <- mean_resultant(deg2rad(f * x))
  R <- mean_resultant(deg2rad(2 * f * x))
  (1 - (R2 / 2)) / (2 * R^2)
}



# Distribution ####
pvm <- function(theta, mu, kappa, acc = 1e-20) {
  theta <- theta %% (2 * pi)
  mu <- mu %% (2 * pi)
  pvm.mu0 <- function(theta, kappa, acc) {
    flag <- "true"
    p <- 1
    sum <- 0
    while (flag == "true") {
      term <- (besselI(x = kappa, nu = p, expon.scaled = FALSE) *
        sin(p * theta)) / p
      sum <- sum + term
      p <- p + 1
      if (abs(term) < acc) {
        flag <- "false"
      }
    }
    theta / (2 * pi) + sum / (pi * besselI(
      x = kappa, nu = 0,
      expon.scaled = FALSE
    ))
  }
  if (mu == 0) {
    result <- pvm.mu0(theta, kappa, acc)
  } else {
    if (theta <= mu) {
      upper <- (theta - mu) %% (2 * pi)
      if (upper == 0) {
        upper <- 2 * pi
      }
      lower <- (-mu) %% (2 * pi)
      result <- pvm.mu0(upper, kappa, acc) - pvm.mu0(
        lower,
        kappa, acc
      )
    } else {
      upper <- theta - mu
      lower <- mu %% (2 * pi)
      result <- pvm.mu0(upper, kappa, acc) + pvm.mu0(
        lower,
        kappa, acc
      )
    }
  }
  result
}


#' Kuiper test of uniformity
#'
#' Hypothesis tests of uniformity for circular data.
#'
#' @param x numeric vector containing the circular data which are expressed in degrees
#' @details Kuiper test of uniformity,
#' @returns the test statistic
#' @examples
#' data(san_andreas)
#' kuiper_statistics(san_andreas$azi)
kuiper_statistics <- function(x) {
  x <- na.omit(x)
  u <- sort(deg2rad(x) %% (2 * pi)) / (2 * pi)
  n <- length(u)
  i <- 1:n
  D.P <- max(i / n - u)
  D.M <- max(u - (i - 1) / n)
  sqrt_n <- sqrt(n)
  V <- D.P + D.M
  V * (sqrt_n + 0.155 + 0.24 / sqrt_n)
}

#' Watson's goodness of fit test for the von Mises or circular uniform distribution, Rayleigh Test of Uniformity
#'
#' Watson's goodness of fit test for the von Mises or circular uniform distribution.
#'
#' @param x numeric vector containing the circular data which are expressed in degrees
#' @param dist Distribution to test for. The default is the uniform distribution (\code{"uniform"}), and \code{"vm"} to test for the von Mises distribution.
#' @param ... optional parameter for [circular_mean()]
#' @returns the test statistic
#' @examples
#' data(san_andreas)
#' watson_statistics(san_andreas$azi)
watson_statistics <- function(x, dist = c("uniform", "vm"), ...) {
  dist <- match.arg(dist)
  x <- na.omit(x)
  n <- length(x)

  if (dist == "uniform") {
    u <- sort(deg2rad(x)) / (2 * pi)
    u.bar <- mean(u)
    i <- 1:n
    sum.terms <- (u - u.bar - (2 * i - 1) / (2 * n) + 0.5)^2
    u2 <- sum(sum.terms) + 1 / (12 * n)
    (u2 - 0.1 / n + 0.1 / (n^2)) * (1 + 0.8 / n)
  } else if (dist == "vm") {
    mu.hat <- tectonicr::circular_mean(x, ...) %>% deg2rad()
    kappa.hat <- est.kappa(deg2rad(x))
    x <- (deg2rad(x) - mu.hat) %% (2 * pi)
    x <- matrix(x, ncol = 1)
    z <- apply(x, 1, pvm, 0, kappa.hat)
    z <- sort(z)
    z.bar <- mean(z)
    i <- 1:n
    sum.terms <- (z - (2 * i - 1) / (2 * n))^2
    sum(sum.terms) - n * (z.bar - 0.5)^2 + 1 / (12 * n)
  }
}
