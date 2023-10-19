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
  x.por <- geographical_to_PoR(x, euler) |>
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
  x.por <- geographical_to_PoR(x, euler) |>
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
#   geographical_to_PoR_sf(data.frame(lat = q.geo[1], lon = q.geo[2]) |> sf::st_as_sf(coords=c('lon', 'lat')), euler = data.frame(lat = ep.geo[1], lon = ep.geo[2]))
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
#' circular_mean_deviation(x, 50)
#' circular_median_deviation(x)
#' circular_mean_error(x)
#'
#' data("san_andreas")
#' circular_quasi_median(san_andreas$azi)
#' @name circle_median
NULL

#' @rdname circle_median
#' @export
circular_mean_deviation <- function(x, y, axial = TRUE, na.rm = TRUE) {
  if (axial) {
    f <- 2
    mod <- 180
  } else {
    f <- 1
    mod <- 360
  }
  x <- (x * f) %% (2 * pi)
  y <- (y * f) %% (2 * pi)


  if (na.rm) {
    x <- as.numeric(na.omit(x))
  }

  k <- abs(
    180 - abs(x - y)
  )
  cmd <- 180 - ((1 / n) * sum(k))
  (cmd / f) %% mod
}

#' @rdname circle_median
#' @export
circular_median_deviation <- function(x, na.rm = TRUE) {
  stopifnot(any(is.numeric(x)), is.logical(na.rm))

  if (na.rm) {
    x <- as.numeric(na.omit(x))
  }
  x <- x %% 180

  cmed <- 180 - abs(180 - abs(x - circular_median(x)))
  stats::median(circular_mean_deviation(x, circular_median(x)))
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





#' Critical value for Rayleigh test
#'
#' Estimate critical values of the Rayleigh test for testing whether the population of
#' circular data from which a sample is drawn differs from randomness
#'
#' @param n integer. the number of data
#' @param alpha numeric. the probability level
#' @references Wilkie (1983): Rayleigh Test for Randomness of Circular Data. Appl. Statist. 32, No. 3, pp. 311-312
#' @returns numeric.
#' @export
#' @examples
#' estimate_critical(n = 500, alpha = 0.050)
estimate_critical <- function(n, alpha = 0.05) {
  -1 * log(alpha) - (2 * log(alpha) + log(alpha)^2) / (4 * n)
}



#' Extract azimuths of line segments
#'
#' @param x sf object of type `"LINESTRING"` or `"MULTILINESTRING"`
#'
#' @return sf object of type `"POINT"` with the columns and entries of the first row of `x`
#'
#' @details
#' It is recommended to perform `line_azimuth()` on single line objects, i.e.
#' type `"LINESTRING"`, instead of `"MULTILINESTRING"`. This is because the azimuth
#' of the last point of a line will be calculated to the first point of the
#' next line otherwise. This will cause a warning message. For `MULTILINESTRING`
#' objects, use `lines_azimuths()`.
#'
#' @importFrom sf st_cast st_coordinates st_as_sf st_crs st_drop_geometry
#' @export
#'
#' @name line_azimuth
#' @examples
#' data("plates", package = "tectonicr")
#' subset(plates, pair == "af-eu") |>
#'   smoothr::densify() |>
#'   line_azimuth()
#'
#' line_azimuth(plates)
NULL

#' @rdname line_azimuth
#' @export
line_azimuth <- function(x) {
  if (nrow(x) > 1) warning("It is recommended to only use single line objects")
  if (any(sf::st_geometry_type(x) == "MULTILINESTRING")) warning("It is recommended to only use single line objects")
  mat <- x |>
    sf::st_cast("POINT") |>
    sf::st_coordinates()

  n <- nrow(mat)

  a <- numeric()
  for (i in 1:(n - 1)) {
    a[i] <- get_azimuth(mat[i, 2], mat[i, 1], mat[i + 1, 2], mat[i + 1, 1])
  }
  data.frame(
    x = mat[1:(n - 1), 1],
    y = mat[1:(n - 1), 2],
    azi = a
  ) |>
    cbind(sf::st_drop_geometry(x[1, ])) |>
    sf::st_as_sf(coords = c("x", "y"), crs = sf::st_crs(x))
}

#' @rdname line_azimuth
#' @export
lines_azimuths <- function(x) {
  for (i in 1:nrow(x)) {
    ai <- line_azimuth(x[i, ])
    if (i == 1) {
      a <- ai
    } else {
      a <- rbind(a, ai)
    }
  }
  return(a)
}



#' SHmax direction resulting from multiple plate boundaries
#'
#' Calculates a \eqn{\sigma_{Hmax}}{SHmax} direction at given coordinates,
#' sourced by multiple plate boundaries. This first-order approximation is the
#' circular mean of the superimposed theoretical directions, weighted by the
#' rotation rates of the underlying PoRs.
#'
#' @param df \code{data.frame} containing the coordinates of the point(s)
#' (\code{lat}, \code{lon}), and the direction of
#' \eqn{\sigma_{Hmax}}{SHmax} \code{azi} (in degrees)
#' @param PoRs data.frame or object `"euler.pole"` that must contain `lat`,
#' `lon` and `angle`
#' @param types character vector with length equal to number of rows in `PoRs`.
#' Type of plate boundary. Can be \code{"out"}, \code{"in"}, \code{"right"}, or
#' \code{"left"} for outward, inward, right-lateral, or left-lateral
#' moving plate boundaries, respectively.
#' @param absolute logical. Whether the resultant azimuth should be weighted
#' using the absolute rotation at the points or the angular rotation of the PoRs.
#' @param PoR_weighting (optional) numeric vector with length equal to number of rows in
#' `PoRs`. Extra weightings for the used PoRs.
#'
#' @seealso [model_shmax()]
#'
#' @return numeric. Resultant azimuth in degrees and geographical CRS
#' @export
#'
#' @examples
#' data(san_andreas)
#' data(nuvel1)
#' pors <- subset(nuvel1, plate.rot %in% c("eu", "na"))
#' superimposed_shmax(san_andreas, pors, types = c("in", "right"), PoR_weighting = c(2, 1))
superimposed_shmax <- function(df, PoRs, types = c("in", "out", "right", "left"), absolute = TRUE, PoR_weighting = NULL) {
  res <- c()
  lats <- c()
  if(is.null(PoR_weighting)){
    PoR_weighting <- rep(1, nrow(PoRs))
  }

  if(!absolute){
    lat_j <- rep(1, nrow(df))
    col <- 1
    while(col <= nrow(PoRs)){
      lats <- cbind(lats, lat_j)
      col <- col+1
    }
  }

  for(i in seq_along(PoRs$lat)){
    res_i <- model_shmax(df, PoRs[i, ])
    if(absolute) lat_i <- PoR_coordinates(df, PoRs[i, ])$lat.PoR
    if(types[i] == "in"){
      azi <- res_i$sc
    } else if(types[i] == "out") {
      azi <- res_i$gc
    } else if (types[i] == "right") {
      azi <- res_i$ld.cw
    } else {
      azi <- res_i$ld.ccw
    }
    res <- cbind(res, azi)
    if(absolute) lats <- cbind(lats, lat_i)
  }

  rot <- PoR_weighting * PoRs$angle * cosd(lats)
  azi_res <- numeric()
  for(j in seq_along(res[, 1])){
      azi_res[j] <- circular_mean(res[j, ], w = rot[j, ])
  }
  return(azi_res)
}
