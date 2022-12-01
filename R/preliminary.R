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
  (azi.por * 180 / pi + 180) %% 180
}

#' Azimuth conversion from PoR to geographical coordinate reference system
#'
#' Helper function to convert PoR azimuths into geographical azimuths
#'
#' @param df \code{data.frame} containing the PoR-equivalent coordinates of the
#' point(s) (\code{lat.PoR}, \code{lon.PoR}) and the PoR-equivalent azimuth
#' (\code{azi.PoR})
#' @param euler \code{data.frame} containing the geographical location of
#' the Euler pole (\code{lat}, \code{lon})
#' @param spherical logical
#' @examples
#' \dontrun{
#' data("nuvel1")
#' # North America relative to Pacific plate:
#' euler <- subset(nuvel1, nuvel1$plate.rot == "na")
#' data("san_andreas")
#' head(san_andreas$azi)
#' azi.PoR <- PoR_shmax(san_andreas, euler)
#' res.PoR <- data.frame(
#'   azi.PoR = azi.PoR,
#'   geographical_to_PoR2(san_andreas, euler,
#'     spherical = FALSE
#'   )
#' )
#' res.geo <- PoR2geo_shmax(res.PoR, euler, spherical = FALSE)
#' head(res.geo)
#' }
PoR2geo_shmax <- function(df, euler, spherical = FALSE) {
  stopifnot(is.data.frame(df))
  northpole.por <- geographical_to_PoR2(
    data.frame(lat = 0, lon = 0),
    euler,
    spherical = spherical
  ) %>%
    as.numeric()

  beta <- c()
  if (spherical) {
    p <- cbind(df$colat.PoR, df$lon.PoR)
  } else {
    p <- cbind(df$lat.PoR, df$lon.PoR)
  }

  for (i in 1:nrow(p)) {
    beta[i] <- get_azimuth(
      p[i, ],
      northpole.por
    )
  }

  (df$azi.PoR - beta + 180) %% 180
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
  for (i in 1:length(x.por[, 1])) {
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

    for (i in 1:length(x.por[, 1])) {
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

    for (i in 1:length(x.por[, 1])) {
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
  Q4_to_QScVec(c(q0, q1, q2, q3, q4), normalize)
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