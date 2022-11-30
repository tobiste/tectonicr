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



#' Conversion between PoR to geographical coordinate system using quaternions
#'
#' Helper function for the transformation from PoR to geographical coordinate
#' system or vice versa
#'
#' @param x,euler two-column vectors containing the lat and lon coordinates
#' @examples
#' ep.geo <- c(20, 33)
#' q.geo <- c(10, 45)
#' q.por <- geographical_to_PoR_quat(q.geo, ep.geo)
#' q.por
#' geographical_to_PoR_vec(q.geo, ep.geo, spherical = FALSE)
#' geographical_to_PoR(data.frame(lat = q.geo[1], lon = q.geo[2]) %>% sf::st_as_sf(coords = c("lon", "lat")), euler = data.frame(lat = ep.geo[1], lon = ep.geo[2]))
#' PoR_to_geographical_quat(q.por, ep.geo)
geographical_to_PoR_quat <- function(x, euler) {
  p <- geographical_to_cartesian(x)
  angle_y <- deg2rad(90 - euler[1])
  angle_z <- deg2rad(180 - euler[2])
  axis_y <- c(0, 1, 0)
  axis_z <- c(0, 0, 1)
  qy <- euler_to_Q4(angle = angle_y, axis = axis_y)
  qz <- euler_to_Q4(angle = angle_z, axis = axis_z)
  qq <- product_Q4(q1 = qz, q2 = qy)
  p_trans = rotation_Q4(q = qq, p = p) %>%
    cartesian_to_geographical()
  p_trans[2] <- longitude_modulo(p_trans[2] + 180)
  return(p_trans)
}

PoR_to_geographical_quat <- function(x, euler) {
  x[2] <- longitude_modulo(x[2] - 180)
  p_trans <- geographical_to_cartesian(x)

  angle_y <- deg2rad(90 - euler[1])
  angle_z <- deg2rad(180 - euler[2])
  axis_y <- c(0, 1, 0)
  axis_z <- c(0, 0, 1)
  qy <- euler_to_Q4(angle = angle_y, axis = -axis_y) #%>% conjugate_Q4()
  qz <- euler_to_Q4(angle = angle_z, axis = -axis_z) #%>% conjugate_Q4()

  product_Q4(q1 = qy, q2 = qz) %>%
    rotation_Q4(p = p_trans) %>%
    cartesian_to_geographical()
}

geographical_to_PoR3 <- function(x, euler) {
  ep.geo <- c(euler$lat, euler$lon)
  lat.PoR <- lon.PoR <- c()

  for (i in seq_along(x$lon)) {
    x_por.i <- geographical_to_PoR_quat(c(x$lat[i], x$lon[i]), euler = ep.geo)
    lat.PoR[i] <- x_por.i[1]
    lon.PoR[i] <- x_por.i[2]
  }
  data.frame(lat.PoR = lat.PoR, lon.PoR = lon.PoR)
}


#' Quaternion from Euler angle-axis representation for rotations
#'
#' @param angle angle in radians
#' @param axis,p three-column vector of unit length
#' @param q,q1,q2 quaternion (list)
# as.quaternion <- function(angle, axis) {
#   cos(angle / 2) + (axis * sin(angle / 2))
# }
euler_to_Q4 <- function(angle, axis, normalize = FALSE) {
  Sc <- cos(angle / 2)
  Vec <- axis * sin(angle / 2)
  q <- (list(Sc = c(Sc), Vec = Vec))
  class(q) <- "quaternion"
  if (normalize) {
    q <- normalize_Q4(q)
  }
  return(q)
}


normalize_Q4 <- function(q) {
  q4 <- Q4_2(q)
  q.norm <- q4 / sqrt(q4[1]^2 + q4[2]^2 + q4[3]^2 + q4[4]^2)
  q.norm <- list(Sc = q.norm[1], Vec = c(q.norm[2], q.norm[3], q.norm[4]))
  class(q.norm) <- "quaternion"
  return(q.norm)
}

#' Product of quaternions
#'
#' Multiplication of two quaternions is not commutative.
#' Concatenation of two rotations R1 followed by R2
#'
#' @param q1,q2 quaternion (list). first rotation R1 expressed by q1 followed by second rotation R2 expressed by q2
#' @returns quaternion
product_Q4 <- function(q1, q2, normalize = FALSE) {
  Sc <- q2$Sc * q1$Sc - (q2$Vec %*% q1$Vec)
  Vec <- q1$Sc * q2$Vec + q2$Sc * q1$Vec + vcross(q2$Vec, q1$Vec)
  q <- list(Sc = c(Sc), Vec = Vec)
  class(q) <- "quaternion"
  if (normalize) {
    q <- normalize_Q4(q)
  }
  return(q)
}

#' Euler angle/axis from quaternion
#' @param q quaternion (list)
#' @returns list with the angle (in radians) and the axis (Cartesian coordinates).
Q4_to_euler <- function(q) {
  stopifnot(is.Q4(q))
  angle <- 2 * acos(q$Sc)
  axis <- q$Vec / sin(angle / 2)
  list(angle = angle, axis = axis)
}

conjugate_Q4 <- function(q, normalize = FALSE) {
  q <- list(Sc = q$Sc, Vec = -q$Vec)
  class(q) <- "quaternion"
  if (normalize) {
    q <- normalize_Q4(q)
  }
  return(q)
}

Q4 <- function(q) {
  stopifnot(is.Q4(q))
  q$Sc + q$Vec
}

is.Q4 <- function(q) {
  class(q) == "quaternion"
}

#' Rotation of a vector by a quaternion
#' @param q quaternion (list)
#' @param p three-column vector (Cartesian coordinates) of unit length
#' @returns three-column vector (Cartesian coordinates) of unit length
rotation_Q4 <- function(q, p) {
  stopifnot(is.Q4(q))
  #Q4(q) * p * Q4(conjugate_Q4(q))

  # Rodrigues Formula
  q.euler <- Q4_to_euler(q)
  #p * cos(q.euler$angle) + vcross(q.euler$axis, p) * sin(q.euler$angle) + q.euler$axis * as.numeric(q.euler$axis %*% p) * (1 - cos(q.euler$axis))
  p + vcross(q.euler$axis, p) * sin(q.euler$angle) + vcross(q.euler$axis, vcross(q.euler$axis, p)) * 2 * (sin(q.euler$angle/2))^2
}

rotmat_to_Q4 <- function(x, normalize = FALSE) {
  stopifnot(is.matrix(x))
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

Q4_to_rotmat <- function(x, normalize = TRUE) {
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

QScVec_to_Q4 <- function(x) {
  stopifnot(x$Sc != 0)
  x.euler <- Q4_to_euler(x)
  q <- c()
  q[1] <- x$Sc
  q[2] <- x.euler$axis[1] * sin(x.euler$angle / 2)
  q[3] <- x.euler$axis[2] * sin(x.euler$angle / 2)
  q[4] <- x.euler$axis[3] * sin(x.euler$angle / 2)
  q
}

Q4_to_QScVec <- function(x, normalize = FALSE){
  stopifnot(x[1] != 0)
  angle <- 2 * acos(x[1])
  q1 <- x[2]/sin(angle/2)
  q2 <- x[3]/sin(angle/2)
  q3 <- x[4]/sin(angle/2)

  q <- list(Sc = x[1], Vec = c(q1, q2, q3))
  class(q) <- "quaternion"
  if (normalize) {
    q <- normalize_Q4(q)
  }
  return(q)
}

microbenchmark::microbenchmark(
  geographical_to_PoR_quat(q.geo, ep.geo),
  geographical_to_PoR_vec(q.geo, ep.geo, spherical = FALSE),
  geographical_to_PoR(data.frame(lat = q.geo[1], lon = q.geo[2]) %>% sf::st_as_sf(coords=c('lon', 'lat')), euler = data.frame(lat = ep.geo[1], lon = ep.geo[2]))
)
