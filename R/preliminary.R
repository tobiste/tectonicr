#' Azimuth relative to PoR
#'
#' Transform azimuth from geographical to PoR system according to Wdowinski 1998
#'
#' @param x \code{"data.frame"}. coordinates of data point (lat, lon) and the azimuth
#' @param ep \code{"data.frame"}. coordinates of Euler pole (PoR)
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
PoR_azimuth <- function(x, ep) {
  # .Deprecated("PoR_shmax")
  # convert latitude to colatitude and to radians
  x <- data.frame(lat = x$lat, lon = x$lon, azi = x$azi) * pi / 180
  ep <- data.frame(lat = ep$lat, lon = ep$lon) * pi / 180

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
  stopifnot(is.data.frame(df) & is.data.frame(euler))
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
#' @param ep \code{data.frame} of the geographical coordinates of the Euler pole
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
#' res2 <- stress_matrix(x = san_andreas, ep = na_pa, tangential = TRUE, positive = FALSE)
#' head(res2)
#' }
NULL

#' @rdname stressstrain
displacement_vector <- function(x, ep, tangential = FALSE, positive = TRUE) {
  if (positive) {
    u <- abs(ep$angle)
  } else {
    u <- -abs(ep$angle)
  }
  x.por <- geographical_to_PoR(x, ep) %>%
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
    crs = PoR_crs(ep)
  )
}

#' @rdname stressstrain
stress_matrix <- function(x, ep, tangential = FALSE, positive = FALSE, v = .25, E = 50) {
  if (positive) {
    u <- abs(ep$angle)
  } else {
    u <- -abs(ep$angle)
  }
  x.por <- geographical_to_PoR(x, ep) %>%
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
    crs = PoR_crs(ep)
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
#' geographical_to_PoR_quat(q.por, ep.geo)
geographical_to_PoR_quat <- function(x, euler) {
  p <- geographical_to_cartesian(x)
  wy <- deg2rad(90-euler[1])
  wz <- deg2rad(180-euler[2])

  y <- c(0, 1, 0)
  z <- c(0, 0, 1)
  qy <- as.quaterion(wy, y) #
  qz <- as.quaterion(wz, z) #
  # qy <- cos(wy / 2) + sin(wy / 2) * y
  # qz <- cos(wz / 2) + sin(wz / 2) * z
  qyt <- as.quaterion(-wy, y) #
  qzt <- as.quaterion(-wz, z) #
  # qyt <- cos(wy / 2) - sin(wy / 2) * y
  # qzt <- cos(wz / 2) - sin(wz / 2) * z

  q <- qy * qz
  qt <- qzt * qyt
  #(q * p * t(q)) %>%
  (q * p * qt) %>%
    cartesian_to_geographical()
}

PoR_to_geographical_quat <- function(x, euler) {
  p_trans <- geographical_to_cartesian(x)
  wy <- deg2rad(90 - euler[1])
  wz <- deg2rad(180 - euler[2])

  y <- c(0, 1, 0)
  z <- c(0, 0, 1)
  qy <- as.quaterion(wy, y)
  qz <- as.quaterion(wz, z)

  q <- qz * qy
  (q * p_trans * t(q)) %>%
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

as.quaterion <- function(angle, vector) {
  cos(angle / 2) + (vector * sin(angle / 2))
}


# microbenchmark::microbenchmark(
#   geographical_to_PoR_quat(q.geo, ep.geo),
#   geographical_to_PoR_vec(q.geo, ep.geo, spherical = FALSE),
#   geographical_to_PoR(data.frame(lat = q.geo[1], lon = q.geo[2]) %>% sf::st_as_sf(coords=c('lon', 'lat')), euler = data.frame(lat = ep.geo[1], lon = ep.geo[2]))
# )
