#' Transformations into PoR spherical Coordinates
#'
#' Helper functions for coordinate transformations
#'
#' @param x angle
#' @return matrix
#' @name transform_matrices
NULL

#' @rdname transform_matrices
tmat1 <- function(x) {
  stopifnot(is.numeric(x))
  x <- as.numeric(x) * pi/180
  cbind(
    c(1, 0, 0),
    c(0, cos(x), sin(x)),
    c(0, -sin(x), cos(x))
  )
}

#' @rdname transform_matrices
tmat2 <- function(x) {
  stopifnot(is.numeric(x))
  x <- as.numeric(x) * pi/180
  cbind(
    c(cos(x), 0, -sin(x)),
    c(0, 1, 0),
    c(sin(x), 0, cos(x))
  )
}

#' @rdname transform_matrices
tmat3 <- function(x){
  stopifnot(is.numeric(x))
  x <- as.numeric(x) * pi/180
  cbind(
    c(cos(x), sin(x), 0),
    c(-sin(x), cos(x), 0),
    c(0, 0, 1)
  )
}

#' Conversion between PoR to geographical coordinate system
#'
#' Transformation from PoR to geographical coordinate system or vice versa
#'
#' @param ep two-colum Numeric vector containing the geographic location of the Euler
#' pole (latitude, longitude)
#' @param x two-colum vector of the Point in the PoR system (colatitude, azimuth) or the geographical coordinate system (latitude, longitude)
#' @return two-colum vector of the transformed latitude, longtitude or colatitude, azimuth in geographical or PoR coordinate system, respectively
#' @name por_conversion
#' @references Wdowinski, S., 1998, A theory of intraplate
#'   tectonics: Journal of Geophysical Research: Solid Earth, v. 103, p.
#'   5037-5059, \doi{10.1029/97JB03390}.
#' @examples
#' \dontrun{
#' ep.geo <- c(52, 50)
#' q.geo <- c(10, 45)
#' q.por <- geographical_to_PoR2(q.por, ep.geo)
#' PoR_to_geographical2(q.por, ep.geo)
#' }
NULL

#' @rdname por_conversion
PoR_to_geographical2 <- function(x, ep) {
  stopifnot(is.numeric(ep) & is.numeric(x))
  x.por.cart <- spherical_to_cartesian(x)
  x.cart <- tmat3(180 - ep[2]) %*%  tmat2(90-ep[1]) %*% x.por.cart
  cartesian_to_geographical(x.cart)
}

#' @rdname por_conversion
geographical_to_PoR2 <- function(x, ep){
  stopifnot(is.numeric(ep) & is.numeric(x))
  x.cart <-  geographical_to_cartesian(x)
  x.cart.por <- tmat2(90-ep[1]) %*% tmat3(180-ep[2]) %*% x.cart
  cartesian_to_spherical(x.cart.por)
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