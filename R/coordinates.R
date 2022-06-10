#' @title Coordinate Correction
#'
#' @description Corrects the longitudes or latitudes to value between -180.0 and
#' 180.0 or -90 and 90 degree
#' @param x Longitude(s) or latitude(s) in degree
#' @return numeric
#' @name coordinate_mod
#' @examples
#' longitude_modulo(-361 + 5 * 360)
#' latitude_modulo(-91 + 5 * 180)
NULL

#' @rdname coordinate_mod
#' @export
longitude_modulo <- function(x) {
  stopifnot(is.numeric(x))
  # longitude.mod <- (longitude %% 360 + 540) %% 360 - 180
  (x + 540) %% 360 - 180
}

#' @rdname coordinate_mod
#' @export
latitude_modulo <- function(x) {
  stopifnot(is.numeric(x))
  asind(sind(x))
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
  # if (abs(lat) > 90) {
  #   lat <- latitude_modulo(lat)
  #   lon <- longitude_modulo(lon + 180)
  # }
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

#' PoR coordinate reference system
#'
#' Helper function to create the reference system transformed in Euler pole
#' coordinate
#'
#' @param x \code{data.frame} of the geographical coordinates of the Euler pole
#' @details The PoR coordinate reference system is oblique transformation of the
#' geographical coordinate system with the Euler pole coordinates being the the
#' translation factors.
PoR_crs <- function(x) {
  stopifnot(is.numeric(x$lat) | is.numeric(x$lon))
  if (x$lat < 0) {
    x$lat <- -x$lat
    x$lon <- longitude_modulo(x$lon + 180)
  }

  sf::st_crs(
    paste0(
      "+proj=ob_tran +o_proj=longlat +datum=WGS84 +o_lat_p=",
      x$lat,
      " +o_lon_p=",
      x$lon
    )
  )
}


#' Conversion between PoR to geographical coordinate system
#'
#' Transform spherical objects from PoR to geographical coordinate system and
#' vice versa.
#'
#' @param x \code{sf} object of the data points in geographical or PoR
#' coordinate system
#' @param ep \code{data.frame} of the geographical coordinates of the Euler
#' pole (\code{lat}, \code{lon})
#' @return \code{sf} object of the data points in the transformed geographical
#' or PoR coordinate system
#' @details The PoR coordinate reference system is oblique transformation of the
#' geographical coordinate system with the Euler pole coordinates being the the
#' translation factors.
#' @importFrom dplyr %>%
#' @importFrom sf st_crs st_as_sf st_set_crs st_transform st_wrap_dateline
#' @examples
#' data("nuvel1")
#' euler <- subset(nuvel1, nuvel1$plate.rot == "na") # North America relative to Pacific plate
#'
#' data("wsm2016")
#' example.geo <- sf::st_set_crs(
#'   sf::st_as_sf(wsm2016[1:10, ], coords = c("lon", "lat")), "EPSG:4326"
#' )
#'
#' example.por <- geographical_to_PoR(example.geo, euler)
#' PoR_to_geographical(example.por, euler)
#' @name por_transformation
NULL

#' @rdname por_transformation
#' @export
PoR_to_geographical <- function(x, ep) {
  stopifnot(inherits(x, "sf") & is.data.frame(ep))

  crs.wgs84 <-
    sf::st_crs(
      "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
    )

  crs.ep <- PoR_crs(ep)

  suppressMessages(
    suppressWarnings(
      x.por <- x %>%
        sf::st_set_crs(crs.wgs84) %>%
        sf::st_transform(crs.ep) %>%
        sf::st_set_crs(crs.wgs84) %>%
        sf::st_wrap_dateline()
    )
  )
  return(x.por)
}

#' @rdname por_transformation
#' @export
geographical_to_PoR <- function(x, ep) {
  stopifnot(inherits(x, "sf") & is.data.frame(ep))

  crs.wgs84 <-
    sf::st_crs(
      "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
    )

  crs.ep <- PoR_crs(ep)

  suppressMessages(
    suppressWarnings(
      x.geo <- x %>%
        sf::st_set_crs(crs.ep) %>%
        sf::st_transform(crs.wgs84) %>%
        sf::st_set_crs(crs.ep) %>%
        sf::st_wrap_dateline(
          options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180")
        )
    )
  )
  return(x.geo)
}
