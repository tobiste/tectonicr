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
  # longitude.mod <- (longitude %% 360 + 540) %% 360 - 180
  (x + 540) %% 360 - 180
}

#' @rdname coordinate_mod
#' @export
latitude_modulo <- function(x) {
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
#' @seealso [cartesian_to_spherical()] and [spherical_to_cartesian()] for
#' conversions to spherical coordinates
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
  stopifnot(length(n) == 3)
  r <- sqrt(n[1]^2 + n[2]^2 + n[3]^2)
  p <- c()
  p[1] <- asind(n[3] / r) # lat
  p[2] <- atan2d(n[2], n[1]) # lon
  # if (abs(lat) > 90) {
  #   lat <- latitude_modulo(lat)
  #   lon <- longitude_modulo(lon + 180)
  # }
  p
}

#' @rdname coordinates
#' @export
geographical_to_cartesian <- function(p) {
  stopifnot(length(p) == 2)
  x <- c()
  x[1] <- cosd(p[1]) * cosd(p[2])
  x[2] <- cosd(p[1]) * sind(p[2])
  x[3] <- sind(p[1])
  x
}

#' @rdname coordinates
#' @export
geographical_to_spherical <- function(p) {
  geographical_to_cartesian(p) |> cartesian_to_spherical()
}


#' @title Coordinate Transformations
#'
#' @description Converts vector between Cartesian and spherical coordinate
#' systems
#' @param n Cartesian coordinates (x, y, z) as three-column vector
#' @param p Spherical coordinates (colatitude, azimuth) as two-column vector
#' @return Functions return a (2- or 3-dimensional) vector representing a
#' point in the requested coordinate system.
#' @seealso [cartesian_to_geographical()] and [geographical_to_cartesian()] for
#' conversions to geographical coordinates
#' @examples
#' n <- c(1, -2, 3)
#' cartesian_to_spherical(n)
#' p <- c(50, 10)
#' spherical_to_cartesian(p)
#' @name coordinates2
NULL

#' @rdname coordinates2
#' @export
cartesian_to_spherical <- function(n) {
  stopifnot(length(n) == 3)
  r <- sqrt(n[1]^2 + n[2]^2 + n[3]^2)
  p <- c()
  p[1] <- acosd(n[3] / r) # colat/inclination
  p[2] <- atan2d(n[2], n[1]) # azimuth
  p
}

#' @rdname coordinates2
#' @export
spherical_to_cartesian <- function(p) {
  stopifnot(length(p) == 2)
  x <- c()
  x[1] <- sind(p[1]) * cosd(p[2])
  x[2] <- sind(p[1]) * sind(p[2])
  x[3] <- cosd(p[1])
  x
}

#' @rdname coordinates2
#' @export
spherical_to_geographical <- function(p) {
  spherical_to_cartesian(p) |> cartesian_to_geographical()
}


#' PoR coordinate reference system
#'
#' Create the reference system transformed in Euler pole
#' coordinate
#'
#' @param x \code{"data.frame"} or \code{"euler.pole"} object containing the
#' geographical coordinates of the Euler pole
#'
#' @details The PoR coordinate reference system is oblique transformation of the
#' geographical coordinate system with the Euler pole coordinates being the the
#' translation factors.
#'
#' @returns Object of class `crs`
#'
#' @importFrom sf st_crs
#'
#' @seealso [sf::st_crs()]
#'
#' @export
#'
#' @examples
#' data("nuvel1")
#' por <- subset(nuvel1, nuvel1$plate.rot == "na") # North America relative to Pacific plate
#' PoR_crs(por)
PoR_crs <- function(x) {
  stopifnot(is.data.frame(x) | is.euler(x))
  x <- as.data.frame(x)
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


#' Conversion between PoR to geographical coordinate system using quaternions
#'
#' Helper function for the transformation from PoR to geographical coordinate
#' system or vice versa
#'
#' @param x,PoR two-column vectors containing the lat and lon coordinates
#'
#' @name por_transformation_quat
#'
#' @returns two-element numeric vector
#'
#' @examples
#' ep.geo <- c(20, 33)
#' q.geo <- c(10, 45)
#' q.por <- geographical_to_PoR_quat(q.geo, ep.geo)
#' q.por
#' PoR_to_geographical_quat(q.por, ep.geo)
NULL

#' @name por_transformation_quat
#' @export
geographical_to_PoR_quat <- function(x, PoR) {
  p <- geographical_to_cartesian(x)
  euler_y <- euler_pole(0, 1, 0, angle = 90 - PoR[1], geo = FALSE)
  euler_z <- euler_pole(0, 0, 1, angle = 180 - PoR[2], geo = FALSE)
  qy <- euler_to_Q4(euler_y)
  qz <- euler_to_Q4(euler_z)
  qq <- product_Q4(q1 = qz, q2 = qy)
  p_trans <- rotation_Q4(q = qq, p = p) |>
    cartesian_to_geographical()
  p_trans[2] <- longitude_modulo(p_trans[2] + 180)
  return(p_trans)
}

#' @name por_transformation_quat
#' @export
PoR_to_geographical_quat <- function(x, PoR) {
  x[2] <- longitude_modulo(x[2] - 180)
  p_trans <- geographical_to_cartesian(x)
  PoR_yt <- euler_pole(0, -1, 0, angle = 90 - PoR[1], geo = FALSE)
  PoR_zt <- euler_pole(0, 0, -1, angle = 180 - PoR[2], geo = FALSE)
  qy <- euler_to_Q4(PoR_yt) # |> conjugate_Q4()
  qz <- euler_to_Q4(PoR_zt) # |> conjugate_Q4()

  product_Q4(q1 = qy, q2 = qz) |>
    rotation_Q4(p = p_trans) |>
    cartesian_to_geographical()
}

#' Conversion between spherical PoR to geographical coordinate system
#'
#' Transformation from spherical PoR to geographical coordinate system and
#' vice versa
#'
#' @param x \code{"data.frame"} containing \code{lat} and \code{lon}
#' coordinates of a point in the geographical CRS or the \code{lat.PoR},
#' \code{lon.PoR}) of the point in the PoR CRS.
#' @param PoR Pole of Rotation. \code{"data.frame"} or object of class \code{"euler.pole"}
#' containing the geographical coordinates of the Euler pole
#' @return \code{"data.frame"} with the transformed coordinates
#' (\code{lat.PoR} and \code{lon.PoR} for PoR CRS,
#' or \code{lat} and \code{lon} for geographical CRS).
#' @name por_transformation_df
#' @examples
#' data("nuvel1")
#' por <- subset(nuvel1, nuvel1$plate.rot == "na") # North America relative to Pacific plate
#' data("san_andreas")
#' san_andreas.por <- geographical_to_PoR(san_andreas, por)
#' head(san_andreas.por)
#' head(PoR_to_geographical(san_andreas.por, por))
NULL

#' @name por_transformation_df
#' @export
geographical_to_PoR <- function(x, PoR) {
  stopifnot(is.data.frame(PoR) | is.euler(PoR))
  ep.geo <- c(PoR$lat, PoR$lon)
  lat.PoR <- lon.PoR <- c()

  for (i in seq_along(x$lon)) {
    x_por.i <- geographical_to_PoR_quat(c(x$lat[i], x$lon[i]), PoR = ep.geo)
    lat.PoR[i] <- x_por.i[1]
    lon.PoR[i] <- x_por.i[2]
  }
  data.frame(lat.PoR = lat.PoR, lon.PoR = lon.PoR)
}

#' @name por_transformation_df
#' @export
PoR_to_geographical <- function(x, PoR) {
  stopifnot(is.data.frame(PoR) | is.euler(PoR))
  ep.geo <- c(PoR$lat, PoR$lon)
  lat <- lon <- c()

  for (i in seq_along(x$lon.PoR)) {
    x_geo.i <- PoR_to_geographical_quat(c(x$lat.PoR[i], x$lon.PoR[i]), PoR = ep.geo)
    lat[i] <- x_geo.i[1]
    lon[i] <- x_geo.i[2]
  }
  data.frame(lat = lat, lon = lon)
}


#' PoR coordinates
#'
#' Retrieve the PoR equivalent coordinates of an object
#'
#' @param x \code{sf} or \code{data.frame} containing lat and lon coordinates
#' (\code{lat}, \code{lon})
#' @param PoR Pole of Rotation. \code{"data.frame"} or object of class \code{"euler.pole"}
#' containing the geographical coordinates of the Euler pole
#' @return \code{data.frame} with the PoR coordinates
#' (\code{lat.PoR}, \code{lon.PoR})
#' @export
#' @examples
#' data("nuvel1")
#' por <- subset(nuvel1, nuvel1$plate.rot == "na") # North America relative to Pacific plate
#' data("san_andreas")
#' san_andreas.por_sf <- PoR_coordinates(san_andreas, por)
#' head(san_andreas.por_sf)
#' san_andreas.por_df <- PoR_coordinates(sf::st_drop_geometry(san_andreas), por)
#' head(san_andreas.por_df)
PoR_coordinates <- function(x, PoR) {
  if (is.data.frame(x)) {
    # x <- sf::st_as_sf(x, coords = c("lon", "lat"))
    geographical_to_PoR(x, PoR)
  } else {
    x |>
      tectonicr::geographical_to_PoR_sf(PoR = PoR) |>
      sf::st_coordinates() |>
      sf::st_drop_geometry() |>
      rename("lon.PoR" = "X", "lat.PoR" = "Y")
  }
}


#' Conversion between PoR to geographical coordinate reference system of raster
#' data
#'
#' Helper function to transform raster data set from PoR to geographical
#' coordinates
#'
#' @param x \code{"SpatRaster"} or \code{"RasterLayer"}
#' @param PoR Pole of Rotation. \code{"data.frame"} or object of class \code{"euler.pole"}
#' containing the geographical coordinates of the Euler pole
#' @returns "SpatRaster"
#' @importFrom terra crs project rast
#' @importFrom methods extends
#' @name raster_transformation
NULL

#' @rdname raster_transformation
geographical_to_PoR_raster <- function(x, PoR) {
  if (methods::extends(class(x), "BasicRaster")) {
    x <- terra::rast(x)
  }
  stopifnot(inherits(x, "SpatRaster"))
  crs.wgs84 <- "epsg:4326"
  crs.ep <- PoR_crs(PoR)
  terra::crs(x) <- crs.ep$wkt
  x.por <- terra::project(x, crs.wgs84)
  terra::crs(x.por) <- crs.ep$wkt
  return(x.por)
}

#' @rdname raster_transformation
PoR_to_geographical_raster <- function(x, PoR) {
  if (methods::extends(class(x), "BasicRaster")) {
    x <- terra::rast(x)
  }
  stopifnot(inherits(x, "SpatRaster"))
  crs.wgs84 <- "epsg:4326"
  crs.PoR <- PoR_crs(PoR)
  terra::crs(x) <- crs.wgs84
  x.geo <- terra::project(x, crs.PoR$wkt)
  terra::crs(x.geo) <- crs.wgs84
  return(x.geo)
}


#' Conversion between PoR to geographical coordinates of spatial data
#'
#' Transform spatial objects from PoR to geographical coordinate reference
#' system and vice versa.
#'
#' @param x \code{sf}, \code{SpatRast}, or \code{Raster*} object of the data
#' points in geographical or PoR coordinate system
#' @param PoR Pole of Rotation. \code{"data.frame"} or object of class \code{"euler.pole"}
#' containing the geographical coordinates of the Euler pole
#' @return \code{sf} or \code{SpatRast} object of the data points in the
#' transformed geographical or PoR coordinate system
#' @details The PoR coordinate reference system is oblique transformation of the
#' geographical coordinate system with the Euler pole coordinates being the
#' translation factors.
#' @importFrom sf st_crs st_as_sf st_set_crs st_transform st_wrap_dateline
#' @importFrom methods extends
#' @examples
#' data("nuvel1")
#' PoR <- subset(nuvel1, nuvel1$plate.rot == "na") # North America relative to Pacific plate
#' data("san_andreas")
#' san_andreas.por <- geographical_to_PoR_sf(san_andreas, PoR)
#' PoR_to_geographical_sf(san_andreas.por, PoR)
#' @name por_transformation_sf
NULL

#' @rdname por_transformation_sf
#' @export
PoR_to_geographical_sf <- function(x, PoR) {
  if (methods::extends(class(x), "BasicRaster") | inherits(x, "SpatRaster")) {
    x.por <- geographical_to_PoR_raster(x, PoR)
  } else {
    crs.wgs84 <- sf::st_crs("epsg:4326")
    crs.PoR <- PoR_crs(PoR)
    suppressMessages(
      suppressWarnings(
        x.por <- x |>
          sf::st_set_crs(crs.wgs84) |>
          sf::st_transform(crs.PoR) |>
          sf::st_set_crs(crs.wgs84) |>
          sf::st_wrap_dateline()
      )
    )
  }
  return(x.por)
}

#' @rdname por_transformation_sf
#' @export
geographical_to_PoR_sf <- function(x, PoR) {
  if (methods::extends(class(x), "BasicRaster") | inherits(x, "SpatRaster")) {
    x.geo <- geographical_to_PoR_raster(x, PoR)
  } else {
    crs.wgs84 <- sf::st_crs("epsg:4326")
    crs.PoR <- PoR_crs(PoR)
    suppressMessages(
      suppressWarnings(
        x.geo <- x |>
          sf::st_set_crs(crs.PoR) |>
          sf::st_transform(crs.wgs84) |>
          sf::st_set_crs(crs.PoR) |>
          sf::st_wrap_dateline(
            options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180")
          )
      )
    )
  }
  return(x.geo)
}
