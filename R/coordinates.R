#' @title Coordinate Correction
#'
#' @description Corrects the longitudes or latitudes to value between -180.0 and
#' 180.0 or -90 and 90 degree
#'
#' @param x Longitude(s) or latitude(s) in degrees
#'
#' @returns numeric
#'
#' @name coordinate_mod
#'
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
#'
#' @param n Cartesian coordinates (x, y, z) as vector
#' @param p Geographical coordinates (latitude, longitude) as vector
#'
#' @returns Functions return a (2- or 3-dimensional) vector representing a
#' point in the requested coordinate system.
#'
#' @seealso [cartesian_to_spherical()] and [spherical_to_cartesian()] for
#' conversions to spherical coordinates
#'
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
  p <- numeric(2)
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
  x <- numeric(3)
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
#'
#' @param n Cartesian coordinates (x, y, z) as three-column vector
#' @param p Spherical coordinates (colatitude, azimuth) as two-column vector
#'
#' @returns Functions return a (2- or 3-dimensional) vector representing a
#' point in the requested coordinate system.
#'
#' @seealso [cartesian_to_geographical()] and [geographical_to_cartesian()] for
#' conversions to geographical coordinates
#'
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
  p <- numeric(2)
  p[1] <- acosd(n[3] / r) # colat/inclination
  p[2] <- atan2d(n[2], n[1]) # azimuth
  p
}

#' @rdname coordinates2
#' @export
spherical_to_cartesian <- function(p) {
  stopifnot(length(p) == 2)
  x <- numeric(3)
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
  # if (x$lat < 0) {
  #   x$lat <- -x$lat
  #   x$lon <- longitude_modulo(x$lon + 180)
  # }
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
#' @keywords internal
#'
#' @returns two-element numeric vector
NULL

#' @name por_transformation_quat
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

geographical_to_PoR_helper <- function(lat, lon, PoR) {
  geographical_to_PoR_quat(c(lat, lon), PoR = PoR)
}

PoR_to_geographical_helper <- function(lat.PoR, lon.PoR, PoR) {
  PoR_to_geographical_quat(c(lat.PoR, lon.PoR), PoR = PoR)
}


#' Conversion between spherical PoR to geographical coordinate system
#'
#' Transformation from spherical PoR to geographical coordinate system and
#' vice versa
#'
#' @param x Can be either a \code{"data.frame"} containing \code{lat} and \code{lon}
#' coordinates of a point in the geographical CRS or the \code{lat.PoR},
#' \code{lon.PoR}) of the point in the PoR CRS,
#' a two-column matrix containing the lat and lon coordinates,
#' a `sf` object, or a `raster` object.
#' @param PoR Pole of Rotation. \code{"data.frame"} or object of class \code{"euler.pole"}
#' containing the geographical coordinates of the Euler pole
#'
#' @returns object of same type of `x` with the transformed coordinates. If `x`
#' was a `data.frame`, transformed coordinates are named \code{lat.PoR} and \code{lon.PoR} for PoR CRS,
#' or \code{lat} and \code{lon} for geographical CRS).
#'
#' @name por_transformation
#'
#' @examples
#' data("nuvel1")
#' por <- subset(nuvel1, nuvel1$plate.rot == "na") # North America relative to Pacific plate
#' data("san_andreas")
#' san_andreas.por <- geographical_to_PoR(san_andreas, por)
#' head(san_andreas.por)
#' head(PoR_to_geographical(san_andreas.por, por))
NULL

#' @rdname por_transformation
#' @export
geographical_to_PoR <- function(x, PoR) {
  stopifnot(is.data.frame(PoR) | is.euler(PoR))

  if (methods::extends(class(x), "BasicRaster") | inherits(x, "SpatRaster")) {
    geographical_to_PoR_raster(x, PoR)
  } else if (inherits(x, "sf")) {
    geographical_to_PoR_sf(x, PoR)
  } else if (is.data.frame(x)) {
    geographical_to_PoR_df(x, PoR)
  } else if (is.matrix(x)) {
    geographical_to_PoR_quat(x, PoR)
  } else {
    NULL
  }
}

#' @rdname por_transformation
#' @export
PoR_to_geographical <- function(x, PoR) {
  stopifnot(is.data.frame(PoR) | is.euler(PoR))

  if (methods::extends(class(x), "BasicRaster") | inherits(x, "SpatRaster")) {
    PoR_to_geographical_raster(x, PoR)
  } else if (inherits(x, "sf")) {
    PoR_to_geographical_sf(x, PoR)
  } else if (is.data.frame(x)) {
    PoR_to_geographical_df(x, PoR)
  } else if (is.matrix(x)) {
    PoR_to_geographical_quat(x, PoR)
  } else {
    NULL
  }
}


#' Conversion between spherical PoR to geographical coordinate system of data.frames
#'
#' Transformation from spherical PoR to geographical coordinate system and
#' vice versa
#'
#' @param x \code{"data.frame"} containing \code{lat} and \code{lon}
#' coordinates of a point in the geographical CRS or the \code{lat.PoR},
#' \code{lon.PoR}) of the point in the PoR CRS.
#' @param PoR Pole of Rotation. \code{"data.frame"} or object of class \code{"euler.pole"}
#' containing the geographical coordinates of the Euler pole
#'
#' @returns \code{"data.frame"} with the transformed coordinates
#' (\code{lat.PoR} and \code{lon.PoR} for PoR CRS,
#' or \code{lat} and \code{lon} for geographical CRS).
#' @keywords internal
#' @name por_transformation_df
NULL


#' @name por_transformation_df
geographical_to_PoR_df <- function(x, PoR) {
  stopifnot(c("lat", "lon") %in% colnames(x))
  ep.geo <- c(PoR$lat, PoR$lon)
  lat.PoR <- lon.PoR <- numeric(nrow(x))


  x_por <- mapply(geographical_to_PoR_helper, x$lat, x$lon, MoreArgs = list(PoR = ep.geo))
  data.frame(lat.PoR = x_por[1, ], lon.PoR = x_por[2, ])
}

#' @name por_transformation_df
PoR_to_geographical_df <- function(x, PoR) {
  stopifnot(c("lat.PoR", "lon.PoR") %in% colnames(x))
  ep.geo <- c(PoR$lat, PoR$lon)
  lat <- lon <- numeric(nrow(x))

  x_geo <- mapply(PoR_to_geographical_helper, x$lat.PoR, x$lon.PoR, MoreArgs = list(PoR = ep.geo))
  data.frame(lat = x_geo[1, ], lon = x_geo[2, ])
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
#'
#' @returns terra "SpatRaster" object
#'
#' @importFrom terra crs project rast
#' @importFrom methods extends
#' @keywords internal
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


#' Conversion between PoR to geographical coordinates of sf data
#'
#' Transform spatial objects from PoR to geographical coordinate reference
#' system and vice versa.
#'
#' @param x \code{sf}, \code{SpatRast}, or \code{Raster*} object of the data
#' points in geographical or PoR coordinate system
#' @param PoR Pole of Rotation. \code{"data.frame"} or object of class \code{"euler.pole"}
#' containing the geographical coordinates of the Euler pole
#'
#' @return \code{sf} or \code{SpatRast} object of the data points in the
#' transformed geographical or PoR coordinate system
#'
#' @details The PoR coordinate reference system is oblique transformation of the
#' geographical coordinate system with the Euler pole coordinates being the
#' translation factors.
#'
#' @importFrom sf st_crs st_as_sf st_set_crs st_transform st_wrap_dateline
#' @importFrom methods extends
#' @keywords internal
#' @name por_transformation_sf
NULL

#' @rdname por_transformation_sf
PoR_to_geographical_sf <- function(x, PoR) {
  crs.wgs84 <- sf::st_crs("epsg:4326")
  crs.PoR <- PoR_crs(PoR)
  suppressMessages(
    suppressWarnings(
      x |>
        sf::st_set_crs(crs.wgs84) |>
        sf::st_transform(crs.PoR) |>
        sf::st_set_crs(crs.wgs84) |>
        sf::st_wrap_dateline()
    )
  )
}

#' @rdname por_transformation_sf
geographical_to_PoR_sf <- function(x, PoR) {
  crs.wgs84 <- sf::st_crs("epsg:4326")
  crs.PoR <- PoR_crs(PoR)
  suppressMessages(
    suppressWarnings(
      x |>
        sf::st_set_crs(crs.PoR) |>
        sf::st_transform(crs.wgs84) |>
        sf::st_set_crs(crs.PoR) |>
        sf::st_wrap_dateline(
          options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180")
        )
    )
  )
}





#' Coordinates of the Pole of Rotation Reference System
#'
#' Retrieve the PoR equivalent coordinates of an object
#'
#' @param x \code{sf} or \code{data.frame} containing lat and lon coordinates
#' (\code{lat}, \code{lon})
#' @param PoR Pole of Rotation. \code{"data.frame"} or object of class \code{"euler.pole"}
#' containing the geographical coordinates of the Euler pole
#'
#' @return [PoR_coordinates()] returns \code{data.frame} with the PoR coordinates
#' (\code{lat.PoR}, \code{lon.PoR}).
#'
#' @export
#'
#' @examples
#' data("nuvel1")
#' por <- subset(nuvel1, nuvel1$plate.rot == "na") # North America relative to Pacific plate
#' data("san_andreas")
#'
#' # coordinates from sf object
#' san_andreas.por_sf <- PoR_coordinates(san_andreas, por)
#' head(san_andreas.por_sf)
#'
#' # coordinates from data.frame
#' san_andreas.por_df <- PoR_coordinates(sf::st_drop_geometry(san_andreas), por)
#' head(san_andreas.por_df)
PoR_coordinates <- function(x, PoR) {
  if (is.data.frame(x)) {
    # x <- sf::st_as_sf(x, coords = c("lon", "lat"))
    geographical_to_PoR_df(x, PoR)
  } else {
    x |>
      geographical_to_PoR_sf(PoR = PoR) |>
      sf::st_coordinates() |>
      sf::st_drop_geometry() |>
      dplyr::rename("lon.PoR" = "X", "lat.PoR" = "Y")
  }
}

#' Distance to Pole of Rotation
#'
#' Retrieve the (angular) distance to the PoR (Euler pole).
#'
#' @inheritParams PoR_coordinates
#' @param FUN function to calculate the great-circle distance.
#' [orthodrome()], [haversine()] (the default), or [vincenty()].
#' @return numeric vector
#'
#' @export
#' @examples
#' data("nuvel1")
#' por <- subset(nuvel1, nuvel1$plate.rot == "na") # North America relative to Pacific plate
#' data("san_andreas")
#'
#' # distance form sf object
#' PoR_distance(san_andreas, por)
#'
#' # distance form data.frame
#' PoR_distance(sf::st_drop_geometry(san_andreas), por)
#' PoR_distance(sf::st_drop_geometry(san_andreas), por, FUN = orthodrome)
#' PoR_distance(sf::st_drop_geometry(san_andreas), por, FUN = vincenty)
PoR_distance <- function(x, PoR, FUN = orthodrome) {
  if (inherits(x, "sf")) {
    x_crds <- sf::st_coordinates(x)
    lat <- x_crds[, 2]
    lon <- x_crds[, 1]
  } else {
    lat <- x$lat
    lon <- x$lon
  }
  do.call(FUN, list(lat1 = deg2rad(lat), lon1 = deg2rad(lon), lat2 = deg2rad(PoR$lat), lon2 = deg2rad(PoR$lon))) |>
    rad2deg()
}
# PoR_distance <- function(x, PoR) {
#   res <- PoR_coordinates(x, PoR)
#   90 - abs(res$lat.PoR)
# }
