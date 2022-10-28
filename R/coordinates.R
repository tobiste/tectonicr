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
  geographical_to_cartesian(p) %>% cartesian_to_spherical()
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
  spherical_to_cartesian(p) %>% cartesian_to_geographical()
}


#' PoR coordinate reference system
#'
#' Create the reference system transformed in Euler pole
#' coordinate
#'
#' @param x \code{data.frame} of the geographical coordinates of the Euler pole
#' @details The PoR coordinate reference system is oblique transformation of the
#' geographical coordinate system with the Euler pole coordinates being the the
#' translation factors.
#' @export
#' @examples
#' data("nuvel1")
#' euler <- subset(nuvel1, nuvel1$plate.rot == "na") # North America relative to Pacific plate
#' PoR_crs(euler)
PoR_crs <- function(x) {
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



#' Conversion between spherical PoR to geographical coordinate system
#'
#' Transformation from spherical PoR to geographical coordinate system and
#' vice versa
#'
#' @param x \code{sf} or \code{data.frame} containing (co)lat and lon coordinates
#' (\code{lat}, \code{lon}) of the points to be transformed
#' @param ep \code{data.frame} of the geographical coordinates of the Euler pole
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
PoR_to_geographical2 <- function(x, ep, spherical = TRUE) {
  ep.geo <- c(ep$lat, ep$lon)
  lat <- lon <- c()
  for(i in seq_along(x$lon)){
    if(spherical){
      x_geo.i <- PoR_to_geographical_vec(c(x$colat.PoR[i], x$lon.PoR[i]), ep = ep.geo, spherical = TRUE)
    } else {
      x_geo.i <- PoR_to_geographical_vec(c(x$lat.PoR[i], x$lon.PoR[i]), ep = ep.geo, spherical = FALSE)
    }
    lat[i] <- x_geo.i[1]
    lon[i] <- x_geo.i[2]
  }
  data.frame(lat = lat, lon = lon)
}

#' @rdname por_conversion_df
#' @export
geographical_to_PoR2 <- function(x, ep, spherical = TRUE) {
  ep.geo <- c(ep$lat, ep$lon)
  colat.PoR <- lon.PoR <- c()
  for(i in seq_along(x$lon)){
    x_por.i <- geographical_to_PoR_vec(c(x$lat[i], x$lon[i]), ep = ep.geo, spherical)
    colat.PoR[i] <- x_por.i[1]
    lon.PoR[i] <- x_por.i[2]
  }
  if(spherical){
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
#' @param x,ep two-column vectors containing the (co)lat and lon coordinates
#' @param spherical logical. Whether x or the return are in spherical
#' coordinates
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
PoR_to_geographical_vec <- function(x, ep, spherical = TRUE) {
  if(spherical){
    x.por.cart <- spherical_to_cartesian(x)
  } else {
    x.por.cart <- geographical_to_cartesian(x)
  }
  x.cart <- t(rotmat_z(180 - ep[2])) %*% t(rotmat_y(90 - ep[1])) %*% x.por.cart
  cartesian_to_geographical(x.cart)

}

#' @rdname por_conversion_vec
geographical_to_PoR_vec <- function(x, ep, spherical = TRUE) {
  x.cart <- geographical_to_cartesian(x)
  x.cart.por <- rotmat_y(90 - ep[1]) %*% rotmat_z(180 - ep[2]) %*% x.cart
  if(spherical){
    cartesian_to_spherical(x.cart.por)
  } else {
    cartesian_to_geographical(x.cart.por)
  }
}


#' PoR coordinates
#'
#' Retrieve the PoR equivalent coordinates of an object
#'
#' @param x \code{sf} or \code{data.frame} containing lat and lon coordinates
#' (\code{lat}, \code{lon})
#' @param ep \code{data.frame} of the geographical coordinates of the Euler pole
#' (\code{lat}, \code{lon})
#' @return \code{data.frame} with the PoR coordinates
#' (\code{lat.PoR}, \code{lon.PoR})
#' @export
#' @examples
#' data("nuvel1")
#' euler <- subset(nuvel1, nuvel1$plate.rot == "na") # North America relative to Pacific plate
#' data("san_andreas")
#' san_andreas.por <- PoR_coordinates(san_andreas, euler)
#' head(san_andreas.por)
PoR_coordinates <- function(x, ep) {
  .Deprecated("geographical_to_PoR2")

  if (is.data.frame(x)) {
    x <- sf::st_as_sf(x, coords = c("lon", "lat"))
  }
  if (is.data.frame("x")) {
    x <- sf::st_as_sf(x, coords = c("lon", "lat"))
  }
  x %>%
    tectonicr::geographical_to_PoR(ep = ep) %>%
    sf::st_coordinates() %>%
    as.data.frame() %>%
    rename("lon.PoR" = "X", "lat.PoR" = "Y")
}


#' Conversion between PoR to geographical coordinate reference system of raster
#' data
#'
#' Helper function to transform raster data set from PoR to geographical
#' coordinates
#'
#' @param x \code{"SpatRaster"} or \code{"RasterLayer"}
#' @param ep \code{data.frame} of the geographical coordinates of the Euler pole
#' (\code{lat}, \code{lon})
#' @returns "SpatRaster"
#' @importFrom terra crs project rast
#' @importFrom methods extends
#' @name raster_transformation
NULL

#' @rdname raster_transformation
geographical_to_PoR_raster <- function(x, ep) {
  if (methods::extends(class(x), "BasicRaster")) {
    x <- terra::rast(x)
  }
  stopifnot(is.data.frame(ep) & inherits(x, "SpatRaster"))
  crs.wgs84 <- "epsg:4326"
  crs.ep <- PoR_crs(ep)
  terra::crs(x) <- crs.ep$wkt
  x.por <- terra::project(x, crs.wgs84)
  terra::crs(x.por) <- crs.ep$wkt
  return(x.por)
}

#' @rdname raster_transformation
PoR_to_geographical_raster <- function(x, ep) {
  if (methods::extends(class(x), "BasicRaster")) {
    x <- terra::rast(x)
  }
  stopifnot(is.data.frame(ep) & inherits(x, "SpatRaster"))
  crs.wgs84 <- "epsg:4326"
  crs.ep <- PoR_crs(ep)
  terra::crs(x) <- crs.wgs84
  x.geo <- terra::project(x, crs.ep$wkt)
  terra::crs(x.geo) <- crs.wgs84
  return(x.geo)
}


#' Conversion between PoR to geographical coordinates
#'
#' Transform spherical objects from PoR to geographical coordinate reference
#' system and vice versa.
#'
#' @param x \code{sf}, \code{SpatRast}, or \code{Raster*} object of the data
#' points in geographical or PoR coordinate system
#' @param ep \code{data.frame} of the geographical coordinates of the Euler
#' pole (\code{lat}, \code{lon})
#' @return \code{sf} or \code{SpatRast} object of the data points in the
#' transformed geographical or PoR coordinate system
#' @details The PoR coordinate reference system is oblique transformation of the
#' geographical coordinate system with the Euler pole coordinates being the
#' translation factors.
#' @importFrom magrittr %>%
#' @importFrom sf st_crs st_as_sf st_set_crs st_transform st_wrap_dateline
#' @importFrom methods extends
#' @examples
#' data("nuvel1")
#' euler <- subset(nuvel1, nuvel1$plate.rot == "na") # North America relative to Pacific plate
#' data("san_andreas")
#' san_andreas.por <- geographical_to_PoR(san_andreas, euler)
#' PoR_to_geographical(san_andreas.por, euler)
#' @name por_transformation
NULL

#' @rdname por_transformation
#' @export
PoR_to_geographical <- function(x, ep) {
  if (methods::extends(class(x), "BasicRaster") | inherits(x, "SpatRaster")) {
    x.por <- geographical_to_PoR_raster(x, ep)
  } else {
    stopifnot(inherits(x, "sf") & is.data.frame(ep))
    crs.wgs84 <- sf::st_crs("epsg:4326")
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
  }
  return(x.por)
}

#' @rdname por_transformation
#' @export
geographical_to_PoR <- function(x, ep) {
  if (methods::extends(class(x), "BasicRaster") | inherits(x, "SpatRaster")) {
    x.geo <- geographical_to_PoR_raster(x, ep)
  } else {
    stopifnot(inherits(x, "sf") & is.data.frame(ep))
    crs.wgs84 <- sf::st_crs("epsg:4326")
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
  }
  return(x.geo)
}
