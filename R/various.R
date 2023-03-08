#' Class for Central Position of Spoke Marker
#'
#' position subclass \code{"center_spoke"} to center \code{ggplot::geom_spoke()}
#' marker at its origin
#'
#' @noRd
position_center_spoke <- function() PositionCenterSpoke #

#' @title  Centrically aligned geom_spoke marker
#'
#' @description \code{"position"} subclass "center_spoke" to center
#' \code{ggplot::geom_spoke()} marker at its origin
#' @export
#' @source \url{https://stackoverflow.com/questions/55474143/how-to-center-geom-spoke-around-their-origin}
#' @importFrom ggplot2 ggproto Position
PositionCenterSpoke <- ggplot2::ggproto("PositionCenterSpoke", ggplot2::Position,
  compute_panel = function(self, data, params, scales) {
    data$x <- 2 * data$x - data$xend
    data$y <- 2 * data$y - data$yend
    data$radius <- 2 * data$radius
    data
  }
)

#' Numerical values to World Stress Map Quality Ranking
#'
#' Assigns numeric values of the precision of each measurement to the
#' categorical quality ranking of the World Stress Map (A, B, C, D).
#'
#' @param regime Either a string or a character vector of WSM quality ranking
#' @return \code{"integer"} or vector of type \code{"integer"}
#' @references Heidbach, O., Barth, A., Müller, B., Reinecker, J.,
#' Stephansson, O., Tingay, M., Zang, A. (2016). WSM quality
#' ranking scheme, database description and analysis guidelines for stress
#' indicator. *World Stress Map Technical Report* **16-01**, GFZ German Research
#' Centre for Geosciences. \doi{10.2312/wsm.2016.001}
#' @export
#' @examples
#' quantise_wsm_quality(c("A", "B", "C", "D", NA))
#' data("san_andreas")
#' quantise_wsm_quality(san_andreas$quality)
quantise_wsm_quality <- function(regime) {
  as.numeric(sapply(X = regime, FUN = regime2unc))
}

regime2unc <- function(x) {
  c(
    "A" = 15,
    "B" = 20,
    "C" = 25,
    "D" = 40
  )[x]
}

#' Helper function to get normalized distance from plate boundary
#'
#' @param x numeric
#' @seealso [distance_from_pb()]
get_distance_mod <- function(x) {
  while (x > 180) {
    x <- 360 - x
  }
  return(x)
}

#' @title Normalize angular distance on a sphere distance
#'
#' @description  Helper function to express angular distance on the sphere in
#' the range of 0 to 180 degrees
#'
#' @param x numeric, angular distance (in degrees)
#' @returns numeric vector
#' @keywords internal
distance_mod <- function(x) {
  sapply(X = abs(x) %% 360, FUN = get_distance_mod)
}

#' Helper function to Distance from plate boundary
#'
#' @param lon,lat numeric vectors
#' @param pb.coords matrix
#' @param tangential,km logical
#'
#' @seealso [distance_from_pb()]
get_distance <- function(lon, lat, pb.coords, tangential, km) {
  delta.lat <- distance_mod(pb.coords[, 2] - lat)
  delta.lon <- distance_mod(pb.coords[, 1] - lon)

  if (tangential) {
    # latitudinal differences for tangential plate boundaries and
    # select the one with the closest longitude
    q <- which.min(delta.lon)
    dist <- delta.lat[q] # latitudinal difference in degree
    if (km) {
      dist <- deg2rad(dist) * earth_radius() # great circle distance
    }
  } else {
    # longitudinal differences for inward/outward plate boundaries and
    # select the one with the closest latitude
    q <- which.min(delta.lat)
    dist <- delta.lon[q] # longitudinal difference in degree
    if (km) {
      dist <- deg2rad(dist) * earth_radius() * cosd(lat) # small circle distance
    }
  }
  dist
}


#' Distance from plate boundary
#'
#' Absolute distance of data points from the nearest plate boundary in degree
#'
#' @param x,pb `sf` objects of the data points and the plate boundary
#' geometries in the geographical coordinate system
#' @param euler \code{"data.frame"} or object of class \code{"euler.pole"}
#' containing the geographical coordinates of the Pole of Rotation
#' @param tangential Logical. Whether the plate boundary is a tangential
#' boundary (`TRUE`) or an inward and outward boundary (`FALSE`, the
#' default).
#' @param km Logical. Whether the distance is expressed in kilometers
#' (`TRUE`) or in degrees (`FALSE`, the default).
#' @param ... optional arguments passed to [smoothr::densify()]
#' @return Numeric vector of the great circle distances
#' @details The distance to the plate boundary is the longitudinal or
#' latitudinal difference between the data point and the plate boundary
#' (along the closest latitude or longitude) for inward/outward or tangential
#' plate boundaries, respectively.
#' @export
#' @importFrom sf st_geometry st_cast st_coordinates
#' @importFrom magrittr %>%
#' @importFrom smoothr densify
#' @examples
#' data("nuvel1")
#' na_pa <- subset(nuvel1, nuvel1$plate.rot == "na")
#'
#' data("plates")
#' plate_boundary <- subset(plates, plates$pair == "na-pa")
#'
#' data("san_andreas")
#' res <- distance_from_pb(
#'   x = san_andreas, euler = na_pa, pb = plate_boundary, tangential = TRUE
#' )
#' head(res)
#'
#' res.km <- distance_from_pb(
#'   x = san_andreas, euler = na_pa, pb = plate_boundary, tangential = TRUE, km = TRUE
#' )
#' head(res.km)
distance_from_pb <- function(x, euler, pb, tangential = FALSE, km = FALSE, ...) {
  stopifnot(
    inherits(x, "sf"),
    inherits(pb, "sf"),
    is.logical(tangential),
    is.logical(km),
    is.data.frame(euler) | is.euler(euler)
  )

  x.coords <- sf::st_geometry(x) %>%
    geographical_to_PoR_sf(euler = euler) %>%
    sf::st_coordinates()
  pb.coords <-
    sf::st_geometry(pb) %>%
    geographical_to_PoR_sf(euler = euler) %>%
    sf::st_cast(to = "MULTILINESTRING", warn = FALSE) %>%
    sf::st_cast(to = "LINESTRING", warn = FALSE) %>%
    smoothr::densify(...) %>%
    sf::st_coordinates()

  mapply(
    FUN = get_distance,
    lon = x.coords[, 1], lat = x.coords[, 2],
    MoreArgs = list(pb.coords = pb.coords, tangential = tangential, km = km)
  )
}

#' Helper function to get Distance from plate boundary
#'
#' @param lon,lat,pb.bearing numeric vectors
#' @param pb.coords matrix
#' @param tangential logical
#'
#' @seealso [projected_pb_strike()]
get_projected_pb_strike <- function(lon, lat, pb.coords, pb.bearing, tangential) {
  delta.lat <- distance_mod(pb.coords[, 2] - lat)
  delta.lon <- distance_mod(pb.coords[, 1] - lon)

  if (tangential) {
    # select the one with the closest longitude
    q <- which.min(delta.lon)
    x_pb.bearing <- pb.bearing[q]
  } else {
    # select the one with the closest latitude
    q <- which.min(delta.lat)
    x_pb.bearing <- pb.bearing[q]
  }
  x_pb.bearing
}


#' Strike of the plate boundary projected on data point
#'
#' The fault's strike in the PoR CRS projected on the data point along the
#' predicted stress trajectories.
#'
#' @param x,pb `sf` objects of the data points and the plate boundary
#' geometries in the geographical coordinate system
#' @param euler \code{"data.frame"} or object of class \code{"euler.pole"}
#' containing the geographical coordinates of the Euler pole
#' @param tangential Logical. Whether the plate boundary is a tangential
#' boundary (`TRUE`) or an inward and outward boundary (`FALSE`, the
#' default).
#' @param ... optional arguments passed to [smoothr::densify()]
#' @details Useful to calculate the beta angle, i.e. the angle
#' between SHmax direction (in PoR CRS!) and the fault's strike (in PoR CRS).
#' The beta angle is the same in geographical and PoR coordinates.
#' @note The algorithm calculates the great circle bearing between line
#' vertices. Since transform plate boundaries represent small circle lines in
#' the PoR system, this great-circle azimuth is only a approximation of the
#' true (small-circle) azimuth.
#' @return Numeric vector of the strike direction of the plate boundary
#' (in degree)
#' @export
#' @importFrom sf st_cast st_coordinates st_geometry
#' @importFrom magrittr %>%
#' @importFrom smoothr densify
#' @examples
#' data("nuvel1")
#' na_pa <- subset(nuvel1, nuvel1$plate.rot == "na")
#'
#' data("plates")
#' plate_boundary <- subset(plates, plates$pair == "na-pa")
#'
#' data("san_andreas")
#' res <- projected_pb_strike(
#'   x = san_andreas, euler = na_pa, pb = plate_boundary, tangential = TRUE
#' )
#' head(res)
#' head(san_andreas$azi - res) # beta angle
projected_pb_strike <- function(x, euler, pb, tangential = FALSE, ...) {
  stopifnot(
    inherits(x, "sf"),
    inherits(pb, "sf"),
    is.logical(tangential),
    is.data.frame(euler) | is.euler(euler)
  )

  x.coords <- sf::st_geometry(x) %>%
    geographical_to_PoR_sf(euler = euler) %>%
    sf::st_coordinates()
  pb.coords <-
    sf::st_geometry(pb) %>%
    geographical_to_PoR_sf(euler = euler) %>%
    sf::st_cast(to = "MULTILINESTRING", warn = FALSE) %>%
    sf::st_cast(to = "LINESTRING", warn = FALSE) %>%
    smoothr::densify(...) %>%
    sf::st_coordinates()

  pb.bearing <- c()
  for (i in 1:nrow(pb.coords)) {
    if (i == nrow(pb.coords)) {
      pb.bearing[i] <- NA
    } else {
      pb.bearing[i] <- get_azimuth(pb.coords[i, 2], pb.coords[i, 1], pb.coords[i + 1, 2], pb.coords[i + 1, 1])
    }
  }

  mapply(
    FUN = get_projected_pb_strike,
    lon = x.coords[, 1], lat = x.coords[, 2],
    MoreArgs = list(pb.coords = pb.coords, pb.bearing = pb.bearing %% 180, tangential = tangential)
  )
}

#' Quick analysis of a stress data set
#'
#' Returns the converted azimuths, distances to the plate boundary,
#' statistics of the model, and some plots.
#'
#' @param x \code{data.frame} or `sf` object containing the coordinates of the point(s)
#' (\code{lat}, \code{lon}), the direction of
#' \eqn{\sigma_{Hmax}}{SHmax} \code{azi} and its standard deviation
#' \code{unc} (optional)
#' @param euler \code{"data.frame"} or object of class \code{"euler.pole"}
#' containing the geographical coordinates of the Euler pole
#' @param type Character. Type of plate boundary (optional). Can be
#' \code{"out"}, \code{"in"}, \code{"right"}, or
#' \code{"left"} for outward, inward, right-lateral, or left-lateral
#' moving plate boundaries, respectively. If \code{"none"} (the default), only
#' the PoR-equivalent azimuth is returned.
#' @param pb (optional) `sf` object of the plate boundary geometries in the geographical
#' coordinate system
#' @param plot (logical). Whether to produce a plot additional to output.
#' @param ... optional arguments to [distance_from_pb()]
#' @return list. results of the coordinate and azimuth conversion, deviation,
#' misfit to predicted stress direction and, if given, distance to tested
#' plate boundary as well as the normalized Chi-squared test statistic.
#' @export
#' @seealso [PoR_shmax()], [distance_from_pb()], [norm_chisq()], [PoR_plot()]
#' @examples
#' data("nuvel1")
#' na_pa <- subset(nuvel1, nuvel1$plate.rot == "na")
#'
#' data("plates")
#' plate_boundary <- subset(plates, plates$pair == "na-pa")
#'
#' data("san_andreas")
#'
#' stress_analysis(san_andreas, na_pa, type = "right", plate_boundary, plot = FALSE)
stress_analysis <- function(x, euler, type = c("none", "in", "out", "right", "left"), pb, plot = TRUE, ...) {
  type <- match.arg(type)
  stopifnot(is.logical(plot))
  tangential <- ifelse(type %in% c("right", "left"), TRUE, FALSE)
  res <- PoR_shmax(x, euler, type)
  res <- cbind(res, PoR_coordinates(x, euler))
  if (!missing(pb)) {
    res$distance <- distance_from_pb(x, euler, pb, tangential, ...)
  }
  nchisq <- norm_chisq(res$azi.PoR, res$prd, x$unc)

  if (plot) {
    PoR_plot(azi = res$azi.PoR, distance = res$distance, unc = x$unc, regime = x$regime, prd = res$prd)
  }

  list(result = res, norm_chisq = nchisq)
}
