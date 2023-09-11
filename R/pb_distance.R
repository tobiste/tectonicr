#' Helper function to get normalized distance from plate boundary
#'
#' @param x numeric
#' @seealso [distance_from_pb()]
get_distance_mod <- function(x) {
  while (abs(x) > 180) {
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
  # sapply(X = abs(x) %% 360, FUN = get_distance_mod)
  sapply(X = x, FUN = get_distance_mod)
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
    q <- which.min(abs(delta.lon))
    dist <- delta.lat[q] # latitudinal difference in degree
    if (km) {
      dist <- deg2rad(dist) * earth_radius() # great circle distance
    }
  } else {
    # longitudinal differences for inward/outward plate boundaries and
    # select the one with the closest latitude
    q <- which.min(abs(delta.lat))
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
#' @param PoR Pole of Rotation. \code{"data.frame"} or object of class \code{"euler.pole"}
#' containing the geographical coordinates of the Pole of Rotation
#' @param tangential Logical. Whether the plate boundary is a tangential
#' boundary (`TRUE`) or an inward and outward boundary (`FALSE`, the
#' default).
#' @param km Logical. Whether the distance is expressed in kilometers
#' (`TRUE`) or in degrees (`FALSE`, the default).
#' @param ... optional arguments passed to [smoothr::densify()]
#'
#' @return Numeric vector of the great circle distances
#'
#' @details The distance to the plate boundary is the longitudinal or
#' latitudinal difference between the data point and the plate boundary
#' (along the closest latitude or longitude) for inward/outward or tangential
#' plate boundaries, respectively.
#'
#' @export
#'
#' @references Wdowinski, S. (1998). A theory of intraplate tectonics. Journal
#' of Geophysical Research: Solid Earth, 103(3), 5037–5059.
#' http://dx.doi.org/10.1029/97JB03390
#'
#' @importFrom sf st_geometry st_cast st_coordinates
#' @importFrom smoothr densify
#'
#' @examples
#' \donttest{
#' data("nuvel1")
#' na_pa <- subset(nuvel1, nuvel1$plate.rot == "na")
#'
#' data("plates")
#' plate_boundary <- subset(plates, plates$pair == "na-pa")
#'
#' data("san_andreas")
#' res <- distance_from_pb(
#'   x = san_andreas, PoR = na_pa, pb = plate_boundary, tangential = TRUE
#' )
#' head(res)
#'
#' res.km <- distance_from_pb(
#'   x = san_andreas, PoR = na_pa, pb = plate_boundary, tangential = TRUE, km = TRUE
#' )
#' range(res.km)
#' }
distance_from_pb <- function(x, PoR, pb, tangential = FALSE, km = FALSE, ...) {
  stopifnot(
    inherits(x, "sf"),
    inherits(pb, "sf"),
    is.logical(tangential),
    is.logical(km),
    is.data.frame(PoR) | is.euler(PoR)
  )

  x.coords <- sf::st_geometry(x) |>
    geographical_to_PoR_sf(PoR = PoR) |>
    sf::st_coordinates()
  pb.coords <-
    sf::st_geometry(pb) |>
    geographical_to_PoR_sf(PoR = PoR) |>
    sf::st_cast(to = "MULTILINESTRING", warn = FALSE) |>
    sf::st_cast(to = "LINESTRING", warn = FALSE) |>
    smoothr::densify(...) |>
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
    q <- which.min(abs(delta.lon))
    x_pb.bearing <- pb.bearing[q]
  } else {
    # select the one with the closest latitude
    q <- which.min(abs(delta.lat))
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
#' @param PoR Pole of rotation. \code{"data.frame"} or object of class \code{"euler.pole"}
#' containing the geographical coordinates of the Euler pole
#' @param tangential Logical. Whether the plate boundary is a tangential
#' boundary (`TRUE`) or an inward and outward boundary (`FALSE`, the
#' default).
#' @param ... optional arguments passed to [smoothr::densify()]
#'
#' @details Useful to calculate the beta angle, i.e. the angle
#' between SHmax direction (in PoR CRS!) and the fault's strike (in PoR CRS).
#' The beta angle is the same in geographical and PoR coordinates.
#'
#' @note The algorithm calculates the great circle bearing between line
#' vertices. Since transform plate boundaries represent small circle lines in
#' the PoR system, this great-circle azimuth is only a approximation of the
#' true (small-circle) azimuth.
#'
#' @return Numeric vector of the strike direction of the plate boundary
#' (in degree)
#'
#' @export
#'
#' @importFrom sf st_cast st_coordinates st_geometry
#' @importFrom smoothr densify
#'
#' @examples
#' data("nuvel1")
#' na_pa <- subset(nuvel1, nuvel1$plate.rot == "na")
#'
#' data("plates")
#' plate_boundary <- subset(plates, plates$pair == "na-pa")
#'
#' data("san_andreas")
#' res <- projected_pb_strike(
#'   x = san_andreas, PoR = na_pa, pb = plate_boundary, tangential = TRUE
#' )
#' head(res)
#' head(san_andreas$azi - res) # beta angle
projected_pb_strike <- function(x, PoR, pb, tangential = FALSE, ...) {
  stopifnot(
    inherits(x, "sf"),
    inherits(pb, "sf"),
    is.logical(tangential),
    is.data.frame(PoR) | is.euler(PoR)
  )

  x.coords <- sf::st_geometry(x) |>
    geographical_to_PoR_sf(PoR = PoR) |>
    sf::st_coordinates()
  pb.coords <-
    sf::st_geometry(pb) |>
    geographical_to_PoR_sf(PoR = PoR) |>
    sf::st_cast(to = "MULTILINESTRING", warn = FALSE) |>
    sf::st_cast(to = "LINESTRING", warn = FALSE) |>
    smoothr::densify(...) |>
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
