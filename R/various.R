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
#' @param x Either a string or a character vector of WSM quality ranking
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
#' quantise_wsm_quality(san_andreas$quality[1:20])
quantise_wsm_quality <- function(x) {
  unc <- c()
  for (i in seq_along(x)) {
    unc[i] <- ifelse(x[i] == "A", 15,
      ifelse(x[i] == "B", 20,
        ifelse(x[i] == "C", 25,
          ifelse(x[i] == "D", 40, NA)
        )
      )
    )
  }
  return(unc)
}

#' @title Normalize angular distance on a sphere distance
#'
#' @description  Helper function to express angular distance on the sphere in
#' the range of
#' -180 to 180 degrees
#'
#' @param x numeric, angular distance (in degrees)
#' @returns numeric vector
#' @keywords internal
distance_mod <- function(x) {
  x <- abs(x)
  for (i in seq_along(x)) {
    while (x[i] > 180) x[i] <- 360 - (x[i] %% 360)
  }
  x
}

#' Distance from plate boundary
#'
#' Absolute distance of data points from the nearest plate boundary in degree
#'
#' @param x,pb `sf` objects of the data points and the plate boundary
#' geometries in the geographical coordinate system
#' @param ep \code{data.frame} of the geographical coordinates of the Euler pole
#' (`lat`, `lon`)
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
#' @importFrom sf st_geometry_type st_cast st_coordinates
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
#'   x = san_andreas, ep = na_pa, pb = plate_boundary, tangential = TRUE
#' )
#' head(res)
#'
#' res.km <- distance_from_pb(
#'   x = san_andreas, ep = na_pa, pb = plate_boundary, tangential = TRUE, km = TRUE
#' )
#' head(res.km)
distance_from_pb <- function(x, ep, pb, tangential = FALSE, km = FALSE, ...) {
  stopifnot(
    inherits(x, "sf") &
      inherits(pb, "sf") & is.data.frame(ep) &
      is.logical(tangential) &
      is.logical(km)
  )

  x.por <- geographical_to_PoR(x, ep)
  pb.por <- geographical_to_PoR(pb, ep) %>%
    sf::st_cast(to = "LINESTRING") %>%
    smoothr::densify(...)

  pb.coords <- sf::st_coordinates(pb.por)
  x.coords <- sf::st_coordinates(x.por)

  dist <- c()
  for (i in seq_along(x.coords[, 1])) {
    delta.lat <- distance_mod(pb.coords[, 2] - x.coords[i, 2])
    delta.lon <- distance_mod(pb.coords[, 1] - x.coords[i, 1])

    if (tangential) {
      # latitudinal differences for tangential plate boundaries and
      # select the one with the closest longitude
      q <- which.min(delta.lon)
      dist[i] <- delta.lat[q] # latitudinal difference in degree
      if (km) {
        dist[i] <- deg2rad(dist[i]) * earth_radius() # great circle distance
      }
    } else {
      # longitudinal differences for inward/outward plate boundaries and
      # select the one with the closest latitude
      q <- which.min(delta.lat)
      dist[i] <- delta.lon[q] # longitudinal difference in degree
      if (km) {
        dist[i] <- deg2rad(dist[i]) * earth_radius() * cosd(x.coords[i, 2]) # small circle distance
      }
    }
  }
  dist
}
