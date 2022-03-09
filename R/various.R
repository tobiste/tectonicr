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

#' Quantize World Stress Map Quality Ranking
#'
#' Relates the categorical quality ranking of the World Stress Map (A, B, C, D)
#' to numeric values of the precision of each measurement.
#'
#' @param x Either a string or a character vector of WSM quality ranking
#' @return \code{"integer"} or vector of type \code{"integer"}
#' @references Heidbach, O.; Barth, A.; Müller, B.; Reinecker, J.;
#' Stephansson, O.; Tingay, M.; Zang, A. (2016). WSM quality
#' ranking scheme, database description and analysis guidelines for stress
#' indicator. World Stress Map Technical Report 16-01, GFZ German Research
#' Centre for Geosciences. \doi{10.2312/wsm.2016.001}
#' @export
#' @examples
#' quantise_wsm_quality(c("A", "B", "C", "D", NA))
#' data("wsm2016")
#' quantise_wsm_quality(wsm2016$quality)
quantise_wsm_quality <- function(x) {
  azi.std <- c()
  for (i in seq_along(x)) {
    azi.std[i] <- ifelse(x[i] == "A", 15,
      ifelse(x[i] == "B", 20,
        ifelse(x[i] == "C", 25,
          ifelse(x[i] == "D", 40, NA)
        )
      )
    )
  }
  return(azi.std)
}



#' Distance from plate boundary
#'
#' Distance of data points from the nearest plate boundary in degree
#'
#' @param x,pb \code{sf} objects of the data points and the plate boundary
#' geometries in the geographical coordinate system
#' @param ep \code{data.frame} of the geographical coordinates of the Euler pole (\code{lat}, \code{lon})
#' @param tangential Logical. \code{TRUE} for tangential boundaries and
#' \code{FALSE} (the default) for inward and outward boundaries.
#' @return Numeric vector of the (great circle) distances in degree
#' @export
distance_from_pb <- function(x, ep, pb, tangential = FALSE){
  stopifnot(inherits(x, "sf") & inherits(pb, "sf") & is.data.frame(ep))

  x.por <- geographical_to_PoR(ep, x)
  pb.por <- geographical_to_PoR(ep, pb)

  pb.coords <- sf::st_coordinates(pb.por)
  colnames(coords) <- c('lon', 'lat')

  x.coords <- sf::st_coordinates(x.por)
  colnames(coords) <- c('lon', 'lat')

  dist <- c()
  for(i in seq_along(x.coords$lon)) {
    if(tangential) {
      # latitudinal distance for tangential plate boundaries
      dist[i] <- min(abs(pb.coords$lat[i] - x.coords$lat[i]))
    } else {
      # azimuthal distance for inward/outward plate boundaries
      dist[i] <- min(abs(pb.coords$lon[i] - x.coords$lon[i]))
    }
  }

  dist
}


