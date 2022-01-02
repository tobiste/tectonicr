#' Azimuth between two points
#'
#' Calculate initial bearing (or forward azimuth) which is followed along a
#' straight line along a great-circle arc between the start point and the end
#' point
#'
#' @author Copyright (C) 2021 Tobias Stephan
#' @param p two-column numeric vector containing latitude and longitude of start point
#' @param q two-column numeric vector containing latitude and longitude of end point
#' @references http://www.movable-type.co.uk/scripts/latlong.html
#' @return azimuth between two points at start point p
#' @export
#' @examples
#' p <- c(35, 45) # Baghdad
#' q <- c(35, 135) # Osaka
#' get_azimuth(p1, p2)
#'
get_azimuth <- function(p, q){
  p.lat <- p[1]; p.lon <- p[2]
  q.lat <- q[1]; q.lon <- q[2]

  delta.lon <- q.lon - p.lon

  #azimuth
  theta <- pracma::atan2d(
    pracma::sind(delta.lon) * pracma::cosd(q.lat),
    pracma::cosd(p.lat) * pracma::sind(q.lat)  -
      pracma::sind(p.lat) * pracma::cosd(q.lat) * pracma::cosd(delta.lon))

  # Normalize result to a compass bearing (0-360)
  theta <- (theta + 360) %% 360
  return(theta)
}
