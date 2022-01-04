#' Azimuth between two points
#'
#' Calculate initial bearing (or forward azimuth/direction) to go from point \code{p1} 
#' to point \code{p2} following great-circle arc on a sphere.
#'
#' @author Tobias Stephan
#' @param p latitude/longitude of start point(s). Can be a vector of two numbers or a 
#' matrix of 2 columns (first one is latitude, second is longitude)
#' @param q as above
#' @references http://www.movable-type.co.uk/scripts/latlong.html
#' C.F.F. Karney, 2013. Algorithms for geodesics, J. Geodesy 87: 43-55. https://dx.doi.org/10.1007/s00190-012-0578-z
#' @return azimuth in degrees
#' @export
#' @examples
#' p <- c(35, 45) # Baghdad
#' q <- c(35, 135) # Osaka
#' get_azimuth(p, q)
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
