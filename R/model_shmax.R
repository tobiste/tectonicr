#' @title Azimuth between two points
#'
#' @description Calculate initial bearing (or forward azimuth/direction) to go
#' from point \code{p1} to point \code{p2} following great-circle arc on a
#' sphere.
#'
#' @param p,q latitude/longitude of start and end point(s).
#' Can be vectors of two numbers or a matrix of 2 columns (latitude, longitude)
#' @details This formula is for the initial bearing (sometimes referred to as
#' forward azimuth) which if followed in a straight line along a great-circle
#' arc will take you from the start point to the end point:
#' \eqn{\theta = \arctan2 (\sin \Delta\lambda *
#' \cos\psi_2, \cos\psi_1 \sin\psi_1-\sin\psi_1 \cos\psi_2 \cos\Delta\lambda)}
#' where  \eqn{\psi_1, \lambda_1} is the start point, \eqn{\psi_2},
#' \eqn{\lambda_2} the end point (\eqn{\Delta\lambda} is the difference in
#' longitude)
#' @references http://www.movable-type.co.uk/scripts/latlong.html
#' @return Azimuth in degrees
#' @export
#' @examples
#' p <- c(35, 45) # Baghdad
#' q <- c(35, 135) # Osaka
#' get_azimuth(p, q)
get_azimuth <- function(p, q) {
  stopifnot(is.numeric(p) & is.numeric(q))

  p.lat <- p[1]
  p.lon <- p[2]
  q.lat <- q[1]
  q.lon <- q[2]

  delta.lon <- q.lon - p.lon

  # azimuth
  theta <- atan2d(
    sind(delta.lon) * cosd(q.lat),
    cosd(p.lat) * sind(q.lat) -
      sind(p.lat) * cosd(q.lat) * cosd(delta.lon)
  )

  # Normalize result to a compass bearing (0-360)
  theta <- (theta + 360) %% 360
  return(theta)
}



#'  Modeling directions of maximum horizontal stress
#'
#' Calculate the direction of maximum horizontal stress along great circles,
#' small circles, and loxodromes around the relative plate motion´s Euler pole
#' at a given point or points
#'
#' @author Tobias Stephan
#' @param df \code{data.frame} containing the coordinates of the point(s)
#' (lat, lon)
#' @param euler \code{data.frame} containing the coordinates of the Euler pole
#' for the plate boundary (lat, lon)
#' @references Wdowinski, S., 1998, A theory of intraplate
#'   tectonics: Journal of Geophysical Research: Solid Earth, v. 103, p.
#'   5037-5059, http://dx.doi.org/10.1029/97JB03390.
#' @return An object of class \code{data.frame}
#' \describe{
#'   \item{gc}{azimuth of the modeled maximum horizontal following stress
#'   great circles}
#'   \item{sc}{small circles}
#'   \item{ld.cw}{clockwise loxodromes}
#'   \item{ld.ccw}{counter-clockwise loxodromes}
#'   }
#' @export
#' @examples
#' data("nuvel1")
#' euler <- subset(nuvel1, nuvel1$ID == "na") # North America relative to
#' # Pacific plate
#' point <- data.frame(lat = 45, lon = 20)
#' model_shmax(point, euler)
model_shmax <- function(df, euler) {
  for (i in seq_along(df$lat)) {
    # great circles
    gc <- get_azimuth(
      c(df$lat[i], df$lon[i]), c(euler$lat[1], euler$lon[1]) + 180
    ) %% 180

    # plate vector is perpendicular to the bearing between the point and the
    # Euler pole
    sc <- (
      get_azimuth(c(df$lat[i], df$lon[i]), c(euler$lat[1], euler$lon[1])) - 90 + 180
    ) %% 180


    # counterclockwise loxodrome
    ld.ccw <- (
      get_azimuth(c(df$lat[i], df$lon[i]), c(euler$lat[1], euler$lon[1])) - 45 + 180
    ) %% 180


    # clockwise loxodrome
    ld.cw <- (
      get_azimuth(c(df$lat[i], df$lon[i]), c(euler$lat[1], euler$lon[1])) + 45 + 180
    ) %% 180

    x.i <- data.frame(
      gc, sc, ld.ccw, ld.cw
    )
    if (i == 1) {
      x <- x.i
    } else {
      x <- rbind(x, x.i)
    }
  }
  return(x)
}
