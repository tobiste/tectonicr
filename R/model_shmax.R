#' @title Azimuth Between two Points
#'
#' @description Calculate initial bearing (or forward azimuth/direction) to go
#' from point \code{p1} to point \code{p2} following great circle arc on a
#' sphere.
#'
#' @param a,b Latitude/longitude of start and end point(s).
#' Can be vectors of two numbers or a matrix of 2 columns (latitude, longitude)
#' @details This formula is for the initial bearing (sometimes referred to as
#' forward azimuth) which if followed in a straight line along a great circle
#' arc will take you from the start point to the end point:
#' \deqn{\theta = \arctan2 (\sin \Delta\lambda *
#' \cos\psi_2, \cos\psi_1 \sin\psi_1-\sin\psi_1 \cos\psi_2 \cos\Delta\lambda)}
#' where  \eqn{\psi_1, \lambda_1} is the start point, \eqn{\psi_2},
#' \eqn{\lambda_2} the end point (\eqn{\Delta\lambda} is the difference in
#' longitude)
#' @references \url{http://www.movable-type.co.uk/scripts/latlong.html}
#' @return Azimuth in degrees
#' @export
#' @examples
#' berlin <- c(52.517, 13.4) # Berlin
#' tokyo <- c(35.7, 139.767) # Tokyo
#' get_azimuth(berlin, tokyo)
get_azimuth <- function(a, b) {
  stopifnot(is.numeric(a) & is.numeric(b))

  # convert deg into rad
  phi1 <- pi / 180 * a[1]
  phi2 <- pi / 180 * b[1]

  d.lambda <- (b[2] - a[2]) * (pi / 180)

  y <- sin(d.lambda) * cos(phi2)
  x <- cos(phi1) * sin(phi2) -
    sin(phi1) * cos(phi2) * cos(d.lambda)
  theta <- atan2(y, x) * 180 / pi

  # Normalize result to a compass bearing (0-360)
  (theta + 360) %% 360
}



#' @title Theoretical Direction of Maximum Horizontal Stress
#'
#' @description Models the direction of maximum horizontal stress \eqn{\sigma_\text{Hmax}}{SHmax} along
#' great circles, small circles, and loxodromes at a given point or points
#' according to the relative plate motion
#'
#' @author Tobias Stephan
#' @param df \code{data.frame} containing the coordinates of the point(s)
#' (lat, lon)
#' @param euler \code{data.frame} containing the coordinates of the Euler pole
#' for the plate boundary (lat, lon)
#' @details \eqn{\sigma_\text{Hmax}}{SHmax} following *great circles* is the
#' (initial) bearing between the given point and the pole of relative plate
#' motion. \eqn{\sigma_\text{Hmax}}{SHmax} along *small circles*, clockwise, and
#' counter-clockwise *loxodromes* is 90\eqn{^{\circ}}{ degree}, +45\eqn{^{\circ}}{ degree}, and 135\eqn{^{\circ}}{ degree} (-45\eqn{^{\circ}}{ degree})
#'  to this great circle bearing, respectively.
#' @return \code{data.frame}
#' \describe{
#'   \item{gc}{Azimuth of modeled \eqn{\sigma_\text{Hmax}}{SHmax} following
#'   great circles}
#'   \item{sc}{Small circles}
#'   \item{ld.cw}{Clockwise loxodromes}
#'   \item{ld.ccw}{Counter-clockwise loxodromes}
#'  }
#' @seealso [misfit_shmax()] to compute the deviation of the modeled direction
#'  from the observed direction of \eqn{\sigma_\text{Hmax}}{SHmax}.
#' @export
#' @examples
#' data("nuvel1")
#' # North America relative to Pacific plate:
#' euler <- subset(nuvel1, nuvel1$ID == "na")
#'
#' # the point where we mant to model the SHmax direction:
#' point <- data.frame(lat = 45, lon = 20)
#'
#' model_shmax(point, euler)
model_shmax <- function(df, euler) {
  sc <- c()
  gc <- c()
  ld.cw <- c()
  ld.ccw <- c()
  for (i in seq_along(df$lat)) {
    # great circles
    gc[i] <- get_azimuth(
      c(df$lat[i], df$lon[i]), c(euler$lat[1], euler$lon[1]) + 180
    ) %% 180

    # plate vector is perpendicular to the bearing between the point and the
    # Euler pole
    sc[i] <- (
      get_azimuth(c(df$lat[i], df$lon[i]), c(euler$lat[1], euler$lon[1])) - 90 + 180
    ) %% 180


    # counterclockwise loxodrome
    ld.ccw[i] <- (
      get_azimuth(c(df$lat[i], df$lon[i]), c(euler$lat[1], euler$lon[1])) - 45 + 180
    ) %% 180


    # clockwise loxodrome
    ld.cw[i] <- (
      get_azimuth(c(df$lat[i], df$lon[i]), c(euler$lat[1], euler$lon[1])) + 45 + 180
    ) %% 180
  }
  data.frame(sc, gc, ld.cw, ld.ccw)
}
