#' @title Azimuth Between two Points
#'
#' @description Calculate initial bearing (or forward azimuth/direction) to go
#' from point \code{a} to point \code{b} following great circle arc on a
#' sphere.
#'
#' @param a,b Coordinates of start and end point(s).
#' Can be vectors of two numbers or a matrix of 2 columns (latitude, longitude)
#' @details This formula is for the initial bearing (sometimes referred to as
#' forward azimuth) which if followed in a straight line along a great circle
#' arc will lead from the start point a to the end point b.
#' \deqn{\theta = \arctan2 (\sin \Delta\lambda *
#' \cos\psi_2, \cos\psi_1 \sin\psi_1-\sin\psi_1 \cos\psi_2 \cos\Delta\lambda)}
#' where  \eqn{\psi_1, \lambda_1} is the start point, \eqn{\psi_2},
#' \eqn{\lambda_2} the end point (\eqn{\Delta\lambda} is the difference in
#' longitude)
#' @source \url{http://www.movable-type.co.uk/scripts/latlong.html}
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


#' @title Theoretical Direction of Maximum Horizontal Stress in the
#' geographical reference system.
#'
#' @description Models the direction of maximum horizontal stress
#' \eqn{\sigma_\text{Hmax}}{SHmax} along great circles, small circles, and
#' loxodromes at a given point or points according to the relative plate motion
#' in the geographical coordinate reference system.
#'
#' @author Tobias Stephan
#' @param df \code{data.frame} containing the coordinates of the point(s)
#' (\code{lat}, \code{lon}).
#' @param euler \code{data.frame} containing the coordinates of the Euler pole
#' for the plate boundary (\code{lat}, \code{lon}).
#' @details \eqn{\sigma_\text{Hmax}}{SHmax} following *great circles* is the
#' (initial) bearing between the given point and the pole of relative plate
#' motion. \eqn{\sigma_\text{Hmax}}{SHmax} along *small circles*, clockwise, and
#' counter-clockwise *loxodromes* is 90\eqn{^{\circ}}{ degree},
#' +45\eqn{^{\circ}}{ degree}, and 135\eqn{^{\circ}}{ degree}
#' (-45\eqn{^{\circ}}{ degree}) to this great circle bearing, respectively.
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
#'  [PoR_shmax()] to calculate the azimuth of \eqn{\sigma_\text{Hmax}}{SHmax}
#'  in the pole of rotation reference system.
#' @export
#' @examples
#' data("nuvel1")
#' # North America relative to Pacific plate:
#' euler <- subset(nuvel1, nuvel1$plate.rot == "na")
#'
#' # the point where we mant to model the SHmax direction:
#' point <- data.frame(lat = 45, lon = 20)
#'
#' model_shmax(point, euler)
model_shmax <- function(df, euler) {
  stopifnot(is.data.frame(df) & is.data.frame(euler))

  beta <- c()
  sc <- c()
  gc <- c()
  ld.cw <- c()
  ld.ccw <- c()
  for (i in seq_along(df$lat)) {
    beta[i] <- get_azimuth(
      c(df$lat[i], df$lon[i]),
      c(euler$lat[1], euler$lon[1])
    )

    # great circles
    gc[i] <- (beta[i] + 180) %% 180

    # small circles
    sc[i] <- (beta[i] + 90 + 180) %% 180

    # counterclockwise loxodrome
    ld.ccw[i] <- (beta[i] + 135 + 180) %% 180

    # clockwise loxodrome
    ld.cw[i] <- (beta[i] + 45 + 180) %% 180
  }
  data.frame(sc, ld.ccw, gc, ld.cw)
}

#' @title Normalize Angle Between Two Directions
#'
#' @description Normalizes the angle between two directions to the acute angle
#' in between, i.e. angles between 0 and 90°
#'
#' @author Tobias Stephan
#' @param x Numeric vector containing angles in degrees
#' @return Numeric vector, acute angles between two directions, i.e. values
#' between 0 and 90°
#' @export
#' @examples
#'
#' deviation_norm(91)
#' deviation_norm(c(-91, NA, 23497349))
deviation_norm <- function(x) {
  # deviation is between 0 and 90
  if (length(x) > 1) {
    for (i in seq_along(x)) {
      if (!is.na(x[i])) {
        while (abs(x[i]) > 90) {
          x[i] <- 180 - abs(x[i])
        }
      }
    }
  } else {
    if (!is.na(x)) {
      while (abs(x) > 90) {
        x <- 180 - abs(x)
      }
    }
  }
  abs(x)
}



#' Deviation of Observed and Predicted Directions of Maximum Horizontal Stress
#'
#' Calculate the angular difference between the observed and modeled direction
#' of maximum horizontal stresses (\eqn{\sigma_\text{Hmax}}{SHmax}) along
#' great circles, small circles, and
#' loxodromes of the relative plate motion's Euler pole
#'
#' @author Tobias Stephan
#' @param prd \code{data.frame} containing the modeled azimuths of
#' \eqn{\sigma_\text{Hmax}}{SHmax}, i.e.
#' the return object from \code{model_shmax()}
#' @param obs Numeric vector containing the observed azimuth of
#' \eqn{\sigma_\text{Hmax}}{SHmax},
#' same length as \code{prd}
#' @return An object of class \code{data.frame}
#' \describe{
#'   \item{dev.gc}{Deviation of observed stress from modeled
#'   \eqn{\sigma_\text{Hmax}}{SHmax} following
#'   great circles}
#'   \item{dev.sc}{Small circles}
#'   \item{dev.ld.cw}{Clockwise loxodromes}
#'   \item{dev.ld.ccw}{Counter-clockwise loxodromes}
#' }
#' @seealso [model_shmax()] to calculate the theoretical direction of
#' \eqn{\sigma_\text{Hmax}}{SHmax}.
#' @export
#' @examples
#' data("nuvel1")
#' # North America relative to Pacific plate:
#' euler <- subset(nuvel1, nuvel1$plate.rot == "na")
#'
#' # the point where we want to model the SHmax direction:
#' point <- data.frame(lat = 45, lon = 20)
#'
#' prd <- model_shmax(point, euler)
#' misfit_shmax(prd, obs = 90)
misfit_shmax <- function(prd, obs) {
  stopifnot(length(obs) == length(seq_along(prd$gc)))

  # normalize azimuth
  obs <- (obs + 180) %% 180

  dev.gc <- prd$gc - obs
  dev.sc <- prd$sc - obs
  dev.ld.cw <- prd$ld.cw - obs
  dev.ld.ccw <- prd$ld.ccw - obs

  data.frame(dev.gc, dev.sc, dev.ld.cw, dev.ld.ccw)
}


#' @title Theoretical Direction of Maximum Horizontal Stress in PoR reference system
#'
#' @description Models the direction of maximum horizontal stress
#' \eqn{\sigma_\text{Hmax}}{SHmax} in the Euler pole (Pole of Rotation)
#' coordinate reference system. When type of plate boundary is given, it also
#' gives the deviation from the theoretically predicted azimuth of
#' \eqn{\sigma_\text{Hmax}}{SHmax}, the deviation, and the normalized
#' \eqn{\chi^2}{chi-squared} statistics.
#' @param df \code{data.frame} containing the coordinates of the point(s)
#' (\code{lat}, \code{lon}), the orientation of
#' \eqn{\sigma_\text{Hmax}}{SHmax} \code{azi} and its standard deviation
#' \code{unc} (optional)
#' @param euler \code{data.frame} containing the coordinates of the Euler pole
#' for the plate boundary  (\code{lat}, \code{lon})
#' @param type Character. Type of plate boundary (optional). Can be
#' \code{"out"}, \code{"in"}, \code{"right"}, or
#' \code{"left"} for outward, inward, right-lateral, or left-lateral
#' moving plate boundaries, respectively.
#' @return Either a numeric vector of the azimuths in the transformed coordinate
#' system, or a \code{"data.frame"} with
#' \describe{
#' \item{\code{"azi.PoR"}}{the transformed azimuths,}
#' \item{\code{"prd"}}{the predicted azimuths,}
#' \item{\code{"dev"}}{the deviation, and}
#' \item{\code{"nchisq"}}{the normalized \eqn{\chi^2}{chi-squared} statistics.}
#' }
#' @seealso [model_shmax()] to compute the theoretical orientation of
#' \eqn{\sigma_\text{Hmax}}{SHmax} in the geographical reference system.
#' [misfit_shmax()] to compute the deviation of the modeled direction
#'  from the observed direction of \eqn{\sigma_\text{Hmax}}{SHmax}.
#'  [norm_chisq()] to calculate the normalized \eqn{\chi^2}{chi-squared}
#'  statistics.
#' @details According to the theory, the azimuth of
#' \eqn{\sigma_\text{Hmax}}{SHmax} in the pole of rotation reference system is
#' approximate 0 (or 180), 45, 90, 135 degrees if the stress is sourced by an
#' outward, sinistral, inward, or dextral moving plate boundary, respectively.
#' directions of \eqn{\sigma_\text{Hmax}}{SHmax} with respect to the four
#' plate boundary types.
#' @export
#' @examples
#' data("nuvel1")
#' # North America relative to Pacific plate:
#' euler <- subset(nuvel1, nuvel1$plate.rot == "na")
#'
#' data("san_andreas")
#' res <- PoR_shmax(san_andreas, euler, type = "right")
#' head(res)
PoR_shmax <- function(df, euler, type = c("in", "out", "right", "left")) {
  stopifnot(is.data.frame(df) & is.data.frame(euler))
  type <- match.arg(type)
  beta <- c()
  for (i in seq_along(df$lat)) {
    beta[i] <- (get_azimuth(
      c(df$lat[i], df$lon[i]),
      c(euler$lat[1], euler$lon[1])
    ) + 180) %% 180
  }
  azi.por <- (df$azi - beta + 180) %% 180

  if (!missing(type) & !is.null(df$unc)) {
    prd <- NA
    prd <- ifelse(type == "out", 180, prd)
    prd <- ifelse(type == "right", 135, prd)
    prd <- ifelse(type == "in", 90, prd)
    prd <- ifelse(type == "left", 45, prd)

    dev <- azi.por - prd
    nchisq.i <- (deviation_norm(dev) / df$unc)^2 / (90 / df$unc)^2

    data.frame("azi.PoR" = azi.por, "prd" = prd, "dev" = dev, "nchisq" = nchisq.i)
  } else {
    azi.por
  }
}
