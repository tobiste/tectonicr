#' @title Theoretical Direction of Maximum Horizontal Stress in the
#' geographical reference system.
#'
#' @description Models the direction of maximum horizontal stress
#' \eqn{\sigma_{Hmax}}{SHmax} along great circles, small circles, and
#' loxodromes at a given point or points according to the relative plate motion
#' in the geographical coordinate reference system.
#'
#' @author Tobias Stephan
#'
#' @param df \code{data.frame} containing the coordinates of the point(s)
#' (\code{lat}, \code{lon}).
#' @param euler \code{"data.frame"} or object of class \code{"euler.pole"}
#' containing the geographical coordinates of the Euler pole
#'
#' @details \eqn{\sigma_{Hmax}}{SHmax} following *great circles* is the
#' (initial) bearing between the given point and the pole of relative plate
#' motion. \eqn{\sigma_{Hmax}}{SHmax} along *small circles*, clockwise, and
#' counter-clockwise *loxodromes* is 90\eqn{^{\circ}}{ degree},
#' +45\eqn{^{\circ}}{ degree}, and 135\eqn{^{\circ}}{ degree}
#' (-45\eqn{^{\circ}}{ degree}) to this great circle bearing, respectively.
#'
#' @returns \code{data.frame}
#' \describe{
#'   \item{gc}{Azimuth of modeled \eqn{\sigma_{Hmax}}{SHmax} following
#'   great circles}
#'   \item{sc}{Small circles}
#'   \item{ld.cw}{Clockwise loxodromes}
#'   \item{ld.ccw}{Counter-clockwise loxodromes}
#'  }
#'
#' @seealso [deviation_shmax()] to compute the deviation of the modeled direction
#'  from the observed direction of \eqn{\sigma_{Hmax}}{SHmax}.
#'  [PoR_shmax()] to calculate the azimuth of \eqn{\sigma_{Hmax}}{SHmax}
#'  in the pole of rotation reference system.
#'
#' @references Stephan, T., Enkelmann, E., and Kroner, U. "Analyzing the
#' horizontal orientation of the crustal stress adjacent to plate boundaries".
#' *Sci Rep* 13. 15590 (2023). \doi{10.1038/s41598-023-42433-2}.
#'
#' @export
#'
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
  stopifnot(is.data.frame(df), is.data.frame(euler) | is.euler(euler))
  theta <- mapply(FUN = get_azimuth, lat_a = df$lat, lon_a = df$lon, lat_b = euler$lat, lon_b = euler$lon)
  # great circles
  gc <- theta %% 180

  # small circles
  sc <- (theta + 90) %% 180

  # counterclockwise loxodrome
  ld.ccw <- (theta + 135) %% 180

  # clockwise loxodrome
  ld.cw <- (theta + 45) %% 180

  data.frame(sc, ld.ccw, gc, ld.cw)
}


#' @title Normalize Angle Between Two Directions
#'
#' @description Normalizes the angle between two directions to the acute angle
#' in between, i.e. angles between 0 and 90\eqn{^\circ}{degree}
#'
#' @author Tobias Stephan
#'
#' @param x,y Minuend and subtrahend. Both numeric vectors of angles in degrees.
#' If `y` is missing, it treats `x` as difference. If not, length of subtrahend
#' `y` is either `1` or equal to `length(x)`.
#'
#' @returns numeric vector, acute angles between two directions, i.e. values
#' between 0 and 90\eqn{^\circ}{degree}
#'
#' @export
#'
#' @examples
#' deviation_norm(175, 5)
#' deviation_norm(c(175, 95, 0), c(5, 85, NA))
#' deviation_norm(c(-5, 85, 95, 175, 185, 265, 275, 355, 365))
deviation_norm <- function(x, y = NULL) {
  nx <- length(x)

  if (is.null(y)) {
    # # deviation is between 0 and 90
    # if (nx > 1) {
    #   for (i in seq_along(x)) {
    #     if (!is.na(x[i])) {
    #       while (abs(x[i]) > 90) {
    #         x[i] <- 180 - abs(x[i])
    #       }
    #     }
    #   }
    # } else {
    #   if (!is.na(x)) {
    #     while (abs(x) > 90) {
    #       x <- 180 - abs(x)
    #     }
    #   }
    # }
    # abs(x)
    y <- rep(0, nx)
  }

  ny <- length(y)
  stopifnot(nx == ny | ny == 1)
  if (ny == 1 & nx > 1) {
    ny <- rep(y, nx)
  }
  d <- (x %% 180) - (y %% 180)
  ifelse(d < 90, d, 180 - d)
}



#' Deviation of Observed and Predicted Directions of Maximum Horizontal Stress
#'
#' Calculate the angular difference between the observed and modeled direction
#' of maximum horizontal stresses (\eqn{\sigma_{Hmax}}{SHmax}) along
#' great circles, small circles, and
#' loxodromes of the relative plate motion's Euler pole
#'
#' @author Tobias Stephan
#'
#' @param prd \code{data.frame} containing the modeled azimuths of
#' \eqn{\sigma_{Hmax}}{SHmax}, i.e.
#' the return object from \code{model_shmax()}
#' @param obs Numeric vector containing the observed azimuth of
#' \eqn{\sigma_{Hmax}}{SHmax},
#' same length as \code{prd}
#' @returns An object of class \code{data.frame}
#'
#' \describe{
#'   \item{dev.gc}{Deviation of observed stress from modeled
#'   \eqn{\sigma_{Hmax}}{SHmax} following
#'   great circles}
#'   \item{dev.sc}{Small circles}
#'   \item{dev.ld.cw}{Clockwise loxodromes}
#'   \item{dev.ld.ccw}{Counter-clockwise loxodromes}
#' }
#'
#' @seealso [model_shmax()] to calculate the theoretical direction of
#' \eqn{\sigma_{Hmax}}{SHmax}.
#'
#' @references Stephan, T., Enkelmann, E., and Kroner, U. "Analyzing the
#' horizontal orientation of the crustal stress adjacent to plate boundaries".
#' *Sci Rep* 13. 15590 (2023). \doi{10.1038/s41598-023-42433-2}.
#'
#' @export
#'
#' @examples
#' data("nuvel1")
#' # North America relative to Pacific plate:
#' PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
#'
#' # the point where we want to model the SHmax direction:
#' point <- data.frame(lat = 45, lon = 20)
#'
#' prd <- model_shmax(point, PoR)
#' deviation_shmax(prd, obs = 90)
deviation_shmax <- function(prd, obs) {
  stopifnot(length(obs) == length(seq_along(prd$gc)))

  # normalize azimuth
  obs <- (obs + 360) %% 180

  dev.gc <- prd$gc - obs
  dev.sc <- prd$sc - obs
  dev.ld.cw <- prd$ld.cw - obs
  dev.ld.ccw <- prd$ld.ccw - obs

  data.frame(dev.gc, dev.sc, dev.ld.cw, dev.ld.ccw)
}


#' @title Direction of Maximum Horizontal Stress in PoR reference
#' system
#'
#' @description Models the direction of maximum horizontal stress
#' \eqn{\sigma_{Hmax}}{SHmax} in the Euler pole (Pole of Rotation)
#' coordinate reference system. When type of plate boundary is given, it also
#' gives the deviation from the theoretically predicted azimuth of
#' \eqn{\sigma_{Hmax}}{SHmax}, the deviation, and the normalized
#' \eqn{\chi^2}{chi-squared} statistics.
#'
#' @param df \code{data.frame} containing the coordinates of the point(s)
#' (\code{lat}, \code{lon}), the direction of
#' \eqn{\sigma_{Hmax}}{SHmax} \code{azi} and its standard deviation
#' \code{unc} (optional)
#' @param PoR \code{"data.frame"} or object of class \code{"euler.pole"}
#' containing the geographical coordinates of the Euler  pole
#' @param type Character. Type of plate boundary (optional). Can be
#' \code{"out"}, \code{"in"}, \code{"right"}, or
#' \code{"left"} for outward, inward, right-lateral, or left-lateral
#' moving plate boundaries, respectively. If \code{"none"} (the default), only
#' the PoR-equivalent azimuth is returned.
#'
#' @returns Either a numeric vector of the azimuths in the transformed coordinate
#' system (in degrees), or a \code{"data.frame"} with
#' \describe{
#' \item{`azi.PoR`}{the transformed azimuths (in degrees),}
#' \item{`prd`}{the predicted azimuths (in degrees),}
#' \item{`dev`}{the deviation between the transformed and the predicted azimuth (in degrees),}
#' \item{`nchisq`}{the Norm \eqn{\chi^2}{chi-squared} test statistic, and}
#' \item{`cdist`}{the angular distance between the transformed and the predicted azimuth.}
#' }
#'
#' @seealso [model_shmax()] to compute the theoretical direction of
#' \eqn{\sigma_{Hmax}}{SHmax} in the geographical reference system.
#' [deviation_shmax()] to compute the deviation of the modeled direction
#'  from the observed direction of \eqn{\sigma_{Hmax}}{SHmax}.
#'  [norm_chisq()] to calculate the normalized \eqn{\chi^2}{chi-squared}
#'  statistics. [circular_distance()] to calculate the angular distance.
#'
#' @details The azimuth of \eqn{\sigma_{Hmax}}{SHmax} in the pole of rotation
#' reference system is
#' approximate 0 (or 180), 45, 90, 135 degrees if the stress is sourced by an
#' outward, sinistral, inward, or dextral moving plate boundary, respectively.
#' directions of \eqn{\sigma_{Hmax}}{SHmax} with respect to the four
#' plate boundary types.
#'
#' @references Stephan, T., Enkelmann, E., and Kroner, U. "Analyzing the
#' horizontal orientation of the crustal stress adjacent to plate boundaries".
#' *Sci Rep* 13. 15590 (2023). \doi{10.1038/s41598-023-42433-2}.
#'
#' @export
#'
#' @examples
#' data("nuvel1")
#' # North America relative to Pacific plate:
#' PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
#'
#' data("san_andreas")
#' res <- PoR_shmax(san_andreas, PoR, type = "right")
#' head(res)
PoR_shmax <- function(df, PoR, type = c("none", "in", "out", "right", "left")) {
  stopifnot(is.data.frame(df), is.data.frame(PoR) | is.euler(PoR))
  type <- match.arg(type)

  theta <- mapply(FUN = get_azimuth, lat_a = df$lat, lon_a = df$lon, lat_b = PoR$lat, lon_b = PoR$lon)
  azi.por <- (df$azi - theta + 180) %% 180

  if (type != "none" && !is.null(df$unc)) {
    prd <- switch(type,
      "none" = NA,
      "out" = 180,
      "right" = 135,
      "in" = 90,
      "left" = 45
    )

    dev <- azi.por - prd
    cdist <- (1 - cosd(2 * dev)) / 2
    nchisq.i <- (deviation_norm(azi.por, prd) / 90)^2

    data.frame(
      "azi.PoR" = azi.por, "prd" = prd,
      "dev" = dev, "nchisq" = nchisq.i, "cdist" = cdist
    )
  } else {
    azi.por
  }
}

#' Azimuth conversion from PoR to geographical coordinate reference system
#'
#' Conversion of PoR azimuths into geographical azimuths
#'
#' @param x \code{data.frame} containing the PoR equivalent azimuths
#' (\code{azi.PoR}), and either the geographical coordinates of the
#' point(s) or the PoR-equivalent coordinates.
#' @param PoR \code{data.frame} containing the geographical location of
#' the Euler pole (\code{lat}, \code{lon})
#'
#' @seealso [PoR_shmax()]
#'
#' @returns numeric vector of transformed azimuths (in degrees)
#'
#' @references Stephan, T., Enkelmann, E., and Kroner, U. "Analyzing the
#' horizontal orientation of the crustal stress adjacent to plate boundaries".
#' *Sci Rep* 13. 15590 (2023). \doi{10.1038/s41598-023-42433-2}.
#'
#' @export
#'
#' @examples
#' data("nuvel1")
#' # North America relative to Pacific plate:
#' PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
#' data("san_andreas")
#' head(san_andreas$azi)
#' san_andreas$azi.PoR <- PoR_shmax(san_andreas, PoR)
#' res.geo <- PoR2Geo_azimuth(san_andreas, PoR)
#' head(res.geo)
PoR2Geo_azimuth <- function(x, PoR) {
  # Northern Hemisphere Euler pole
  if (PoR$lat < 0) {
    PoR$lat <- -PoR$lat
    PoR$lon <- longitude_modulo(180 + PoR$lon)
  }

  if (unique(c("lat.PoR", "lon.PoR") %in% colnames(x))) {
    northpole <- geographical_to_PoR(
      data.frame(lat = 90, lon = 0),
      PoR
    )
    beta <- mapply(FUN = get_azimuth, lat_a = x$lat.PoR, lon_a = x$lon.PoR, lat_b = northpole$lat.PoR, lon_b = northpole$lon.PoR)
    azi.geo <- x$azi.PoR - beta
  } else {
    beta <- mapply(FUN = get_azimuth, lat_a = x$lat, lon_a = x$lon, lat_b = PoR$lat, lon_b = PoR$lon)
    azi.geo <- x$azi.PoR + beta
  }
  azi.geo %% 180
}

#' SHmax direction resulting from multiple plate boundaries
#'
#' Calculates a \eqn{\sigma_{Hmax}}{SHmax} direction at given coordinates,
#' sourced by multiple plate boundaries. This first-order approximation is the
#' circular mean of the superimposed theoretical directions, weighted by the
#' rotation rates of the underlying PoRs.
#'
#' @param df `data.frame` containing the coordinates of the point(s)
#' (`lat`, `lon`), and the direction of
#' \eqn{\sigma_{Hmax}}{SHmax} `azi` (in degrees)
#' @param PoRs multirow `data.frame` or `"euler.pole"` object that must contain `lat`,
#' `lon` and `angle`
#' @param types character vector with length equal to number of rows in `PoRs`.
#' Type of plate boundary. Must be `"out"`, `"in"`, `"right"`, or
#' `"left"` for outward, inward, right-lateral, or left-lateral
#' moving plate boundaries, respectively.
#' @param absolute logical. Whether the resultant azimuth should be weighted
#' using the absolute rotation at the points or the angular rotation of the PoRs.
#' @param PoR_weighting (optional) numeric vector with length equal to number of rows in
#' `PoRs`. Extra weightings for the used PoRs.
#'
#' @seealso [model_shmax()]
#'
#' @return numeric. Resultant azimuth in degrees and geographical CRS
#' @export
#'
#' @examples
#' data(san_andreas)
#' data(nuvel1)
#' pors <- subset(nuvel1, plate.rot %in% c("eu", "na"))
#' superimposed_shmax(san_andreas, pors, types = c("in", "right"), PoR_weighting = c(2, 1))
superimposed_shmax <- function(df, PoRs, types, absolute = TRUE, PoR_weighting = NULL) {
  res <- c()
  lats <- c()
  stopifnot(is.character(types))
  if (is.null(PoR_weighting)) {
    PoR_weighting <- rep(1, nrow(PoRs))
  }

  if (!absolute) {
    lat_j <- rep(1, nrow(df))
    col <- 1
    while (col <= nrow(PoRs)) {
      lats <- cbind(lats, lat_j)
      col <- col + 1
    }
  }

  for (i in seq_along(PoRs$lat)) {
    res_i <- model_shmax(df, PoRs[i, ])
    if (absolute) lat_i <- PoR_coordinates(df, PoRs[i, ])$lat.PoR
    if (types[i] == "in") {
      azi <- res_i$sc
    } else if (types[i] == "out") {
      azi <- res_i$gc
    } else if (types[i] == "right") {
      azi <- res_i$ld.cw
    } else if (types[i] == "left") {
      azi <- res_i$ld.ccw
    } else {
      stop("`types` must be one of in, out, right, and left.")
    }

    res <- cbind(res, azi)
    if (absolute) lats <- cbind(lats, lat_i)
  }

  rot <- PoR_weighting * PoRs$angle * cosd(lats)
  azi_res <- numeric(length(res[, 1]))

  for (j in seq_along(res[, 1])) {
    azi_res[j] <- circular_mean(res[j, ], w = rot[j, ])
  }
  return(azi_res)
}
