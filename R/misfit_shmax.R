#' @title Normalize angle between two directions
#'
#' @description Normalizes the angle between two directions to the acute angle
#' in between, i.e. angles between 0 and 90 degrees.
#'
#' @author Tobias Stephan
#' @param x numeric vector containing angles in degrees
#' @return numeric vector, acute angles between two directions, i.e. values
#' between 0° and 90°
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
  return(abs(x))
}


#' Deviation of observed and predicted directions of SHmax
#'
#' Calculate the angular difference between the observed and modeled direction
#' of maximum horizontal stresses along great circles, small circles, and
#' loxodromes of the relative plate motion´s Euler pole
#'
#' @author Tobias Stephan
#' @param prd \code{data.frame} containing the modeled azimuths of SHmax, i.e.
#' the return object from \code{model_shmax()}
#' @param obs numeric vector containing the observed azimuth of SHmax,
#' same length as \code{prd}
#' @return An object of class \code{data.frame}
#' \describe{
#'   \item{dev.gc}{deviation of observed stress from modeled stress following
#'   great circles}
#'   \item{dev.sc}{small circles}
#'   \item{dev.ld.cw}{clockwise loxodromes}
#'   \item{dev.ld.ccw}{counter-clockwise loxodromes}
#'   }
#' @export
#' @examples
#' data("nuvel1")
#' euler <- subset(nuvel1, nuvel1$ID == "na") # North America relative to
#' # Pacific plate
#' point <- data.frame(lat = 45, lon = 20)
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

  dev.gc <- ifelse(dev.gc <= obs, -dev.gc, dev.gc)
  dev.sc <- ifelse(dev.sc <= obs, -dev.sc, dev.sc)
  dev.ld.cw <- ifelse(dev.ld.cw <= obs, -dev.ld.cw, dev.ld.cw)
  dev.ld.ccw <- ifelse(dev.ld.ccw <= obs, -dev.ld.ccw, dev.ld.ccw)


  dev.gc <- deviation_norm(dev.gc)
  dev.sc <- deviation_norm(dev.sc)
  dev.ld.cw <- deviation_norm(dev.ld.cw)
  dev.ld.ccw <- deviation_norm(dev.ld.ccw)

  dev.df <- data.frame(dev.gc, dev.sc, dev.ld.cw, dev.ld.ccw)
  return(dev.df)
}

