#'  Modeling directions of maximum horizontal stress
#'
#' Calculate the direction of maximum horizontal stress along great circles,
#' small circles, and loxodromes around the relative plate motion´s Euler pole
#' at a given point or points
#'
#' @author Copyright (C) 2021 Tobias Stephan
#' @param df \code{"data.frame"} containing the coordinates of the points (lat,
#' lon)
#' @param euler \code{"data.frame"} containing the coordinates of the Euler pole
#' for the plate boundary (lat, lon)
#' @references Wdowinski, S., 1998, A theory of intraplate
#'   tectonics: Journal of Geophysical Research: Solid Earth, v. 103, p.
#'   5037-5059, http://dx.doi.org/10.1029/97JB03390.
#' @return An object of class \code{"data.frame"}
#' \describe{
#'   \item{gc} {azimuth of the modeled maximum horizontal following stress
#'   great circles}
#'   \item{sc}{small circles}
#'   \item{ld.cw}{clockwise loxodromes}
#'   \item{ld.ccw}{counter-clockwise loxodromes}
#'   }
#' @export
#' @examples
#' data("nuvel1")
#' euler <- subset(nuvel1, nuvel1$ID=='na') #North America relative to Pacific plate
#' point <- data.frame(lat = 45, lon = 20)
#' model_shmax(point, euler)
#'
model_shmax <- function(df, euler) {


  azi <- df$azi %% 180

  for(i in seq_along(df$lat)){
    # great circles
    gc <- get_azimuth(
      c(df$lat[i], df$lon[i]), c(euler$lat[1], euler$lon[1])
    ) %% 180

    # plate vector is perpendicular to the bearing between the point and the
    # Euler pole
    sc <- (
      get_azimuth(c(df$lat[i], df$lon[i]), c(euler$lat[1], euler$lon[1])) - 90
    ) %% 180


    # counterclockwise loxodrome
    ld.ccw <- (
      get_azimuth(c(df$lat[i], df$lon[i]), c(euler$lat[1], euler$lon[1])) - 45
    ) %% 180


    # clockwise loxodrome
    ld.cw <- (
      get_azimuth(c(df$lat[i], df$lon[i]), c(euler$lat[1], euler$lon[1])) + 45
    ) %% 180

    x.i <- data.frame(
      gc, sc,ld.ccw, ld.cw
    )
    if(i == 1){
      x  <- x.i
    } else {
      x <- rbind(x, x.i)
    }
  }
  return(x)
}
