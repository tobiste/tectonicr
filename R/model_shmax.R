#'  Modeling directions of maximum horizontal stress
#'
#' Calculate the direction of maximum horizontal stress along great circles,
#' small circles, and loxodromes around the relative plate motion´s Euler pole
#' at a given point or points
#'
#' @author Copyright (C) 2021 Tobias Stephan
#' @param df data.frame containing the coordinates of the points (lat, lon)
#' @param euler data.frame containing the coordinates of the Euler pole for the
#' responsible plate boundary (lat, lon)
#' @references Wdowinski, S., 1998, A theory of intraplate
#'   tectonics: Journal of Geophysical Research: Solid Earth, v. 103, p.
#'   5037-5059, http://dx.doi.org/10.1029/97JB03390.
#' @return data.frame containing the azimuth of the modeled maximum horizontal
#' stress great circles (gc), small circles (sc), and loxodromes (ld.cw, ld.ccw)
#' at the given point
#' @export
#' @examples
#' point <- data.frame(lat = 45, lon = 20)
#' euler <- data.frame(lat = 90, lon = 0)
#' plate_vector(point, euler)
model_shmax <- function(df, euler) {


  azi <- df$azi %% 180

  for(i in seq_along(df$lat)){
    # great circles
    gc <- get_azimuth(
      c(df$lat[i], df$lon[i]), c(euler$lat[1], euler$lon[1])
    ) %% 180

    # plate motion vector is perpendicular to the bearing between the point and
    # the Euler pole
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
