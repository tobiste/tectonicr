#' @title Small circle grid dummy
#'
#' @description create a dummy
#'
#' @param x angle between small circles
#' @return data.frame
#' @export
#' @examples
#' smallcircle_dummy(30)
smallcircle_dummy <- function(x) {
  sm_range <- seq(0, 180, x)
  lons <- seq(-180, 180, x)

  sm.df <- data.frame(
    "lon" = as.numeric(),
    "lat" = as.numeric(),
    "small_circle" = as.numeric()
  )

  for (i in sm_range) {
    # loop through all small circles
    sm <- sm_range[sm_range == i]
    lat <- sm - 90
    sm.l <- data.frame(
      "lat" = rep(lat, length(lons)),
      "lon" = lons,
      "small_circle" = sm
    )
    sm.df <- rbind(sm.df, sm.l)
  }
  return(sm.df)
}

#' @title Wrap Spatial object at dateline
#'
#' @description Wrap a Spatial object at date line
#'
#' @param x Spatial object
#' @return Spatial object
#' @importFrom sf st_as_sf as_Spatial st_wrap_dateline
#' @export
wrap_dateline <- function(x) {
  x.wrapped.sf <- sf::st_wrap_dateline(
    sf::st_as_sf(x),
    options =  c("WRAPDATELINE=YES", "DATELINEOFFSET=180"),
    quiet = TRUE
  )
  x.wrapped <- sf::as_Spatial(x.wrapped.sf)

  return(x.wrapped)
}

#' @title Draws small-circles around Euler pole
#'
#' @description Calculate the direction of maximum horizontal stress along great circles,
#' small circles, and loxodromes around the relative plate motion´s Euler pole
#' at a given point or points
#'
#' @author Tobias Stephan
#' @param x \code{data.frame} containing coordinates of Euler pole in lat, lon, and rotation angle (optional)
#' @param sm numeric, angle between small circles (in degree); default: 10
#' @return An object of class \code{SpatialLinesDataFrame}
#' @details if angle is given: output additionally gives absolute velocity on small-circle (degree/Myr -> km/Myr)
#' @importFrom dplyr "%>%" mutate select
#' @importFrom pracma cross
#' @importFrom sp Line Lines SpatialLines SpatialLinesDataFrame CRS
#' @export
#' @examples
#' data("nuvel1")
#' euler <- subset(nuvel1, nuvel1$ID == "na") # North America relative to
#' euler$angle <- euler$rate
#' # Pacific plate
#' eulerpole_smallcircles(euler)
eulerpole_smallcircles <- function(x, sm) {
  if(missing(sm)){
    sm <- 10
    }
  sm.df <- smallcircle_dummy(sm)
  sm_range <- unique(sm.df$small_circle)

  if (is.null(x$angle)) {
    sm_range.df <- sm_range
  } else {
    velocity <- sm.df %>%
      dplyr::mutate(abs_vel = abs_vel(x$angle, sm.df$small_circle)) %>%
      dplyr::select(small_circle, abs_vel) %>%
      unique()
    sm_range.df <- velocity
  }

  SL.list <- list()

  for (s in unique(sm.df$small_circle)) {
    # loop through all small circles
    sm.subset <- subset(sm.df, sm.df$small_circle == s)
    sm.subset.rot <- sm.subset

    for (i in seq_along(sm.subset$lat)) {
      # loop through all coordinates

      pt <- c(sm.subset$lat[i], sm.subset$lon[i])

      # Rotation matrix: axis is axis perpendicular to north pole 0/90 and Euler
      # pole, rotation-angle is angle between Northpole and Euler pole
      rot.axis <- pracma::cross(
        geographical_to_cartesian(c(90, 0)),
        geographical_to_cartesian(c(x$lat, x$lon))
      )
      rot.angle <-
        angle_vectors(
          geographical_to_cartesian(c(90, 0)),
          geographical_to_cartesian(c(x$lat, x$lon))
        )
      if (is.nan(rot.angle)) {
        rot.angle <- 0
      } # if there is no angle between,
      # than angle=0

      rot.matrix <- rotation_matrix(rot.axis, rot.angle)

      pt1 <- cartesian_to_geographical(rot.matrix %*% geographical_to_cartesian(pt))

      sm.subset.rot$lat[i] <- pt1[1]
      sm.subset.rot$lon[i] <- pt1[2]
    }

    l.i <- sp::Lines(slinelist = sp::Line(cbind(
      "lon" = sm.subset.rot$lon,
      "lat" = sm.subset.rot$lat
    )), ID = as.character(s))
    SL.list[as.character(s)] <- l.i

    if (s == 0) {
      sm.rot <- sm.subset.rot
    } else {
      sm.rot <- rbind(sm.rot, sm.subset.rot)
    }
  }

  SL.t <- sp::SpatialLinesDataFrame(
    sp::SpatialLines(SL.list, proj4string = sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")),
    data.frame(
      "small_circle" = sm_range.df,
      row.names = sm_range
    )
  )
  SL.t <- wrap_dateline(SL.t)

  return(SL.t)
}

#' @title Great circle grid dummy
#'
#' @description create a dummy
#'
#' @param x number of great circles
#' @return data.frame
#' @importFrom dplyr "%>%" mutate
#' @importFrom  geosphere greatCircleBearing
#' @export
#' @examples
#' greatcircle_dummy(4)
greatcircle_dummy <- function(x) {
  angle <- 360 / x
  i <- 0
  while (i <= 360) {
    if (i %% 180 == 0) {
      great.circle.i <- data.frame(
        lon = rep(0, 181),
        lat = seq(-90, 90, 1),
        great.circle = i
      )
    } else {
      great.circle.i <-
        geosphere::greatCircleBearing(c(0, 0), i) %>%
        as.data.frame() %>%
        dplyr::mutate(great.circle = i)
    }

    if (i == 0) {
      great.circle <- great.circle.i
    } else {
      great.circle <- rbind(great.circle, great.circle.i)
    }
    i <- i + angle
  }
  return(great.circle)
}

#' @title Draws great-circles of Euler pole
#'
#' @description Calculate the direction of maximum horizontal stress along great circles,
#' small circles, and loxodromes around the relative plate motion´s Euler pole
#' at a given point or points
#'
#' @author Tobias Stephan
#' @param x \code{data.frame} containing coordinates of Euler poles in lat, lon, and rotation angle (optional)
#' @param gm numeric, angle between great circles
#' @param n numeric; number of equally spaced great circles (angle between great circles (number of great circles n = 360 / gm); default = 12
#' @return An object of class \code{SpatialLinesDataFrame}
#' @details if angle is given: output additionally gives absolute velocity on small-circle (degree/Myr -> km/Myr)
#' @importFrom dplyr "%>%" first
#' @importFrom pracma cross
#' @importFrom sp Line Lines SpatialLines SpatialLinesDataFrame CRS
#' @export
#' @examples
#' data("nuvel1")
#' euler <- subset(nuvel1, nuvel1$ID == "na") # North America relative to
#' # Pacific plate
#' euler$angle <- euler$rate
#' eulerpole_greatcircles(euler)
eulerpole_greatcircles <- function(x, gm, n) {
  if (missing(n) & !missing(gm)) {
    n <- round(360 / gm, 0)
  } else if (missing(n) & missing(gm)) {
    n <- 12
  } else if (!missing(n) & !missing(gm)) {
    warning("Both gm and n are given. Only n is considered")
    n <- 12
  }

  gm.df <- greatcircle_dummy(n)
  gm_range <- unique(gm.df$great.circle)
  SL.list <- list()

  for (g in unique(gm.df$great.circle)) {
    # loop through all great circles
    gm.subset <- subset(gm.df, gm.df$great.circle == g)
    gm.subset.rot <- gm.subset

    for (i in seq_along(gm.subset$lat)) {
      # loop through all coordinates
      pt <- c(gm.subset$lat[i], gm.subset$lon[i])

      # Rotation matrix: axis is axis perpendicular to gm origin pole 0/0 and Euler
      # pole, rotation-angle is angle between gm origin and euler pole
      rot.axis <- pracma::cross(
        geographical_to_cartesian(c(0, 0)),
        geographical_to_cartesian(c(x$lat, x$lon))
      )
      rot.angle <-
        angle_vectors(
          geographical_to_cartesian(c(0, 0)),
          geographical_to_cartesian(c(x$lat, x$lon))
        )
      if (is.nan(rot.angle)) {
        rot.angle <- 0
      } # if there is no angle between,
      # than angle=0

      rot.matrix <- rotation_matrix(rot.axis, rot.angle)

      pt1 <- cartesian_to_geographical(rot.matrix %*% geographical_to_cartesian(pt))

      gm.subset.rot$lat[i] <- pt1[1]
      gm.subset.rot$lon[i] <- pt1[2]
    }

    l.i <- suppressWarnings(
      sp::Lines(
        slinelist = sp::Line(cbind(
          "lon" = gm.subset.rot$lon,
          "lat" = gm.subset.rot$lat
        )), ID = as.character(g)
      )
    )

    SL.list[as.character(g)] <- l.i

    if (g == dplyr::first(unique(gm.df$great.circle))) {
      gm.rot <- gm.subset.rot
    } else {
      gm.rot <- rbind(gm.rot, gm.subset.rot)
    }
  }

  SL.t <- sp::SpatialLines(SL.list, proj4string = sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

  SL.t.df <- suppressWarnings(
    sp::SpatialLinesDataFrame(
      SL.t,
      data.frame("great_circle" = as.character(gm_range), row.names = gm_range)
    )
  )

  SL.t.df <- wrap_dateline(SL.t.df)

  return(SL.t.df)
}
