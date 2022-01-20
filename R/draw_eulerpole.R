#' @title Small circle grid dummy
#'
#' @description create a dummy
#'
#' @param x angle between small circles
#' @return data.frame
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
#' @description Construct a line that has a constant distance to the Euler pole.
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
eulerpole_smallcircles <- function(x, sm, small_circle = NULL) {
  if (missing(sm)) {
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
#' @description Construct curves cutting small circles perpendicular.
#'
#' @param x number of great circles
#' @return data.frame
#' @importFrom dplyr "%>%" mutate
#' @importFrom  geosphere greatCircleBearing
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
#' @description .
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




#' @title Rotate lines
#' @description Rotate a set of lines around a angle
#'
#' @param theta  Angle of rotation (in degree)
#' @param p Points of lines
#' @param centre Center point of rotation
#' @return \code{'matrix'}
#' @importFrom pracma cosd sind
rotate_lines <- function(theta, p, centre) {
  new_x <- pracma::cosd(theta) * (p[, 1] - centre[1]) - pracma::sind(theta) * (p[, 2] - centre[2]) + centre[1]
  new_y <- pracma::sind(theta) * (p[, 1] - centre[1]) + pracma::cosd(theta) * (p[, 2] - centre[2]) + centre[2]
  return(matrix(c(new_x, new_y), ncol = 2))
}


#' @title Loxodrome dummy
#' @description create a dummy for loxodromes
#'
#' @param angle Direction of loxodromes (in degree)
#' @param n number of loxodromes  (in  degree)
#' @param sense Sense of loxodromes 'sinistral' or 'dextral' for 'clockwise' or 'counterclockwise' loxodromes, respectively
#' @return data.frame
#' @importFrom dplyr "%>%" filter mutate
loxodrome_dummy <- function(angle = 45, n = 10, sense, lon = NULL, lat = NULL) {
  if (sense == "sinistral" | sense == "clockwise") {
    s <- -1
  } else if (sense == "dextral" | sense == "counterclockwise") {
    s <- 1
  } else {
    stop("sense must be sinistral, dextral, clockwise, or clounterclockwise")
  }

  lats <- seq(-180, 180, 1)

  line.dummy <- data.frame(
    lon = rep(0, length(lats)),
    lat = lats,
    line = 0
  )
  for (j in seq_along(line.dummy$lon)) {
    line.dummy.rot <- rotate_lines(theta = s * 45, p = cbind(line.dummy$lon[j], line.dummy$lat[j]), centre = c(line.dummy$lon[j], 0))
    loxodrome.dummy.j <- data.frame(lon = line.dummy.rot[, 1], lat = line.dummy.rot[, 2], loxodrome = line.dummy$line[j])

    if (j == 1) {
      loxodrome.dummy <- loxodrome.dummy.j
    } else {
      loxodrome.dummy <- rbind(loxodrome.dummy, loxodrome.dummy.j)
    }
  }

  for (i in seq(-360, 360, n)) {
    line.i <- loxodrome.dummy %>%
      mutate(
        lon = lon - i,
        loxodrome = i
      )

    if (i == -360) {
      loxodromes <- line.i
    } else {
      loxodromes <- rbind(loxodromes, line.i)
    }
  }

  # summary(loxodromes)
  #
  loxodromes.filt <- loxodromes %>%
    unique() %>%
    filter(abs(lat) <= 90) %>%
    filter(abs(lon) <= 180)

  return(loxodromes.filt)
}


#' @title Draws loxodromes directed towards an Euler pole
#'
#' @description Construct curves cutting small circles at a constant angle
#'
#' @author Tobias Stephan
#' @param x \code{data.frame} containing coordinates of Euler pole in lat and lon
#' @param angle Direction of loxodromes; default = 45
#' @param ld numeric, angle between loxodromes (in degree); default: 10
#' @param sense Sense of loxodromes  'sinistral' or 'dextral' for 'clockwise' or 'counterclockwise' loxodromes, respectively
#' @return An object of class \code{SpatialLinesDataFrame}
#' @importFrom dplyr first filter
#' @importFrom sp Line Lines SpatialLines SpatialLinesDataFrame CRS
#' @import rgdal
#' @export
#' @examples
#' data("nuvel1")
#' euler <- subset(nuvel1, nuvel1$ID == "na") # North America relative to Pacific plate
#' eulerpole_loxodromes(x = euler, angle = 45, ld = 10, sense = "sinistral")
eulerpole_loxodromes <- function(x, angle = 45, ld = 10, sense, loxodrome = NULL) {
  if (missing(sense)) {
    stop("sense missing")
  }

  ld.df <- loxodrome_dummy(angle = abs(angle), n = 360 / ld, sense = sense)


  ld_range <- unique(ld.df$loxodrome)
  SL.list <- list()

  # idee:
  # 1. transformiere ld.df data.frame in spatiallinesdataframe
  # 2. definiere koordinatesystem: oblique transmercator mit euler pole als fake north pole
  # 3. transformiere in wgs84
  for (l in unique(ld.df$loxodrome)) {
    # loop through all great circles
    ld.subset <- dplyr::filter(ld.df, loxodrome == l)

    l.i <- suppressWarnings(
      sp::Lines(
        slinelist = sp::Line(cbind(
          "lon" = ld.subset$lon,
          "lat" = ld.subset$lat
        )), ID = as.character(l)
      )
    )

    suppressWarnings(
      SL.list[as.character(l)] <- l.i
    )

    if (l == dplyr::first(unique(ld.df$loxodrome))) {
      ld <- ld.subset
    } else {
      ld <- rbind(ld, ld.subset)
    }
  }

  SL.wgs84 <- sp::SpatialLines(
    SL.list,
    proj4string =  sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  )


  SL.wgs84.df <- suppressWarnings(
    sp::SpatialLinesDataFrame(
      SL.wgs84,
      data.frame("loxodrome" = as.character(ld_range), row.names = ld_range)
    )
  )

  # transform into mercator
  SL.merc.df <- sp::spTransform(SL.wgs84.df, sp::CRS("+proj=merc"))

  # define lines in oblique mercator coordinate system with Euler pole as 'North pole'
  merc.obl <- sp::CRS(paste0("+proj=omerc +alpha=", x$lat, " +lon_c=", x$lon))
  suppressWarnings(
    sp::proj4string(SL.merc.df) <- merc.obl
  )

  # wrap at dateline
  SL.merc.df <- wrap_dateline(SL.merc.df)

  # transform back to WGS84
  SL.df <- sp::spTransform(SL.merc.df, sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

  # wrap at dateline
  SL.df <- wrap_dateline(SL.df)

  return(SL.df)
}
