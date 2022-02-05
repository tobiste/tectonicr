#' @title Small circle grid dummy
#'
#' @description create a dummy
#'
#' @param x angle between small circles
#' @return data.frame
#' @export
#' @examples
#' smallcircle_dummy(10)
smallcircle_dummy <- function(x) {
  sm_range <- seq(0, 180, x)
  lons <- seq(-180, 180, x)

  sm.df <- data.frame("lon" = as.numeric(),
                      "lat" = as.numeric(),
                      "small_circle" = as.numeric())

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


#' @title Small-circles around Euler pole
#'
#' @description Construct a line that has a constant distance to the Euler pole.
#'
#' @author Tobias Stephan
#' @param x \code{data.frame} containing coordinates of Euler pole in lat, lon,
#' and rotation angle (optional)
#' @param sm numeric, angle between small circles (in degree); default: 10
#' @param sf logical. Export object type of the lines.
#' TRUE (the default) simple feature (\code{sf}) objects.
#' FALSE for \code{"SpatialLInesDataFrame"} object.
#' @return An object of class \code{sf} or \code{SpatialLinesDataFrame}
#' @details if angle is given: output additionally gives absolute velocity on
#' small circle (degree/Myr -> km/Myr)
#' @importFrom dplyr "%>%" mutate select
#' @importFrom sp Line Lines SpatialLines SpatialLinesDataFrame
#' @importFrom sf st_crs st_as_sf st_set_crs st_transform st_wrap_dateline as_Spatial
#' @export
#' @examples
#' data("nuvel1")
#' euler <- subset(nuvel1, nuvel1$ID == "na") # North America relative to
#' euler$angle <- euler$rate
#' # Pacific plate
#' eulerpole_smallcircles(euler)
eulerpole_smallcircles <- function(x, sm, sf = TRUE) {
  small_circle <- NULL
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

    l.i <- sp::Lines(slinelist = sp::Line(cbind(
      "lon" = sm.subset$lon,
      "lat" = sm.subset$lat
    )), ID = as.character(s))

    suppressWarnings(SL.list[as.character(s)] <- l.i)
  }

  wgs84 <- sf::st_crs("+proj=longlat +datum=WGS84")
  ep <- sf::st_crs(paste0(
    "+proj=ob_tran +o_proj=longlat +datum=WGS84 +o_lat_p=",
    x$lat,
    " +o_lon_p=",
    x$lon)
  )

  SL.wgs84 <- sp::SpatialLines(SL.list)
  SL.wgs84.df <- sp::SpatialLinesDataFrame(
    SL.wgs84,
    data.frame("small_circle" = sm_range.df, row.names = sm_range)
  )

  suppressMessages(
    suppressWarnings(
      SL <- SL.wgs84.df %>%
        sf::st_as_sf() %>%
        sf::st_set_crs(wgs84) %>%
        sf::st_transform(ep) %>%
        sf::st_set_crs(wgs84) %>%
        sf::st_wrap_dateline(options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))
    )
  )

  if (!sf) {
    SL <- sf::as_Spatial(SL)
  }
  return(SL)
}



#' @title Great-circles of Euler pole
#'
#' @description Get points on a great circle as defined by the shortest distance
#'  between two specified points
#'
#' @author Tobias Stephan
#' @param x \code{data.frame} containing coordinates of Euler poles in lat, lon,
#'  and rotation angle (optional)
#' @param gm numeric, angle between great circles
#' @param n numeric; number of equally spaced great circles (angle between great
#' circles (number of great circles n = 360 / gm); default = 12
#' @param sf logical. Export object type of the lines.
#' TRUE (the default) simple feature (\code{sf}) objects.
#' FALSE for \code{"SpatialLInesDataFrame"} object.
#' @return An object of class \code{sf} or \code{SpatialLinesDataFrame}
#' @export
#' @examples
#' data("nuvel1")
#' euler <- subset(nuvel1, nuvel1$ID == "na") # North America relative to Pacific plate
#' euler$angle <- euler$rate
#' eulerpole_greatcircles(euler)
eulerpole_greatcircles <- function(x, gm, n, sf=TRUE) {
  if (missing(n) & !missing(gm)) {
    n <- round(360 / gm, 0)
  } else if (missing(n) & missing(gm)) {
    n <- 12
  } else if (!missing(n) & !missing(gm)) {
    warning("Both gm and n are given. Only n is considered\n")
    n <- 12
  }

  SL <- eulerpole_loxodromes(x, angle = 0, ld = n, sense = "dextral", sf = TRUE)
  if (sf) {
    names(SL)[1] <- "greatcircle"
  } else {
    names(SL) <- "greatcircle"
  }
  return(SL)
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
  new_x <-
    pracma::cosd(theta) * (p[, 1] - centre[1]) - pracma::sind(theta) *
    (p[, 2] - centre[2]) + centre[1]
  new_y <-
    pracma::sind(theta) * (p[, 1] - centre[1]) + pracma::cosd(theta) *
    (p[, 2] - centre[2]) + centre[2]
  return(matrix(c(new_x, new_y), ncol = 2))
}


#' @title Loxodrome dummy
#' @description create a dummy for loxodromes
#'
#' @param angle Direction of loxodromes (in degree)
#' @param n number of loxodromes  (in  degree)
#' @param sense Sense of loxodromes 'sinistral' or 'dextral' for 'clockwise' or
#' 'counterclockwise' loxodromes, respectively
#' @return data.frame
#' @importFrom dplyr "%>%" filter mutate
#' @export
#' @examples
#' loxodrome_dummy(angle = 45, n = 10, sense = "sinistral")
loxodrome_dummy <- function(angle, n, sense) {
  lon <- lat <- NULL
  if (sense == "sinistral" | sense == "clockwise") {
    s <- -1
  } else if (sense == "dextral" | sense == "counterclockwise") {
    s <- 1
  } else {
    stop("sense must be sinistral, dextral, clockwise, or clounterclockwise\n")
  }

  lats <- seq(-180, 180, 1)

  line.dummy <- data.frame(lon = rep(0, length(lats)),
                           lat = lats,
                           line = 0)
  for (j in seq_along(line.dummy$lon)) {
    line.dummy.rot <-
      rotate_lines(
        theta = s * angle,
        p = cbind(line.dummy$lon[j],
                  line.dummy$lat[j]),
        centre = c(line.dummy$lon[j], 0)
      )
    loxodrome.dummy.j <- data.frame(lon = line.dummy.rot[, 1],
                                    lat = line.dummy.rot[, 2],
                                    loxodrome = line.dummy$line[j])

    if (j == 1) {
      loxodrome.dummy <- loxodrome.dummy.j
    } else {
      loxodrome.dummy <- rbind(loxodrome.dummy, loxodrome.dummy.j)
    }
  }

  for (i in seq(-360, 360, n)) {
    line.i <- loxodrome.dummy %>%
      mutate(lon = lon - i,
             loxodrome = i)

    if (i == -360) {
      loxodromes <- line.i
    } else {
      loxodromes <- rbind(loxodromes, line.i)
    }
  }

  loxodromes.filt <- loxodromes %>%
    unique() %>%
    filter(abs(lat) <= 90) %>%
    filter(abs(lon) <= 180)

  return(loxodromes.filt)
}


#' @title Loxodromes directed towards an Euler pole
#'
#' @description Construct curves cutting small circles at a constant angle
#'
#' @author Tobias Stephan
#' @param x \code{data.frame} containing coordinates of Euler pole in lat and
#' lon
#' @param angle Direction of loxodromes; default = 45
#' @param ld numeric, angle between loxodromes (in degree); default: 10
#' @param sense Sense of loxodromes  'sinistral' or 'dextral' for 'clockwise'
#' or 'counterclockwise' loxodromes, respectively
#' @param sf logical. Export object type of the lines.
#' TRUE (the default) simple feature (\code{sf}) objects.
#' FALSE for \code{"SpatialLInesDataFrame"} object.
#' @return An object of class \code{sf} or \code{SpatialLinesDataFrame}
#' @importFrom dplyr first filter
#' @importFrom sp Line Lines SpatialLines SpatialLinesDataFrame
#' @importFrom sf st_crs st_as_sf st_set_crs st_transform st_wrap_dateline as_Spatial
#' @import rgdal
#' @export
#' @examples
#' data("nuvel1")
#' euler <- subset(nuvel1, nuvel1$ID == "na") # North America relative to Pacific plate
#' eulerpole_loxodromes(x = euler, angle = 45, ld = 10, sense = "sinistral")
eulerpole_loxodromes <- function(x, angle = 45, ld = 10, sense, sf = TRUE) {
  loxodrome <- NULL
  if (missing(sense)) {
    stop("sense missing\n")
  }

  ld.df <-
    loxodrome_dummy(angle = abs(angle),
                    n = 360 / ld,
                    sense = sense)

  ld_range <- unique(ld.df$loxodrome)
  SL.list <- list()

  for (l in unique(ld.df$loxodrome)) {
    # loop through all circles
    ld.subset <- dplyr::filter(ld.df, loxodrome == l)

    l.i <- suppressWarnings(sp::Lines(slinelist = sp::Line(
      cbind("lon" = ld.subset$lon,
            "lat" = ld.subset$lat)
    ), ID = as.character(l)))

    suppressWarnings(SL.list[as.character(l)] <- l.i)
  }

  wgs84 <- sf::st_crs("+proj=longlat +datum=WGS84")
  ep <- sf::st_crs(paste0(
    "+proj=ob_tran +o_proj=longlat +datum=WGS84 +o_lat_p=",
    x$lat,
    " +o_lon_p=",
    x$lon)
  )

  SL.wgs84 <- sp::SpatialLines(SL.list)
  SL.wgs84.df <- sp::SpatialLinesDataFrame(
    SL.wgs84,
    data.frame("loxodrome" = as.character(ld_range), row.names = ld_range)
  )

  suppressMessages(
    suppressWarnings(
      SL <- SL.wgs84.df %>%
        sf::st_as_sf() %>%
        sf::st_set_crs(wgs84) %>%
        sf::st_transform(ep) %>%
        sf::st_set_crs(wgs84) %>%
        sf::st_wrap_dateline(options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))
    )
  )

  if (!sf) {
    SL <- sf::as_Spatial(SL)
  }
  return(SL)
}
