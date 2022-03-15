#' @title Rotate Lines
#' @description Rotates a set of straight lines around an angle
#'
#' @param theta Angle of rotation (in degree)
#' @param p Coordinates of the lines end points
#' @param centre Coordinates of the center point of rotation
#' @return \code{matrix}
rotate_lines <- function(theta, p, centre) {
  new_x <-
    cosd(theta) * (p[, 1] - centre[1]) - sind(theta) *
      (p[, 2] - centre[2]) + centre[1]
  new_y <-
    sind(theta) * (p[, 1] - centre[1]) + cosd(theta) *
      (p[, 2] - centre[2]) + centre[2]
  matrix(c(new_x, new_y), ncol = 2)
}


#' @title Plate Stress Dummy Grid
#' @description Creates a dummy grid for small circles, great circles, and
#' loxodromes of an Euler pole
#' @param n Number of curves
#' @param angle Direction of loxodromes (in degree)
#' @param cw logical. Sense of loxodromes: \code{TRUE} for clockwise
#' loxodromes (right-lateral displaced plate boundaries). \code{FALSE} for
#' counterclockwise loxodromes (left-lateral displaced plate boundaries).
#' @return \code{data.frame}
#' @importFrom dplyr "%>%" filter mutate
#' @name dummy
#' @examples
#' smallcircle_dummy(10)
#' greatcircle_dummy(10)
#' loxodrome_dummy(10, 45, cw = FALSE)
NULL

#' @rdname dummy
#' @export
smallcircle_dummy <- function(n) {
  sm_range <- seq(0, 180, 180 / n)
  lons <- seq(-180, 180, 180 / n)

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
      "small_circle" = i
    )
    sm.df <- rbind(sm.df, sm.l)
  }
  return(sm.df)
}

#' @rdname dummy
#' @export
greatcircle_dummy <- function(n) {
  loxodrome_dummy(n, angle = 180, cw = FALSE)
}

#' @rdname dummy
#' @export
loxodrome_dummy <- function(n, angle, cw) {
  stopifnot(is.logical(cw))
  lon <- lat <- NULL
  if (cw) {
    s <- 1
  } else {
    s <- -1
  }
  lats <- seq(-180, 180, 1)

  line.dummy <- data.frame(
    lon = rep(0, length(lats)),
    lat = lats,
    line = 0
  )
  for (j in seq_along(line.dummy$lon)) {
    line.dummy.rot <-
      rotate_lines(
        theta = s * angle,
        p = cbind(
          line.dummy$lon[j],
          line.dummy$lat[j]
        ),
        centre = c(line.dummy$lon[j], 0)
      )
    loxodrome.dummy.j <- data.frame(
      lon = line.dummy.rot[, 1],
      lat = line.dummy.rot[, 2],
      loxodrome = line.dummy$line[j]
    )

    if (j == 1) {
      loxodrome.dummy <- loxodrome.dummy.j
    } else {
      loxodrome.dummy <- rbind(loxodrome.dummy, loxodrome.dummy.j)
    }
  }

  for (i in seq(-360, 360, 360 / n)) {
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

  loxodromes.filt <- loxodromes %>%
    unique() %>%
    filter(abs(lat) <= 90) %>%
    filter(abs(lon) <= 180)

  return(loxodromes.filt)
}


#' @title Theoretical Plate Tectonic Stress Paths
#'
#' @description Construct \eqn{\sigma_\text{Hmax}}{SHmax} lines that are
#' following small circles, great circles, or loxodromes of an Euler pole for
#' the relative plate motion.
#'
#' @author Tobias Stephan
#' @param x \code{data.frame} containing coordinates of Euler pole in lat, lon,
#' and rotation angle (optional)
#' @param n Number of equally spaced curves
#' @param angle Direction of loxodromes; default = 45
#' @param cw logical. Sense of loxodromes: \code{TRUE} for clockwise
#' loxodromes (right-lateral displaced plate boundaries). \code{FALSE} for
#' counterclockwise loxodromes (left-lateral displaced plate boundaries).
#' @param type Character string specifying the type of curves to export. Either
#' \code{"sm"} for small circles (default), \code{"gc"} for great circles, or
#' \code{"ld"} for loxodromes.
#' @param returnclass "sf" (default) for simple features or "sp" for spatial objects
#' @return \code{sf} or \code{SpatialLinesDataFrame}
#' @details Maximum horizontal stress can be aligned to three types of curves
#' related to relative plate motion:
#' \describe{
#' \item{Small circles}{Lines that have a constant distance to the Euler pole.
#' If x contains \code{angle}, output additionally gives absolute
#' velocity on small circle (degree/Myr -> km/Myr).}
#' \item{Great circles}{Paths of the the shortest distance between the Euler
#' pole and its antipodal position.}
#'  \item{Loxodromes}{Lines of constant bearing, i.e. curves cutting small
#'  circles at a constant angle.}
#'  }
#' @importFrom dplyr "%>%" mutate select
#' @importFrom sp Line Lines SpatialLines SpatialLinesDataFrame proj4string CRS
#' @importFrom sf st_crs st_as_sf st_set_crs st_transform st_wrap_dateline as_Spatial
#' @name stress_paths
#' @examples
#' data("nuvel1")
#' euler <- subset(nuvel1, nuvel1$plate.rot == "na") # North America relative to Pacific plate
#'
#' eulerpole_smallcircles(euler)
#' eulerpole_greatcircles(euler)
#' eulerpole_loxodromes(x = euler, angle = 45, n = 10, cw = FALSE)
#' eulerpole_loxodromes(x = euler, angle = 30, cw = TRUE)
NULL

#' @rdname stress_paths
#' @export
eulerpole_paths <- function(x, type = c("sc", "gc", "ld"), n = 10, angle, cw, returnclass = c("sf", "sp")) {
  type <- match.arg(type)
  if (type == "gc") {
    eulerpole_greatcircles(x, n, returnclass)
  } else if (type == "ld") {
    eulerpole_loxodromes(x, n, angle, cw, returnclass)
  } else {
    eulerpole_smallcircles(x, n, returnclass)
  }
}

#' @rdname stress_paths
#' @export
eulerpole_smallcircles <-
  function(x, n = 10, returnclass = c("sf", "sp")) {
    returnclass <- match.arg(returnclass)
    small_circle <- NULL
    sm.df <- smallcircle_dummy(n)
    sm_range <- unique(sm.df$small_circle)

    if (is.null(x$angle)) {
      sm_range.df <- sm_range
    } else {
      velocity <- sm.df %>%
        dplyr::mutate(abs_vel = abs_vel(w = x$angle, alpha = sm.df$small_circle)) %>%
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

    SL.wgs84 <- sp::SpatialLines(SL.list)
    SL.wgs84.df <- sp::SpatialLinesDataFrame(
      SL.wgs84,
      data.frame("n" = sm_range.df, row.names = sm_range)
    )

    SL <- PoR_to_geographical(x = sf::st_as_sf(SL.wgs84.df), ep = x)

    if (returnclass == "sp") {
      SL <- sf::as_Spatial(SL)
      suppressMessages(suppressWarnings(
        sp::proj4string(SL) <-
          sp::CRS(
            "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
          )
      ))
    }
    return(SL)
  }


#' @rdname stress_paths
#' @export
eulerpole_greatcircles <- function(x, n = 10, returnclass = c("sf", "sp")) {
  returnclass <- match.arg(returnclass)

  SL <-
    eulerpole_loxodromes(
      x,
      angle = 0,
      n = n,
      cw = TRUE,
      returnclass = returnclass
    )
  return(SL)
}

#' @rdname stress_paths
#' @export
eulerpole_loxodromes <- function(x, n = 10, angle = 45, cw, returnclass = c("sf", "sp")) {
  stopifnot(is.logical(cw))
  returnclass <- match.arg(returnclass)
  loxodrome <- NULL

  ld.df <-
    loxodrome_dummy(
      angle = abs(angle),
      n = n,
      cw = cw
    )

  ld_range <- unique(ld.df$loxodrome)
  SL.list <- list()

  for (l in unique(ld.df$loxodrome)) {
    # loop through all circles
    ld.subset <- subset(ld.df, ld.df$loxodrome == l)

    l.i <- suppressWarnings(sp::Lines(slinelist = sp::Line(
      cbind(
        "lon" = ld.subset$lon,
        "lat" = ld.subset$lat
      )
    ), ID = as.character(l)))

    suppressWarnings(SL.list[as.character(l)] <- l.i)
  }

  SL.wgs84 <- sp::SpatialLines(SL.list)
  SL.wgs84.df <- sp::SpatialLinesDataFrame(
    SL.wgs84,
    data.frame("n" = as.character(ld_range), row.names = ld_range)
  )

  SL <- PoR_to_geographical(x = sf::st_as_sf(SL.wgs84.df), ep = x)

  if (returnclass == "sp") {
    SL <- sf::as_Spatial(SL)
    suppressMessages(suppressWarnings(
      sp::proj4string(SL) <-
        sp::CRS(
          "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
        )
    ))
  }
  return(SL)
}
