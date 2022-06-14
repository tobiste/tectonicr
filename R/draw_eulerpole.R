#' @title Rotate Lines
#' @description Rotates a set of straight lines around an angle
#'
#' @param theta Angle of rotation (in degree)
#' @param p Coordinates of the lines end points
#' @param centre Coordinates of the center point of rotation
#' @keywords internal
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
#' @description Helper functions to create a dummy grid for small circles,
#' great circles, and loxodromes of an Euler pole
#' @param n Number of curves
#' @param angle Direction of loxodromes (in degree)
#' @param cw logical. Sense of loxodromes: \code{TRUE} for clockwise
#' loxodromes (right-lateral displaced plate boundaries). \code{FALSE} for
#' counterclockwise loxodromes (left-lateral displaced plate boundaries).
#' @return \code{data.frame}
#' @keywords internal
#' @importFrom dplyr filter mutate
#' @importFrom magrittr %>%
#' @name dummy
NULL

#' @rdname dummy
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
greatcircle_dummy <- function(n) {
  loxodrome_dummy(n, angle = 180, cw = FALSE)
}

#' @rdname dummy
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

  loxodromes %>%
    unique() %>%
    filter(abs(lat) <= 90) %>%
    filter(abs(lon) <= 180)
}

#' @title Theoretical Plate Tectonic Stress Paths
#'
#' @description Construct \eqn{\sigma_\text{Hmax}}{SHmax} lines that are
#' following small circles, great circles, or loxodromes of an Euler pole for
#' the relative plate motion.
#'
#' @author Tobias Stephan
#' @param x Either an object of class \code{"euler.pole"} or \code{data.frame}
#' containing coordinates of Euler pole in lat, lon, and rotation angle
#' (optional).
#' @param n Number of equally spaced curves; n = 10 by default (angular distance
#' between curves: 180 / n)
#' @param angle Direction of loxodromes; angle = 45 by default.
#' @param cw logical. Sense of loxodromes: \code{TRUE} for clockwise
#' loxodromes (right-lateral displaced plate boundaries). \code{FALSE} for
#' counterclockwise loxodromes (left-lateral displaced plate boundaries).
#' @param type Character string specifying the type of curves to export. Either
#' \code{"sm"} for small circles (default), \code{"gc"} for great circles, or
#' \code{"ld"} for loxodromes.
#' @return \code{sf} object
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
#' @importFrom dplyr mutate select summarise group_by rename
#' @importFrom magrittr %>%
#' @importFrom sf st_crs st_as_sf st_set_crs st_transform as_Spatial st_cast
#' @note If package "smoothr" is installed, the sf objects will be "densified"
#' via [smoothr::densify()].
#' @name stress_paths
#' @examples
#' data("nuvel1")
#' euler <- subset(nuvel1, nuvel1$plate.rot == "na") # North America relative to
#' # Pacific plate
#'
#' eulerpole_smallcircles(euler)
#' eulerpole_greatcircles(euler)
#' eulerpole_loxodromes(x = euler, angle = 45, n = 10, cw = FALSE)
#' eulerpole_loxodromes(x = euler, angle = 30, cw = TRUE)
#' eulerpole_smallcircles(data.frame(lat = 30, lon = 10))
NULL

#' @rdname stress_paths
#' @export
eulerpole_paths <- function(x, type = c("sc", "gc", "ld"), n = 10, angle, cw) {
  stopifnot(is.data.frame(x))
  stopifnot(dim(x)[1] > 0)
  type <- match.arg(type)
  if (type == "gc") {
    eulerpole_greatcircles(x, n)
  } else if (type == "ld") {
    eulerpole_loxodromes(x, n, angle, cw)
  } else {
    eulerpole_smallcircles(x, n)
  }
}

#' @rdname stress_paths
#' @export
eulerpole_smallcircles <-
  function(x, n = 10) {
    stopifnot(is.data.frame(x))
    stopifnot(dim(x)[1] > 0)
    small_circle <- d <- NULL
    sm.df <- smallcircle_dummy(n)

    sm.sf <- sm.df %>%
      sf::st_as_sf(coords = c("lon", "lat")) %>%
      dplyr::group_by(small_circle) %>%
      dplyr::summarise(do_union = FALSE) %>%
      sf::st_cast("MULTILINESTRING")

    # If "smoothr" is installed, the object will be densified
    if (requireNamespace("smoothr", quietly = TRUE)) {
      sm.sf <- smoothr::densify(sm.sf)
    }
    sm.sf <- dplyr::mutate(sm.sf, d = ifelse(
      small_circle < 90, -1 * small_circle, 180 - small_circle
    ))

    if ("angle" %in% colnames(x)) {
      if (!is.na(x$angle)) {
        sm.sf <- sm.sf %>%
          dplyr::mutate(abs_vel = abs_vel(w = x$angle, alpha = small_circle))
      }
    }

    sm.sf <- sm.sf %>% dplyr::select(-small_circle)

    PoR_to_geographical(x = sf::st_as_sf(sm.sf), ep = x) %>%
      sf::st_wrap_dateline(
        options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"),
        quiet = TRUE
      )
  }


#' @rdname stress_paths
#' @export
eulerpole_greatcircles <- function(x, n = 10) {
  eulerpole_loxodromes(
    x,
    angle = 0,
    n = n,
    cw = TRUE
  )
}

#' @rdname stress_paths
#' @export
eulerpole_loxodromes <- function(x, n = 10, angle = 45, cw) {
  stopifnot(is.data.frame(x))
  stopifnot(dim(x)[1] > 0)
  stopifnot(abs(angle) != 90)

  stopifnot(is.logical(cw))
  loxodrome <- NULL

  ld.df <-
    loxodrome_dummy(
      angle = abs(angle),
      n = n,
      cw = cw
    )

  ld.sf <- ld.df %>%
    sf::st_as_sf(coords = c("lon", "lat")) %>%
    dplyr::group_by(loxodrome) %>%
    dplyr::summarise(do_union = FALSE) %>%
    sf::st_cast("MULTILINESTRING")

  # If "smoothr" is installed, the object will be densified
  if (requireNamespace("smoothr", quietly = TRUE)) {
    ld.sf <- smoothr::densify(ld.sf)
  }

  ld.sf <- ld.sf %>%
    dplyr::mutate(loxodrome = loxodrome %% 180) %>%
    dplyr::rename(d = loxodrome)

  PoR_to_geographical(x = sf::st_as_sf(ld.sf), ep = x) %>%
    sf::st_wrap_dateline(
      options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"),
      quiet = TRUE
    )
}
