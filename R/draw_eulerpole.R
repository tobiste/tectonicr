#' @title Rotate Lines
#'
#' @description Rotates a set of straight lines around an angle
#'
#' @param theta Angle of rotation (in degree)
#' @param p Coordinates of the lines end points
#' @param centre Coordinates of the center point of rotation
#'
#' @keywords internal
#'
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

#' @keywords internal
get_loxodromes <- function(lon, lat, line, theta) {
  line.dummy.rot <-
    rotate_lines(
      theta = theta,
      p = cbind(
        lon,
        lat
      ),
      centre = c(lon, 0)
    )

  data.frame(
    lon = line.dummy.rot[, 1],
    lat = line.dummy.rot[, 2],
    loxodrome = line
  )
}


#' @title Plate Stress Dummy Grid
#'
#' @description Helper functions to create a dummy grid for small circles,
#' great circles, and loxodromes of an Euler pole
#'
#' @param n Number of curves
#' @param angle Direction of loxodromes (in degree)
#' @param cw logical. Sense of loxodromes: \code{TRUE} for clockwise
#' loxodromes (right-lateral displaced plate boundaries). \code{FALSE} for
#' counterclockwise loxodromes (left-lateral displaced plate boundaries).
#'
#' @returns \code{data.frame}
#'
#' @keywords internal
#'
#' @importFrom dplyr filter mutate
#'
#' @name dummy
NULL

#' @rdname dummy
smallcircle_dummy <- function(n) {
  sm_range <- seq(0, 180, 180 / n)
  lons <- seq(-180, 180, 180 / n)

  sm.df <- cbind(
    "lon" = as.numeric(),
    "lat" = as.numeric(),
    "small_circle" = as.numeric()
  )

  lapply(sm_range, function(x) {
    lat <- x - 90
    sm.l <- data.frame(
      "lat" = rep(lat, length(lons)),
      "lon" = lons,
      "small_circle" = x
    )
  }) |>
    dplyr::bind_rows()
}

#' @rdname dummy
greatcircle_dummy <- function(n) {
  loxodrome_dummy(n, angle = 180, cw = FALSE)
}

#' @rdname dummy
loxodrome_dummy <- function(n, angle, cw) {
  stopifnot(is.logical(cw))
  lon <- lat <- numeric()
  s <- ifelse(cw, -1, 1)
  lats <- seq(-180, 180, 1)

  line.dummy <- data.frame(
    lon = rep(0, length(lats)),
    lat = lats,
    line = 0
  )

  lx <- mapply(FUN = get_loxodromes, lon = line.dummy$lon, lat = line.dummy$lat, line = line.dummy$line, theta = s * angle)
  loxodrome.dummy <- data.frame(lon = as.numeric(lx[1, ]), lat = as.numeric(lx[2, ]), loxodrome = as.numeric(lx[3, ]))


  # for (i in seq(-360, 360, 360 / n)) {
  #   line.i <- loxodrome.dummy |>
  #     mutate(
  #       lon = lon - i,
  #       loxodrome = i
  #     )
  #
  #   if (i == -360) {
  #     loxodromes <- line.i
  #   } else {
  #     loxodromes <- rbind(loxodromes, line.i)
  #   }
  # }
  lapply(seq(-360, 360, 360 / n), function(x) {
    mutate(loxodrome.dummy,
      lon = lon - x,
      loxodrome = x
    )
  }) |>
    dplyr::bind_rows() |>
    unique() |>
    dplyr::filter(abs(lat) <= 90) |>
    dplyr::filter(abs(lon) <= 180)
}

#' @title Theoretical Plate Tectonic Stress Paths
#'
#' @description Construct \eqn{\sigma_{Hmax}}{SHmax} lines that are
#' following small circles, great circles, or loxodromes of an Euler pole for
#' the relative plate motion.
#'
#' @author Tobias Stephan
#'
#' @param x Either an object of class \code{"euler.pole"} or \code{"data.frame"}
#' containing coordinates of Euler pole in lat, lon, and rotation angle
#' (optional).
#' @param n Number of equally spaced curves; `n = 10` by default (angular
#' distance between curves: `180 / n`)
#' @param angle Direction of loxodromes; `angle = 45` by default.
#' @param cw logical. Sense of loxodromes: `TRUE` for clockwise
#' loxodromes (left-lateral displaced plate boundaries). `FALSE` for
#' counterclockwise loxodromes (right-lateral displaced plate boundaries).
#' @param type Character string specifying the type of curves to export. Either
#' \code{"sm"} for small circles (default), \code{"gc"} for great circles, or
#' \code{"ld"} for loxodromes.
#'
#' @returns `sf` object
#'
#' @details Maximum horizontal stress can be aligned to three types of curves
#' related to relative plate motion:
#' \describe{
#' \item{Small circles}{Lines that have a constant distance to the Euler pole.
#' If `x` contains `angle`, output additionally gives absolute
#' velocity on small circle (degree/Myr -> km/Myr).}
#' \item{Great circles}{Paths of the shortest distance between the Euler
#' pole and its antipodal position.}
#'  \item{Loxodromes}{Lines of constant bearing, i.e. curves cutting small
#'  circles at a constant angle.}
#'  }
#'
#' @importFrom dplyr mutate select summarise group_by rename
#' @importFrom sf st_crs st_as_sf st_set_crs st_transform as_Spatial st_cast
#' @importFrom smoothr densify
#'
#' @name stress_paths
#'
#' @examples
#' data("nuvel1")
#' por <- subset(nuvel1, nuvel1$plate.rot == "na") # North America relative to
#' # Pacific plate
#'
#' eulerpole_smallcircles(por)
#' eulerpole_greatcircles(por)
#' eulerpole_loxodromes(x = por, angle = 45, n = 10, cw = FALSE)
#' eulerpole_loxodromes(x = por, angle = 30, cw = TRUE)
#' eulerpole_smallcircles(data.frame(lat = 30, lon = 10))
NULL

#' @rdname stress_paths
#' @export
eulerpole_paths <- function(x, type = c("sc", "gc", "ld"), n = 10, angle = 45, cw) {
  stopifnot(is.data.frame(x), dim(x)[1] > 0)
  type <- match.arg(type)
  # if (type == "gc") {
  #   eulerpole_greatcircles(x, n)
  # } else if (type == "ld") {
  #   eulerpole_loxodromes(x, n, angle, cw)
  # } else {
  #   eulerpole_smallcircles(x, n)
  # }
  switch(type,
    "gc" = eulerpole_greatcircles(x, n),
    "ld" = eulerpole_loxodromes(x, n, angle, cw),
    "sc" = eulerpole_smallcircles(x, n)
  )
}

#' @rdname stress_paths
#' @export
eulerpole_smallcircles <-
  function(x, n = 10) {
    stopifnot(is.data.frame(x), dim(x)[1] > 0)
    small_circle <- numeric()
    # d <- NULL
    sm.df <- smallcircle_dummy(n)

    sm.sf <- sm.df |>
      sf::st_as_sf(coords = c("lon", "lat")) |>
      dplyr::group_by(small_circle) |>
      dplyr::summarise(do_union = FALSE) |>
      sf::st_cast("MULTILINESTRING", warn = FALSE) |>
      smoothr::densify()

    sm.sf <- dplyr::mutate(sm.sf, d = ifelse(
      small_circle < 90, -1 * small_circle, 180 - small_circle
    ))

    if ("angle" %in% colnames(x)) {
      if (!is.na(x$angle)) {
        sm.sf <- sm.sf |>
          dplyr::mutate(abs_vel = abs_vel(w = x$angle, alpha = small_circle))
      }
    }

    sm.sf <- sm.sf |> dplyr::select(-small_circle)

    PoR_to_geographical_sf(x = sf::st_as_sf(sm.sf), PoR = x) |>
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
  stopifnot(is.data.frame(x), dim(x)[1] > 0, abs(angle) != 90, is.logical(cw))
  loxodrome <- numeric()

  ld.df <-
    loxodrome_dummy(
      angle = abs(angle),
      n = n,
      cw = cw
    )

  ld.sf <- ld.df |>
    sf::st_as_sf(coords = c("lon", "lat")) |>
    dplyr::group_by(loxodrome) |>
    dplyr::summarise(do_union = FALSE) |>
    sf::st_cast("MULTILINESTRING", warn = FALSE) |>
    smoothr::densify()

  ld.sf <- ld.sf |>
    dplyr::mutate(loxodrome = loxodrome %% 180) |>
    dplyr::rename(d = loxodrome)

  PoR_to_geographical_sf(x = sf::st_as_sf(ld.sf), PoR = x) |>
    sf::st_wrap_dateline(
      options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"),
      quiet = TRUE
    )
}
