#' @title Euler pole object
#' @description Creates an object of the orientation of the Euler pole axis
#' @param x latitude or x coordinate of Euler pole axis
#' @param y longitude or y
#' @param z z coordinate
#' @param geo logical,\code{TRUE} (the default) if Euler pole axis is given in
#' geographical coordinates (latitude, longitude). \code{FALSE} if given in
#' Cartesian coordinates (x, y, z)
#' @param angle (optional) Angle of rotation in degrees (CCW rotation if angle
#' is positive)
#' @return An object of class \code{"euler.pole"} containing the Euler pole
#' axis in both geographical and Cartesian coordinates and the angle of rotation
#' in radians.
#' @export
#' @examples
#' euler_pole(90, 0, angle = 45)
#' euler_pole(0, 0, 1, geo = FALSE)
euler_pole <- function(x, y, z = NA, geo = TRUE, angle = NA) {
  stopifnot(is.logical(geo))
  if (!is.na(angle)) angle <- deg2rad(angle)

  if (geo) {
    cart <- geographical_to_cartesian(c(x, y))
    lat <- x
    lon <- y
    x <- cart[1]
    y <- cart[2]
    z <- cart[3]
  } else {
    pol <- cartesian_to_geographical(c(x, y, z))
    lat <- pol[1]
    lon <- pol[2]
    x <- x
    y <- y
    z <- z
  }

  ep <- data.frame(lat, lon, x, y, z, angle)
  class(ep) <- append(class(ep), "euler.pole")
  return(ep)
}

#' Vector cross product
#'
#' Vector or cross product
#'
#' @param x,y numeric vectors of length 3
#' @returns numeric vector of length 3
#' @export
#' @examples
#' vcross(c(1, 2, 3), c(4, 5, 6))
vcross <- function(x, y) {
  stopifnot(is.numeric(x) & is.numeric(y))
  stopifnot(length(x) == length(y) && length(x) == 3)

  c(x[2] * y[3] - x[3] * y[2], x[3] * y[1] -
    x[1] * y[3], x[1] * y[2] - x[2] * y[1])
}



#' Relative rotation between two rotations
#'
#' Calculates the relative rotation between two rotations, i.e. the
#' difference from rotation 1 to rotation 2.
#'
#' @param r1,r2 Objects of class \code{"euler.pole"}. First rotation is
#' \code{r1}, followed rotation \code{r2}.
#' @references Schaeben, H., Kroner, U. and Stephan, T. (2021). Euler Poles of
#' Tectonic Plates. In B. S. Daza Sagar, Q. Cheng, J. McKinley and F. Agterberg
#' (Eds.), *Encyclopedia of Mathematical Geosciences. Encyclopedia of Earth Sciences Series*
#' (pp. 1--7). Springer Nature Switzerland AG 2021.
#' \doi{10.1007/978-3-030-26050-7_435-1}
#' @returns \code{list}. Euler axes
#' (geographical coordinates) and Euler angles (in degrees)
#' @aliases rotation quaternion
#' @seealso [euler_pole()] for class \code{"euler.pole"}
#' @export
#' @examples
#' a <- euler_pole(90, 0, angle = 45)
#' b <- euler_pole(0, 0, 1, geo = FALSE, angle = -15)
#' relative_rotation(a, b)
#' relative_rotation(b, a)
relative_rotation <- function(r1, r2) {
  w1 <- r1$angle
  w2 <- r2$angle

  e1 <- c(r1$x, r1$y, r1$z)
  e2 <- c(r2$x, r2$y, r2$z)

  angle <- as.numeric(
    2 * acos(
      cos(w2 / 2) * cos(w1 / 2) + (sin(w2 / 2) * e2) %*% (sin(w1 / 2) * e1)
    )
  )

  a <- 1 / sin(angle / 2)
  b <- -cos(w2 / 2) * sin(w1 / 2) * e1 + cos(w1 / 2) * sin(w2 / 2) *
    e2 - vcross(sin(w2 / 2) * e2, sin(w1 / 2) * e1)

  axis <- a * b

  list(
    axis = cartesian_to_geographical(axis),
    angle = rad2deg(angle)
  )
}


#' @title Equivalent rotation
#' @description Transforms a sequence of rotations into a new reference system
#' @param x Object of class \code{"data.frame"} containing the Euler poles of
#' plate rotations:
#' \describe{
#'   \item{\code{plate.rot}}{Moving plate}
#'   \item{\code{lat}, \code{lon}}{coordinates of Euler pole}
#'   \item{\code{angle}}{Angle of rotation}
#'   \item{\code{plate.fix}}{Fixed plate}
#'   }
#' @param fixed plate that will be regarded as fixed. Has to be one out of
#' \code{x$plate.fix}
#' @return sequence of plate rotations in new reference system. Same object
#' class as \code{x}
#' @seealso [relative_rotation()]
#' @export
#' @examples
#' data(nuvel1) # load the NUVEL1 rotation parameters
#'
#' # all nuvel1 rotation equivalent to fixed Africa:
#' equivalent_rotation(nuvel1, fixed = "af")
equivalent_rotation <- function(x, fixed) {
  stopifnot(fixed %in% x$plate.rot)

  lat.eq <- c()
  lon.eq <- c()
  angle.eq <- c()

  fixed.plate <- subset(x, x$plate.rot == fixed)
  fixed.ep <- euler_pole(
    fixed.plate$lat,
    fixed.plate$lon,
    angle = fixed.plate$angle
  )

  if (exists(paste0("fixed.plate$plate.fix"))) {
    fixed0.plate <- subset(x, x$plate.rot == fixed.plate$plate.fix)
    fixed0.ep <- euler_pole(
      fixed0.plate$lat,
      fixed0.plate$lon,
      angle = fixed0.plate$angle
    )

    temp <- relative_rotation(fixed0.ep, fixed.ep)
    fixed.ep <- euler_pole(temp$axis[1], temp$axis[2], angle = temp$angle)
  }

  for (i in seq_along(x$plate.rot)) {
    if (x$plate.rot[i] == fixed) {
      # fixed plate has no rotation
      lat.eq[i] <- 90
      lon.eq[i] <- 0
      angle.eq[i] <- 0
    } else {
      xi.ep <- euler_pole(x$lat[i], x$lon[i], angle = x$angle[i])
      equivalent <- relative_rotation(fixed.ep, xi.ep)

      lat.eq[i] <- equivalent$axis[1]
      lon.eq[i] <- equivalent$axis[2]
      angle.eq[i] <- equivalent$angle
    }
  }
  data.frame(
    plate.rot = x$plate.rot,
    lat = lat.eq,
    lon = lon.eq,
    angle = angle.eq,
    plate.fix = fixed
  )
}

#' @title Absolute Plate Velocity
#'
#' @description Calculates the absolute angular velocity of plate motion
#' @param w Angular velocity or rate or angle of rotation
#' @param alpha Angular distance to Euler pole or small circle around Euler pole
#' @param r Radius. Default is WGS84 Earth's radius (6371.009 km)
#' @return Number (in km/Myr)
#' @seealso [earth_radius()]
#' @export
#' @examples
#' abs_vel(0.21, 0)
#' abs_vel(0.21, 45)
#' abs_vel(0.21, 90)
abs_vel <- function(w, alpha, r = earth_radius()) {
  stopifnot(is.numeric(r))
  deg2rad(w) * r * sind(alpha)
}
