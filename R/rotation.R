#' @title Rotation Matrix
#'
#' @description Calculates the rotation matrix using the rotation axis and the
#' angle of rotation
#' @param n Rotation axis (in Cartesian coordinates). Can be a vector of three
#' numbers or a matrix of 3 columns (x, y, z)
#' @param alpha Rotation angle in degrees
#' @return \code{matrix}
#' @seealso [rotation_axis()] and [rotation_angle()] to extract the axis and the
#' angle from the rotation matrix.
#' @export
#' @examples
#' w <- c(0, 1, 0)
#' rotation_matrix(w, 90)
rotation_matrix <- function(n, alpha) {
  stopifnot(is.numeric(n) & is.numeric(alpha))

  n <- n / sqrt(sum(n^2)) # unit vector
  R <- matrix(nrow = 3, ncol = 3)
  R[1, 1] <- n[1]^2 * (1 - cosd(alpha)) + cosd(alpha)
  R[1, 2] <- n[1] * n[2] * (1 - cosd(alpha)) - n[3] * sind(alpha)
  R[1, 3] <- n[1] * n[3] * (1 - cosd(alpha)) + n[2] * sind(alpha)
  R[2, 1] <- n[2] * n[1] * (1 - cosd(alpha)) + n[3] * sind(alpha)
  R[2, 2] <- n[2]^2 * (1 - cosd(alpha)) + cosd(alpha)
  R[2, 3] <- n[2] * n[3] * (1 - cosd(alpha)) - n[1] * sind(alpha)
  R[3, 1] <- n[3] * n[1] * (1 - cosd(alpha)) - n[2] * sind(alpha)
  R[3, 2] <- n[3] * n[2] * (1 - cosd(alpha)) + n[1] * sind(alpha)
  R[3, 3] <- n[3]^2 * (1 - cosd(alpha)) + cosd(alpha)
  R
}


#' @title Rotation angle from rotation matrix
#' @description Extracts the rotation angle from rotation matrix
#' @param A 3x3 matrix
#' @note Infinitesimal small rotation (i.e. small angles) will cause an Error
#' due to round-off errors. In order to avoid the error, these rotations will
#' be treated as equal rotations. The function will print a warning message when
#' this is the case.
#' @return numeric angle in degree
#' @seealso [rotation_axis()] to extract the axis of the rotation matrix.
#' @export
#' @examples
#' w <- c(0, 1, 0)
#' rot <- rotation_matrix(w, 90)
#'
#' rotation_angle(rot)
rotation_angle <- function(A) {
  stopifnot(is.matrix(A))

  a <- (sum(diag(A)) - 1) / 2
  if (a >= 1) {
    warning("introduces round-off")
    a <- 1
  }
  psi <- acos(a)
  rad2deg(psi)
}

#' @title Rotation axis from rotation matrix
#' @description Extracts the rotation axis from rotation matrix
#' @param A 3x3 matrix
#' @return vector
#' @seealso [rotation_angle()] to extract the angle of the rotation matrix.
#' @export
#' @examples
#' w <- c(0, 1, 0)
#' rot <- rotation_matrix(w, 90)
#'
#' rotation_axis(rot)
rotation_axis <- function(A) {
  stopifnot(is.matrix(A))

  psi <- rotation_angle(A)
  e1 <- (A[3, 2] - A[2, 3]) / 2 * sind(psi)
  e2 <- (A[1, 3] - A[3, 1]) / 2 * sind(psi)
  e3 <- (A[2, 1] - A[1, 2]) / 2 * sind(psi)
  c(e1, e2, e3)
}

#' @title Euler pole object
#' @description Creates an object of the orientation of the Euler pole axis
#' @param x latitude or x coordinate of Euler pole axis
#' @param y longitude or y
#' @param z z coordinate
#' @param geo logical,\code{TRUE} (the default) if Euler pole axis is given in
#' geographical coordinates (latitude, longitude). \code{FALSE} if given in
#' cartesian coordinates (x, y, z)
#' @return An object of class \code{"euler.pole"}
#' containing the Euler pole axis in both geographical and cartesian coordinates.
#' @export
#' @examples
#' euler_pole(90, 0)
#' euler_pole(0, 0, 1, geo = FALSE)
euler_pole <- function(x, y, z = NA, geo = TRUE) {
  stopifnot(is.logical(geo))

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

  ep <- data.frame(lat, lon, x, y, z)
  class(ep) <- append(class(ep), "euler.pole")
  return(ep)
}

#' @title Euler rotation matrix
#' @description Creates a matrix from the given set of values.
#' @param ep An object of class \code{"euler.pole"}, or a data.frame
#' (containing lon and lat)
#' @param psi Angle of rotation in degree
#' @return \code{matrix}
#' @references Greiner, B. (1999). Euler rotations in plate-tectonic
#' reconstructions. *Computers and Geosciences*, 25(3), 209--216.
#' \doi{10.1016/S0098-3004(98)00160-5}
#' @export
#' @examples
#' ep <- euler_pole(90, 0)
#' euler_rot(ep, psi = 45)
euler_rot <- function(ep, psi) {
  stopifnot(is.numeric(psi))

  mat <- matrix(nrow = 3, ncol = 3)
  mat[1, 1] <- sind(ep$lat) * cosd(ep$lon)
  mat[1, 2] <- sind(ep$lat) * sind(ep$lon)
  mat[1, 3] <- -cosd(ep$lat)
  mat[2, 1] <- -sind(ep$lon)
  mat[2, 2] <- cosd(ep$lon)
  mat[2, 3] <- 0
  mat[3, 1] <- cosd(ep$lat) * cosd(ep$lon)
  mat[3, 2] <- cosd(ep$lat) * sind(ep$lon)
  mat[3, 3] <- sind(ep$lat)

  R <- matrix(nrow = 3, ncol = 3)
  R[1, 1] <- cosd(psi)
  R[1, 2] <- -sind(psi)
  R[1, 3] <- 0
  R[2, 1] <- sind(psi)
  R[2, 2] <- cosd(psi)
  R[2, 3] <- 0
  R[3, 1] <- 0
  R[3, 2] <- 0
  R[3, 3] <- 1

  solve(mat) %*% R %*% mat
}


#' @title Euler axis and angle from Euler matrix
#' @description Extracts the coordinates of an Euler pole and the angle from a
#' Euler matrix
#' @param A 3x3 matrix
#' @details If there is no rotation (i.e. angle i zero), the coordinates of the
#' axis are equal to Earth's spin axis (according to the GPLATES convention).
#' @return \code{list} with the following objects
#' \describe{
#' \item{pole}{object of class \code{"euler.pole"}}
#' \item{psi}{numeric Euler angle in degree}
#' }
#' @export
#' @examples
#' ep <- euler_pole(90, 0)
#' er <- euler_rot(ep, psi = 45)
#'
#' euler_from_rot(er)
euler_from_rot <- function(A) {
  stopifnot(is.matrix(A))

  psi <- rotation_angle(A)
  if (psi != 0) {
    ra <- rotation_axis(A)
  } else {
    ra <- c()
    ra[1] <- 0
    ra[2] <- 0
    ra[3] <- 1
  }
  ep <- euler_pole(ra[1], ra[2], ra[3], geo = FALSE)
  list(pole = ep, psi = psi)
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
  fixed.ep <- euler_pole(fixed.plate$lat, fixed.plate$lon)

  if (exists(paste0("fixed.plate$plate.fix"))) {
    fixed0.plate <- subset(x, x$plate.rot == fixed.plate$plate.fix)
    fixed0.ep <- euler_pole(fixed0.plate$lat, fixed0.plate$lon)
    fixed0.rot <- euler_rot(fixed0.ep, -fixed0.plate$angle)

    fixed.rot <-
      euler_rot(fixed.ep, fixed.plate$angle) %*% fixed0.rot
  } else {
    fixed.rot <- euler_rot(fixed.ep, -fixed.plate$angle)
  }

  for (i in seq_along(x$plate.rot)) {
    if (x$plate.rot[i] == fixed) {
      # fixed plate has no rotation
      lat.eq[i] <- 90
      lon.eq[i] <- 0
      angle.eq[i] <- 0
    } else {
      # Composition of finite rotations
      equivalent.rot <- fixed.rot %*%
        euler_rot(euler_pole(x$lat[i], x$lon[i]), x$angle[i])

      equivalent.ep <- euler_from_rot(equivalent.rot)

      lat.eq[i] <- equivalent.ep$pole$lat
      lon.eq[i] <- equivalent.ep$pole$lon
      angle.eq[i] <- equivalent.ep$psi
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
#' @param r Radius. Default is Earth's radius (6371.00887714 km)
#' @return Number (in km/Myr)
#' @export
#' @examples
#' abs_vel(0.21, 0)
#' abs_vel(0.21, 45)
#' abs_vel(0.21, 90)
abs_vel <- function(w, alpha, r = 6371.00887714) {
  stopifnot(is.numeric(w) & is.numeric(alpha))
  deg2rad(w) * r * sind(alpha)
}
