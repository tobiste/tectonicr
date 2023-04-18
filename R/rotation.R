#' @title Euler pole object
#' @description Creates an object of the orientation of the Euler pole axis
#' @param x latitude or x coordinate of Euler pole axis
#' @param y longitude or y
#' @param z z coordinate
#' @param geo logical, `TRUE` (the default) if Euler pole axis is given in
#' geographical coordinates (latitude, longitude). `FALSE` if given in
#' Cartesian coordinates (`x`, `y`, `z`)
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

#' Check if object is euler.pole
#' @param x object of class \code{"euler.pole"}
#' @returns logical
is.euler <- function(x) {
  inherits(x, "euler.pole")
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
  stopifnot(is.numeric(x), is.numeric(y), length(x) == length(y), length(x) == 3)

  c(x[2] * y[3] - x[3] * y[2], x[3] * y[1] -
    x[1] * y[3], x[1] * y[2] - x[2] * y[1])
}



#' Relative rotation between two rotations
#'
#' Calculates the relative rotation between two rotations, i.e. the
#' difference from rotation 1 to rotation 2.
#'
#' @param r1,r2 Objects of class \code{"euler.pole"}. First rotation is
#' `r1`, followed rotation `r2`.
#' @references Schaeben, H., Kroner, U. and Stephan, T. (2021). Euler Poles of
#' Tectonic Plates. In B. S. Daza Sagar, Q. Cheng, J. McKinley and F. Agterberg
#' (Eds.), *Encyclopedia of Mathematical Geosciences. Encyclopedia of Earth Sciences Series*
#' (pp. 1--7). Springer Nature Switzerland AG 2021.
#' doi: 10.1007/978-3-030-26050-7_435-1.
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

#' Helper function to Equivalent rotation
#'
#' @param plate.rot,fixed character or numeric
#' @param lat,lon,angle numeric
#' @param fixed.ep data.frame
#'
#' @seealso [equivalent_rotation()]
get_relrot <- function(plate.rot, lat, lon, angle, fixed, fixed.ep) {
  if (plate.rot == fixed) {
    # fixed plate has no rotation
    lat.eq <- 90
    lon.eq <- 0
    angle.eq <- 0
  } else {
    xi.ep <- euler_pole(lat, lon, angle = angle)
    equivalent <- relative_rotation(fixed.ep, xi.ep)

    lat.eq <- equivalent$axis[1]
    lon.eq <- equivalent$axis[2]
    angle.eq <- equivalent$angle
  }
  return(c(lat.eq, lon.eq, angle.eq))
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
#' @param rot (optional) plate that will be regarded as rotating. Has to be one out of
#' \code{x$plate.rot}.
#' @return sequence of plate rotations in new reference system. Same object
#' class as \code{x}
#' @seealso [relative_rotation()]
#' @export
#' @examples
#' data(nuvel1) # load the NUVEL1 rotation parameters
#'
#' # all nuvel1 rotation equivalent to fixed Africa:
#' equivalent_rotation(nuvel1, fixed = "af")
#' # relative plate motion between Eurasia and India:
#' equivalent_rotation(nuvel1, "eu", "in")
equivalent_rotation <- function(x, fixed, rot) {
  stopifnot(fixed %in% x$plate.rot)
  plate.rot <- NULL

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

  eq <- mapply(get_relrot, plate.rot = x$plate.rot, lat = x$lat, lon = x$lon, angle = x$angle, MoreArgs = list(fixed = fixed, fixed.ep = fixed.ep))

  res <- data.frame(
    plate.rot = x$plate.rot,
    lat = eq[1, ],
    lon = eq[2, ],
    angle = eq[3, ],
    plate.fix = fixed
  )
  if (!missing(rot)) {
    stopifnot(rot %in% res$plate.rot)

    res <- subset(res, plate.rot == rot)
  }
  return(res)
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


#' Quaternion from Euler angle-axis representation for rotations
#'
#' @param x \code{"euler.pole"} object
#' @param normalize logical. Whether a quaternion normalization should be applied (TRUE) or not (FALSE, the default).
#' @return object of class \code{"quaternion"}
euler_to_Q4 <- function(x, normalize = FALSE) {
  axis <- c(x$x, x$y, x$z)
  Sc <- cos(x$angle / 2)
  Vec <- axis * sin(x$angle / 2)
  q <- (list(Sc = c(Sc), Vec = Vec))
  class(q) <- "quaternion"
  if (normalize) {
    q <- normalize_Q4(q)
  }
  return(q)
}

#' Euler angle/axis from quaternion
#' @param q object of class \code{"quaternion"}
#' @returns \code{"euler.pole"} object
Q4_to_euler <- function(q) {
  stopifnot(is.Q4(q))
  Vec <- q$Vec
  Sc <- q$Sc

  # angle <- 2 * acos(Sc)
  # numerically more stable:
  Vec_l <- sqrt(Vec[1]^2 + Vec[2]^2 + Vec[3]^2)
  angle <- 2 * atan2(Vec_l, Sc)

  axis <- Vec / sin(angle / 2)
  euler_pole(x = axis[1], y = axis[2], z = axis[3], angle = rad2deg(angle), geo = FALSE)
}

QScVec_to_Q4 <- function(x) {
  stopifnot(is.Q4(x), x$Sc != 0)
  x.euler <- Q4_to_euler(x)
  q <- c()
  q[1] <- x$Sc
  q[2] <- x.euler$axis[1] * sin(x.euler$angle / 2)
  q[3] <- x.euler$axis[2] * sin(x.euler$angle / 2)
  q[4] <- x.euler$axis[3] * sin(x.euler$angle / 2)
  q
}

Q4_to_QScVec <- function(x, normalize = FALSE) {
  stopifnot(x[1] != 0, is.logical(normalize))
  angle <- 2 * acos(x[1])
  q1 <- x[2] / sin(angle / 2)
  q2 <- x[3] / sin(angle / 2)
  q3 <- x[4] / sin(angle / 2)

  q <- list(Sc = x[1], Vec = c(q1, q2, q3))
  class(q) <- "quaternion"
  if (normalize) {
    q <- normalize_Q4(q)
  }
  return(q)
}


#' Quaternion normalization
#' @param q quaternion
#' @return object of class \code{"quaternion"}
normalize_Q4 <- function(q) {
  q4 <- QScVec_to_Q4(q)
  q.norm <- q4 / sqrt(q4[1]^2 + q4[2]^2 + q4[3]^2 + q4[4]^2)
  q.norm <- list(Sc = q.norm[1], Vec = c(q.norm[2], q.norm[3], q.norm[4]))
  class(q.norm) <- "quaternion"
  return(q.norm)
}

#' Product of quaternions
#'
#' Helper function for multiplication of two quaternions.
#' Concatenation of two rotations R1 followed by R2
#'
#' @note Multiplication is not commutative.
#'
#' @param q1,q2 two objects of class \code{"quaternion"}. first rotation R1 expressed by q1 followed by second rotation R2 expressed by q2
#' @param normalize logical. Whether a quaternion normalization should be applied (TRUE) or not (FALSE, the default).
#' @returns object of class \code{"quaternion"}
product_Q4 <- function(q1, q2, normalize = FALSE) {
  stopifnot(is.Q4(q1), is.Q4(q2), is.logical(normalize))
  Sc <- q2$Sc * q1$Sc - (q2$Vec %*% q1$Vec)
  Vec <- q1$Sc * q2$Vec + q2$Sc * q1$Vec + vcross(q2$Vec, q1$Vec)
  q <- list(Sc = c(Sc), Vec = Vec)
  class(q) <- "quaternion"
  if (normalize) {
    q <- normalize_Q4(q)
  }
  return(q)
}


#' Conjugation of a Quaternion
#'
#' Inverse rotation given by conjugated quaternion
#'
#' @param q object of class \code{"quaternion"}
#' @param normalize logical. Whether a quaternion normalization should be applied (TRUE) or not (FALSE, the default).
#' @return object of class \code{"quaternion"}
conjugate_Q4 <- function(q, normalize = FALSE) {
  stopifnot(is.Q4(q), is.logical(normalize))
  q <- list(Sc = q$Sc, Vec = -q$Vec)
  class(q) <- "quaternion"
  if (normalize) {
    q <- normalize_Q4(q)
  }
  return(q)
}

#' Check if object is quaternion
#' @param x object of class \code{"quaternion"}
#' @returns logical
is.Q4 <- function(x) {
  class(x) == "quaternion"
}

#' Rotation of a vector by a quaternion
#' @param q object of class \code{"quaternion"}
#' @param p three-column vector (Cartesian coordinates) of unit length
#' @returns three-column vector (Cartesian coordinates) of unit length
rotation_Q4 <- function(q, p) {
  stopifnot(is.Q4(q))

  q.euler <- Q4_to_euler(q)
  # Rodrigues Formula
  # p * cos(q.euler$angle) + vcross(q.euler$axis, p) * sin(q.euler$angle) + q.euler$axis * as.numeric(q.euler$axis %*% p) * (1 - cos(q.euler$axis))

  # faster:
  axis <- c(q.euler$x, q.euler$y, q.euler$z)
  angle <- q.euler$angle
  vap <- vcross(axis, p)
  p + vap * sin(angle) + vcross(axis, vap) * 2 * (sin(angle / 2))^2
}
