#' Mean Cosine and Sine
#'
#' @param x angles in degrees
#' @param w weightings
#' @param na.rm logical
#'
#' @return named two element vector
#'
#' @examples
#' \dontrun{
#' x <- rvm(100, 0, 5)
#' mean_SC(x)
#' }
mean_SC <- function(x, w = NULL, na.rm = TRUE) {
  stopifnot(any(is.numeric(x)), is.logical(na.rm))

  w <- if (is.null(w)) {
    rep(1, times = length(x))
  } else {
    w
  }

  data <- cbind(x = x, w = w)
  if (na.rm) {
    data <- data[stats::complete.cases(data), ] # remove NA values
  }

  x <- deg2rad(data[, "x"])
  w <- data[, "w"]

  Z <- sum(w)

  sinx <- w * sin(x)
  cosx <- w * cos(x)
  # sumsin <- sum(sinx)
  # sumcos <- sum(cosx)
  # meansin <- sumsin / Z
  # meancos<- sumcos / Z
  # cbind(C = meancos, S = meansin)
  #
  # sums <- colSums(cbind(cosx, sinx))
  sums <- c(sum(cosx), sum(sinx))
  setNames(sums / Z, nm = c("C", "S"))
}


#' Mean Resultant Length
#'
#' Measure of spread around the circle. It should be noted that:
#' If R=0, then the data is completely spread around the circle.
#' If R=1, the data is completely concentrated on one point.
#'
#' @param x numeric vector. Values in degrees, for which the
#' mean, median or standard deviation are required.
#' @param w (optional) Weights. A vector of positive numbers, of the same length as
#' \code{x}.
#' @param na.rm logical value indicating whether \code{NA} values in \code{x}
#' should be stripped before the computation proceeds.
#'
#' @returns numeric.
#'
#' @export
#'
#' @references Mardia, K.V. (1972). Statistics of Directional Data: Probability
#' and Mathematical Statistics. London: Academic Press.
#'
#' @examples
#' # Example data from Davis (1986), pp. 316
#' finland_stria <- c(
#'   23, 27, 53, 58, 64, 83, 85, 88, 93, 99, 100, 105, 113,
#'   113, 114, 117, 121, 123, 125, 126, 126, 126, 127, 127, 128, 128, 129, 132,
#'   132, 132, 134, 135, 137, 144, 145, 145, 146, 153, 155, 155, 155, 157, 163,
#'   165, 171, 172, 179, 181, 186, 190, 212
#' )
#' mean_resultant_length(finland_stria, w = NULL, na.rm = FALSE) # 0.800
mean_resultant_length <- function(x, w = NULL, na.rm = TRUE) {
  m <- mean_SC(x, w, na.rm)
  R <- sqrt(m["C"]^2 + m["S"]^2)
  abs(unname(R))
}

#' @title Summary Statistics of Circular Data
#'
#' @description Calculate the (weighted median) and standard deviation
#' of orientation data.
#'
#' @param x numeric vector. Values in degrees.
#' @param w (optional) Weights. A vector of positive numbers and of the same
#' length as \code{x}.
#' @param na.rm logical value indicating whether \code{NA} values in \code{x}
#' should be stripped before the computation proceeds.
#' @param axial logical. Whether the data are axial, i.e. pi-periodical
#' (`TRUE`, the default) or directional, i.e. \eqn{2 \pi}-periodical (`FALSE`).
#'
#' @importFrom stats complete.cases runif setNames
#'
#' @returns numeric vector
#'
#' @note Weighting may be the reciprocal of the data uncertainties.
#'
#' Weightings have no effect on quasi-median and quasi-quantiles if
#' `length(x) %% 2 != 1` and `length(x) %% 4 == 0`, respectively.
#'
#' @references
#' Mardia, K.V. (1972). Statistics of Directional Data: Probability and
#' Mathematical Statistics. London: Academic Press.
#'
#' Mardia, K.V., and Jupp, P.E (1999). Directional Statistics,
#' Wiley Series in Probability and Statistics. John Wiley & Sons, Inc.,
#' Hoboken, NJ, USA. \doi{10.1002/9780470316979}
#'
#' Ziegler, M. O.; Heidbach O. (2019). Manual of the Matlab Script
#' Stress2Grid v1.1. *WSM Technical Report* 19-02,
#' GFZ German Research Centre for Geosciences. \doi{10.2312/wsm.2019.002}
#'
#' Heidbach, O., Tingay, M., Barth, A., Reinecker, J., Kurfess, D., & Mueller,
#' B. (2010). Global crustal stress pattern based on the World Stress Map
#' database release 2008. *Tectonophysics* **482**, 3<U+2013>15,
#' \doi{10.1016/j.tecto.2009.07.023}
#'
#' @examples
#' x <- rvm(10, 0, 100) %% 180
#' unc <- stats::runif(100, 0, 10)
#' circular_mean(x, 1 / unc)
#' circular_var(x, 1 / unc)
#' sample_circular_dispersion(x, 1 / unc)
#' circular_sd(x, 1 / unc)
#' circular_median(x, 1 / unc)
#' circular_quantiles(x, 1 / unc)
#' circular_IQR(x, 1 / unc)
#'
#' data("san_andreas")
#' circular_mean(san_andreas$azi)
#' circular_mean(san_andreas$azi, 1 / san_andreas$unc)
#' circular_median(san_andreas$azi)
#' circular_median(san_andreas$azi, 1 / san_andreas$unc)
#' circular_quantiles(san_andreas$azi)
#' circular_quantiles(san_andreas$azi, 1 / san_andreas$unc)
#' circular_var(san_andreas$azi)
#' circular_var(san_andreas$azi, 1 / san_andreas$unc)
#' sample_circular_dispersion(san_andreas$azi, 1 / san_andreas$unc)
#'
#' data("nuvel1")
#' PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
#' sa.por <- PoR_shmax(san_andreas, PoR, "right")
#' circular_mean(sa.por$azi.PoR, 1 / san_andreas$unc)
#' circular_median(sa.por$azi.PoR, 1 / san_andreas$unc)
#' circular_var(sa.por$azi.PoR, 1 / san_andreas$unc)
#' sample_circular_dispersion(sa.por$azi.PoR, 1 / san_andreas$unc)
#' circular_quantiles(sa.por$azi.PoR, 1 / san_andreas$unc)
#' @name circle_stats
NULL

#' @rdname circle_stats
#' @export
circular_mean <- function(x, w = NULL, axial = TRUE, na.rm = TRUE) {
  if (axial) {
    f <- 2
  } else {
    f <- 1
  }

  x <- x * f
  m <- mean_SC(x, w, na.rm)
  meanx_rad <- atan2(m["S"], m["C"]) / f
  meanx_deg <- rad2deg(meanx_rad) %% (360 / f)
  unname(meanx_deg)
}
#' @rdname circle_stats
#' @export
circular_var <- function(x, w = NULL, axial = TRUE, na.rm = TRUE) {
  if (axial) x <- ax2dir(x)

  R <- mean_resultant_length(x = x, w = w, na.rm = na.rm)
  1 - R
}


var_to_sd <- function(v) {
  s <- sqrt(-2 * log(1 - v))
  rad2deg(s)
}

sd_to_var <- function(s) {
  s_rad <- deg2rad(s)
  1 - exp(-s_rad^2 / 2)
}



#' @rdname circle_stats
#' @export
circular_sd <- function(x, w = NULL, axial = TRUE, na.rm = TRUE) {
  f <- as.numeric(axial) + 1
  # mod <- 360 / f
  x <- (x * f) %% 360

  R <- mean_resultant_length(x = x, w = w, na.rm = na.rm)
  sd <- sqrt(-2 * log(R))
  rad2deg(sd / f)
}


#' @rdname circle_stats
#' @export
circular_median <- function(x, w = NULL, axial = TRUE, na.rm = TRUE) {
  if (is.null(w)) {
    w <- rep(1, times = length(x)) |> unname()
  }

  f <- as.numeric(axial) + 1
  mod <- 360 / f
  x <- deg2rad(x * f) %% (2 * pi)
  data <- cbind(x = x, w = w)
  if (na.rm) {
    data <- data[stats::complete.cases(data), ] # remove NA values
  }

  data <- data[order(data[, "x"]), ]
  x <- data[, "x"]
  w <- data[, "w"]

  n <- length(x)

  if (n %% 2 != 0) { # if odd
    m <- (n - 1) / 2
    sumsin2 <- sin(x[m + 1])
    sumcos2 <- cos(x[m + 1])
  } else { # if even
    m <- n / 2
    Z <- (w[m] + w[m + 1])
    sumsin2 <- (w[m] * sin(x[m]) + w[m + 1] * sin(x[m + 1])) / Z
    sumcos2 <- (w[m] * cos(x[m]) + w[m + 1] * cos(x[m + 1])) / Z
  }
  (atan2d(sumsin2, sumcos2) / f) %% mod
}

#' @rdname circle_stats
#' @export
circular_quantiles <- function(x, w = NULL, axial = TRUE, na.rm = TRUE) {
  f <- as.numeric(axial) + 1
  mod <- 360 / f
  x <- deg2rad(f * x) %% (2 * pi)

  if (is.null(w)) {
    w <- rep(1, times = length(x))
  }

  data <- cbind(x = x, w = w)
  if (na.rm) {
    data <- data[stats::complete.cases(data), ] # remove NA values
  }
  data <- data[order(data[, "x"]), ]


  x <- data[, "x"]
  w <- data[, "w"]
  n <- length(x)

  if (n > 3) {
    # median:
    if (n %% 2 != 0) { # if odd
      m <- (n - 1) / 2
      sumsin2 <- sin(x[m + 1])
      sumcos2 <- cos(x[m + 1])
    } else { # if even
      m <- n / 2
      Z <- w[m] + w[m + 1]
      sumsin2 <- (w[m] * sin(x[m]) + w[m + 1] * sin(x[m + 1])) / Z
      sumcos2 <- (w[m] * cos(x[m]) + w[m + 1] * cos(x[m + 1])) / Z
    }
    med <- atan2d(sumsin2, sumcos2)

    if (n %% 4 == 0) {
      m <- n / 4
      sum.sin.lq <- sin(x[m + 1])
      sum.cos.lq <- cos(x[m + 1])

      sum.sin.uq <- sin(x[3 * m + 1])
      sum.cos.uq <- cos(x[3 * m + 1])

      Zu <- Zl <- 1
    } else if (n %% 4 == 1) {
      m <- (n - 1) / 4
      sum.sin.lq <- 3 * w[m] * sin(x[m]) + w[m + 1] * sin(x[m + 1])
      sum.cos.lq <- 3 * w[m] * cos(x[m]) + w[m + 1] * cos(x[m + 1])

      sum.sin.uq <- 3 * w[3 * m] * sin(x[3 * m]) + w[3 * m + 1] * sin(x[3 * m + 1])
      sum.cos.uq <- 3 * w[3 * m] * cos(x[3 * m]) + w[3 * m + 1] * cos(x[3 * m + 1])

      Zl <- w[m] + w[m + 1]
      Zu <- w[3 * m] + w[3 * m + 1]
    } else if (n %% 4 == 2) {
      m <- (n - 2) / 4
      sum.sin.lq <- w[m] * sin(x[m]) + w[m + 1] * sin(x[m + 1])
      sum.cos.lq <- w[m] * cos(x[m]) + w[m + 1] * cos(x[m + 1])

      sum.sin.uq <- w[3 * m] * sin(x[3 * m]) + w[3 * m + 1] * sin(x[3 * m + 1])
      sum.cos.uq <- w[3 * m] * cos(x[3 * m]) + w[3 * m + 1] * cos(x[3 * m + 1])

      Zl <- w[m] + w[m + 1]
      Zu <- w[3 * m] + w[3 * m + 1]
    } else { # if (n %% 4 == 3) {
      m <- (n - 2) / 4
      sum.sin.lq <- w[m] * sin(x[m]) + 3 * w[m + 1] * sin(x[m + 1])
      sum.cos.lq <- w[m] * cos(x[m]) + 3 * w[m + 1] * cos(x[m + 1])

      sum.sin.uq <- w[3 * m] * sin(x[3 * m]) + 3 * w[3 * m + 1] * sin(x[3 * m + 1])
      sum.cos.uq <- w[3 * m] * cos(x[3 * m]) + 3 * w[3 * m + 1] * cos(x[3 * m + 1])

      Zl <- w[m] + w[m + 1]
      Zu <- w[3 * m] + w[3 * m + 1]
    }
    mean.sin.lq <- sum.sin.lq / Zl
    mean.cos.lq <- sum.cos.lq / Zl

    mean.sin.uq <- sum.sin.uq / Zu
    mean.cos.uq <- sum.cos.uq / Zu

    lq <- atan2d(mean.sin.lq, mean.cos.lq)
    uq <- atan2d(mean.sin.uq, mean.cos.uq)

    res <- c(lq, med, uq) / f
    setNames(res %% mod, nm = c("25%", "50%", "75%"))
  } else {
    message("x needs more than 3 values")
    return(NULL)
  }
}

#' @rdname circle_stats
#' @export
circular_IQR <- function(x, w = NULL, axial = TRUE, na.rm = TRUE) {
  quantiles <- circular_quantiles(x, w, axial, na.rm)
  res <- deviation_norm(quantiles[3], quantiles[1])
  unname(res)
}

#' @rdname circle_stats
#' @export
sample_circular_dispersion <- function(x, w = NULL, axial = TRUE, na.rm = TRUE) {
  if (axial) x <- ax2dir(x)
  Rbar2 <- mean_resultant_length(2 * x, w = w)
  Rbar <- mean_resultant_length(x, w = w)
  (1 - Rbar2) / (2 * Rbar^2)
}


#' Circular Distance and Dispersion
#'
#' Circular distance between two angles and circular dispersion of angles
#' about a specified angle.
#'
#' @param x,y vectors of numeric values in degrees. `length(y)` is either
#' `1` or `length(x)`
#' @param w,w.y (optional) Weights. A vector of positive numbers and of the same
#' length as \code{x}. `w.y` is the (optional) weight of `y`.
#' @param norm logical. Whether the dispersion should be normalized by the
#' maximum possible angular difference.
#' @param axial logical. Whether the data are axial, i.e. pi-periodical
#' (`TRUE`, the default) or directional, i.e. \eqn{2 \pi}-periodical (`FALSE`).
#' @param na.rm logical. Whether \code{NA} values in \code{x}
#' should be stripped before the computation proceeds.
#'
#' @details
#' [circular_distance_alt()] and [circular_dispersion_alt()] are the alternative
#' versions in Mardia and Jupp (2000), pp. 19-20.
#' The alternative dispersion has a minimum at the sample median.
#'
#'
#' @references Mardia, K.V. (1972). Statistics of Directional Data: Probability
#' and Mathematical Statistics. London: Academic Press.
#'
#' Mardia, K.V., and Jupp, P.E (1999). Directional Statistics,
#' Wiley Series in Probability and Statistics. John Wiley & Sons, Inc.,
#' Hoboken, NJ, USA. \doi{10.1002/9780470316979}
#'
#' @importFrom stats complete.cases
#'
#' @returns `circular_distance`returns a numeric vector of positive numbers,
#' `circular_dispersion`returns a positive number.
#'
#' @note
#' If `y` is `NULL`, than the circular variance is returned.
#'
#' @seealso [circular_mean()], [circular_var()].
#'
#' @name dispersion
#'
#' @examples
#' a <- c(0, 2, 359, 6, 354)
#' circular_distance(a, 10) # distance to single value
#'
#' b <- a + 90
#' circular_distance(a, b) # distance to multiple values
#'
#' data("nuvel1")
#' PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
#' sa.por <- PoR_shmax(san_andreas, PoR, "right")
#' circular_dispersion(sa.por$azi.PoR, y = 135)
#' circular_dispersion(sa.por$azi.PoR, y = 135, w = 1 / san_andreas$unc)
NULL

#' @rdname dispersion
#' @export
circular_distance <- function(x, y, axial = TRUE, na.rm = TRUE) {
  f <- as.numeric(axial) + 1

  stopifnot(length(y) == 1 | length(y) == length(x))
  if (length(y) == 1) {
    y <- rep(y, length(x))
  }
  if (length(x) == 1) {
    x <- rep(x, length(y))
  }

  if (length(x) > 1) {
    data <- cbind(x = x, y = y)
    if (na.rm) {
      data <- data[stats::complete.cases(data), ] # remove NA values
    }

    x <- data[, "x"]
    y <- data[, "y"]
  }

  diff <- x - y
  (1 - cosd(f * diff)) / f
}

#' @rdname dispersion
#' @export
circular_dispersion <- function(x, y = NULL, w = NULL, w.y = NULL, norm = FALSE, axial = TRUE, na.rm = TRUE) {
  if (is.null(y)) {
    circular_var(x, w, axial, na.rm)
  } else {
    stopifnot(length(y) == 1 | length(y) == length(x))

    if (is.null(w)) {
      w <- rep(1, times = length(x))
    }
    if (is.null(w.y)) {
      w.y <- rep(1, times = length(x))
    }
    if (length(y) == 1) {
      y <- rep(y, times = length(x))
    }

    data <- cbind(x = x, w = w, y = y, w.y = w.y)
    if (na.rm) {
      data <- data[stats::complete.cases(data), ] # remove NA values
    }

    x <- data[, "x"]
    w.x <- data[, "w"]
    y <- data[, "y"]
    w.y <- data[, "w.y"]

    w <- w.x * w.y

    Z <- sum(w)

    md <- ifelse(norm, 2, 1)

    cdists <- circular_distance(x, y, axial, na.rm = FALSE)
    sum(w * cdists) / (Z * md)
  }
}

#' @rdname dispersion
#' @export
circular_distance_alt <- function(x, y, axial = TRUE, na.rm = TRUE) {
  f <- as.numeric(axial) + 1

  stopifnot(length(y) == 1 | length(y) == length(x))
  if (length(y) == 1) {
    y <- rep(y, length(x))
  }

  if (length(x) > 1) {
    data <- cbind(x = x, y = y)
    if (na.rm) {
      data <- data[stats::complete.cases(data), ] # remove NA values
    }

    x <- data[, "x"]
    y <- data[, "y"]
  }

  diff <- x - y
  # min(c(diff, 360 - (diff)))
  (180 - abs(180 - abs(diff))) / f
}

#' @rdname dispersion
#' @export
circular_dispersion_alt <- function(x, y = NULL, w = NULL, w.y = NULL, norm = FALSE, axial = TRUE, na.rm = TRUE) {
  if (is.null(y)) {
    circular_var(x, w, axial, na.rm)
  } else {
    stopifnot(length(y) == 1 | length(y) == length(x))

    if (is.null(w)) {
      w <- rep(1, times = length(x))
    }
    if (is.null(w.y)) {
      w.y <- rep(1, times = length(x))
    }
    if (length(y) == 1) {
      y <- rep(y, times = length(x))
    }

    data <- cbind(x = x, w = w, y = y, w.y = w.y)
    if (na.rm) {
      data <- data[stats::complete.cases(data), ] # remove NA values
    }

    x <- data[, "x"]
    w.x <- data[, "w"]
    y <- data[, "y"]
    w.y <- data[, "w.y"]

    w <- w.x * w.y

    Z <- sum(w)

    md <- ifelse(norm, 2, 1)

    cdists <- circular_distance_alt(x, y, axial, na.rm = FALSE)
    sum(w * cdists) / (Z * md)
  }
}



#' Circular Mean Difference
#'
#' @param x numeric vector. Values in degrees.
#' @param w (optional) Weights. A vector of positive numbers and of the same
#' length as \code{x}.
#' @param na.rm logical value indicating whether \code{NA} values in \code{x}
#' should be stripped before the computation proceeds.
#' @param axial logical. Whether the data are axial, i.e. pi-periodical
#' (`TRUE`, the default) or directional, i.e. \eqn{2 \pi}-periodical (`FALSE`).
#'
#' @references
#' Mardia, K.V., and Jupp, P.E (1999). Directional Statistics,
#' Wiley Series in Probability and Statistics. John Wiley & Sons, Inc.,
#' Hoboken, NJ, USA. \doi{10.1002/9780470316979}
#'
#' @return numeric
#'
#' @examples
#' data("san_andreas")
#' circular_mean_difference(san_andreas$azi)
#' circular_mean_difference(san_andreas$azi, 1 / san_andreas$unc)
#'
#' circular_mean_difference_alt(san_andreas$azi)
#' circular_mean_difference_alt(san_andreas$azi, 1 / san_andreas$unc)
#' @name circle_mean_diff
NULL

#' @rdname circle_mean_diff
#' @export
circular_mean_difference <- function(x, w = NULL, axial = TRUE, na.rm = TRUE) {
  if (axial) x <- ax2dir(x)
  Rbar2 <- mean_resultant_length(2 * x, w = w)
  1 - Rbar2
}

#' @rdname circle_mean_diff
#' @export
circular_mean_difference_alt <- function(x, w = NULL, axial = TRUE, na.rm = TRUE) {
  f <- 1
  if (axial) {
    x <- ax2dir(x)
    f <- 2
  }
  w <- if (is.null(w)) {
    rep(1, times = length(x))
  } else {
    unname(w)
  }

  data <- cbind(x = x, w = w)
  if (na.rm) {
    data <- data[stats::complete.cases(data), ] # remove NA values
  }

  x <- data[, "x"]
  w <- data[, "w"]

  Z <- sum(w)

  n <- length(x)


  d <- matrix(nrow = n, ncol = n)
  for (j in seq_along(x)) {
    for (i in seq_along(x)) {
      diff <- x[i] - x[j]
      cdists <- (180 - abs(180 - abs(diff))) / f
      d[i, j] <- (w[i] * w[j]) * sum(cdists)
    }
  }
  (sum(d) / Z^2)
}


#' Circular Range
#'
#' Length of the smallest arc which contains all the observations.
#'
#' @param x numeric vector. Values in degrees.
#' @param na.rm logical value indicating whether \code{NA} values in \code{x}
#' should be stripped before the computation proceeds.
#' @param axial logical. Whether the data are axial, i.e. pi-periodical
#' (`TRUE`, the default) or directional, i.e. \eqn{2 \pi}-periodical (`FALSE`).
#'
#' @return numeric. angle in degrees
#' @export
#'
#' @references
#' Mardia, K.V., and Jupp, P.E (1999). Directional Statistics,
#' Wiley Series in Probability and Statistics. John Wiley & Sons, Inc.,
#' Hoboken, NJ, USA. \doi{10.1002/9780470316979}
#'
#' @examples
#' roulette <- c(43, 45, 52, 61, 75, 88, 88, 279, 357)
#' circular_range(roulette, axial = FALSE)
#'
#' data("san_andreas")
#' circular_range(san_andreas$azi)
circular_range <- function(x, axial = TRUE, na.rm = TRUE) {
  f <- as.numeric(axial) + 1
  mod <- 360 / f

  if (na.rm) x <- na.omit(x)
  x <- (x * f) %% 360
  x <- sort(x)
  n <- length(x)

  t <- numeric(n)
  for (i in 1:(n - 1)) {
    t[i] <- x[i + 1] - x[i]
  }
  t[n] <- 360 - x[n] - x[1]

  w <- 360 - max(t)
  w / f
}



cdist2angle <- function(x, axial = TRUE) {
  f <- as.numeric(axial) + 1
  acosd(1 - f * x) / f
}




#' Error of Model's Prediction
#'
#' The maximum error in the model's predicted azimuth given the Pole of
#' rotations uncertainty and distance of the data point to the pole.
#'
#' @param dist_PoR Distance to Euler pole (great circle distance, in degree)
#' @param sigma_PoR uncertainty of the position of the Pole of rotation
#' (in degree).
#'
#' @references Ramsay, J.A. Folding and fracturing of rocks. McGraw-Hill, New York, 1967.
#'
#' @returns numeric vector. The maximum error for azimuths prediction (in degree)
#'
#' @seealso  [PoR_shmax()] and [model_shmax()] for the model's prediction, and
#' [orthodrome()] for great circle distances.
#'
#' @export
#'
#' @examples
#' prd_err(67, 1)
prd_err <- function(dist_PoR, sigma_PoR = 1) {
  x <- 2 * sind(sigma_PoR)^2
  y <- 1 + cosd(dist_PoR)
  acos_beta <- sqrt(1 - x / (sind(dist_PoR)^2) * y)
  acosd(acos_beta) / 2
}


z_score <- function(conf.level) {
  stats::qnorm(1 - (1 - conf.level) / 2)
}

#' Standard Error of Mean Direction of Circular Data
#'
#' Measure of the chance variation expected from sample to sample in estimates
#' of the mean direction.
#' It is a parametric estimate of the the circular standard error of the mean direction
#' by the particular form of the standard error for the von Mises distribution.
#' The approximated standard error of the mean direction is computed by the mean
#' resultant length and the MLE concentration parameter \eqn{\kappa}.
#'
#' @inheritParams circular_mean
#'
#' @returns numeric
#'
#' @seealso [mean_resultant_length()], [circular_mean()]
#'
#' @references
#' N.I. Fisher (1993) Statistical Analysis of Circular Data, Cambridge University Press.
#'
#' Davis (1986) Statistics and data analysis in geology. 2nd ed., John Wiley & Sons.
#'
#' @export
#'
#' @examples
#' # Example data from Davis (1986), pp. 316
#' finland_stria <- c(
#'   23, 27, 53, 58, 64, 83, 85, 88, 93, 99, 100, 105, 113,
#'   113, 114, 117, 121, 123, 125, 126, 126, 126, 127, 127, 128, 128, 129, 132,
#'   132, 132, 134, 135, 137, 144, 145, 145, 146, 153, 155, 155, 155, 157, 163,
#'   165, 171, 172, 179, 181, 186, 190, 212
#' )
#' circular_sd_error(finland_stria, axial = FALSE)
#'
#' data(san_andreas)
#' data("nuvel1")
#' PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
#' sa.por <- PoR_shmax(san_andreas, PoR, "right")
#' circular_sd_error(sa.por$azi.PoR, w = 1 / san_andreas$unc)
circular_sd_error <- function(x, w = NULL, axial = TRUE, na.rm = TRUE) {
  if (axial) {
    f <- 2
    # mod <- 90
  } else {
    f <- 1
    # mod <- 180
  }

  if (is.null(w)) {
    w <- rep(1, times = length(x))
  }
  data <- cbind(x = x, w = w)
  if (na.rm) {
    data <- data[stats::complete.cases(data), ] # remove NA values
  }
  x <- data[, "x"]
  w <- data[, "w"]

  n <- length(x)
  # n <- sum(w)

  kappa <- est.kappa(x, w = w, axial = axial, na.rm = FALSE)

  x <- (x * f) %% 360
  R <- mean_resultant_length(x, w = w, na.rm = FALSE)

  sde <- 1 / sqrt(n * R * kappa)
  return(sde)
}

#' Confidence Interval around the Mean Direction of Circular Data
#'
#' Probabilistic limit on the location of the true or population mean direction,
#' assuming that the estimation errors are normally distributed.
#'
#' @inheritParams circular_mean
#' @param conf.level Level of confidence: \eqn{(1 - \alpha \%)/100}.
#' (`0.95` by default).
#'
#' @returns Angle in degrees
#'
#' @seealso [mean_resultant_length()], [circular_sd_error()]
#'
#' @references
#' * Davis (1986) Statistics and data analysis in geology. 2nd ed., John Wiley
#' & Sons.
#' * Jammalamadaka, S. Rao and Sengupta, A. (2001). Topics in Circular
#' Statistics, Sections 3.3.3 and 3.4.1, World Scientific Press, Singapore.
#'
#' @details
#' The confidence angle gives the interval, i.e. plus and minus the confidence angle,
#' around the mean direction of a particular sample, that contains the true
#' mean direction under a given level of confidence.
#'
#' @examples
#' # Example data from Davis (1986), pp. 316
#' finland_stria <- c(
#'   23, 27, 53, 58, 64, 83, 85, 88, 93, 99, 100, 105, 113,
#'   113, 114, 117, 121, 123, 125, 126, 126, 126, 127, 127, 128, 128, 129, 132,
#'   132, 132, 134, 135, 137, 144, 145, 145, 146, 153, 155, 155, 155, 157, 163,
#'   165, 171, 172, 179, 181, 186, 190, 212
#' )
#' confidence_angle(finland_stria, axial = FALSE)
#' confidence_interval(finland_stria, axial = FALSE)
#'
#' data(san_andreas)
#' data("nuvel1")
#' PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
#' sa.por <- PoR_shmax(san_andreas, PoR, "right")
#' confidence_angle(sa.por$azi.PoR, w = 1 / san_andreas$unc)
#' confidence_interval(sa.por$azi.PoR, w = 1 / san_andreas$unc)
#' @name confidence
NULL

#' @rdname confidence
#' @export
confidence_angle <- function(x, conf.level = .95, w = NULL, axial = TRUE, na.rm = TRUE) {
  if (axial) {
    f <- 2
  } else {
    f <- 1
  }

  Z_alpha <- z_score(conf.level)
  sde <- circular_sd_error(x, w, axial, na.rm)
  asind(Z_alpha * sde) * f
}

#' @rdname confidence
#' @export
confidence_interval <- function(x, conf.level = .95, w = NULL, axial = TRUE, na.rm = TRUE) {
  conf.angle <- confidence_angle(x, conf.level, w, axial, na.rm)
  mu <- circular_mean(x, w = w, axial = axial, na.rm = na.rm)

  list(
    mu = mu,
    conf.angle = conf.angle,
    conf.interval = c(mu - conf.angle, mu + conf.angle) # %% 360
  )
}

#' Confidence Interval around the Mean Direction of Circular Data after Fisher (1993)
#'
#' For large samples (`n >=25`) i performs are parametric estimate based on
#' [sample_circular_dispersion()]. For smaller size samples, it returns a
#' bootstrap estimate.
#'
#' @inheritParams circular_mean
#' @param conf.level Level of confidence: \eqn{(1 - \alpha \%)/100}.
#' (`0.95` by default).
#' @param boot logical. Force bootstrap estimation
#' @param R integer. number of bootstrap replicates
#' @param quiet logical. Prints the used estimation (parametric or bootstrap).
#'
#' @references N.I. Fisher (1993) Statistical Analysis of Circular Data, Cambridge University Press.
#'
#' @importFrom boot boot boot.ci
#'
#' @return list
#' @export
#'
#' @examples
#' # Example data from Davis (1986), pp. 316
#' finland_stria <- c(
#'   23, 27, 53, 58, 64, 83, 85, 88, 93, 99, 100, 105, 113,
#'   113, 114, 117, 121, 123, 125, 126, 126, 126, 127, 127, 128, 128, 129, 132,
#'   132, 132, 134, 135, 137, 144, 145, 145, 146, 153, 155, 155, 155, 157, 163,
#'   165, 171, 172, 179, 181, 186, 190, 212
#' )
#' confidence_interval_fisher(finland_stria, axial = FALSE)
#' confidence_interval_fisher(finland_stria, axial = FALSE, boot = TRUE)
#'
#' data(san_andreas)
#' data("nuvel1")
#' PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
#' sa.por <- PoR_shmax(san_andreas, PoR, "right")
#' confidence_interval_fisher(sa.por$azi.PoR, w = 1 / san_andreas$unc)
#' confidence_interval_fisher(sa.por$azi.PoR, w = 1 / san_andreas$unc, boot = TRUE)
confidence_interval_fisher <- function(x, conf.level = 0.95, w = NULL, axial = TRUE, na.rm = TRUE, boot = FALSE, R = 1000L, quiet = FALSE) {
  n <- ifelse(na.rm, length(na.omit(x)), length(x))

  if (n < 25 | boot) {
    print_message <- paste("Bootstrap estimate based on", R, "replicates")
    bs_result <- boot::boot(data = x, statistic = circular_mean, R = R, w = w, axial = axial, na.rm = na.rm)
    ci <- boot::boot.ci(bs_result, conf.level, "perc")
    sde <- sd(bs_result$t)
    conf.angle <- z_score(conf.level) * sde
    conf.interval <- c(ci$percent[4], ci$percent[5])
    mu <- mean(bs_result$t)
  } else {
    print_message <- "Parametric estimate"
    disp <- sample_circular_dispersion(x = x, w = w, axial = axial, na.rm = na.rm)
    sde <- sqrt(disp / n)
    conf.angle <- asind(z_score(conf.level) * sde)
    mu <- circular_mean(x = x, w = w, axial = axial, na.rm = na.rm)
    conf.interval <- c(mu - conf.angle, mu + conf.angle)
  }
  if (!quiet) message(print_message)
  list(mu = mu, conf.angle = conf.angle, conf.interval = conf.interval)
}


circular_dispersion_i <- function(x, id, ...) {
  circular_dispersion(x$x[id], y = x$mean[id], w = x$w[id], w.y = x$w.y[id], ...)
}


#' Bootstrapped Estimates for Circular Dispersion
#'
#' Calculates bootstrapped estimates of the circular dispersion,
#' its standard error and its confidence interval.
#'
#' @param x numeric values in degrees.
#' @param y numeric. The angle(s) about which the angles `x` disperse (in degrees).
#' @param w,w.y (optional) Weights for `x` and `y`, respectively. A vector of
#' positive numbers and of the same length as \code{x}.
#' @param R The number of bootstrap replicates. positive integer
#' (1000 by default).
#' @param conf.level Level of confidence: \eqn{(1 - \alpha \%)/100}.
#' (`0.95` by default).
#' @param ... optional arguments passed to [boot::boot()]
#'
#' @importFrom boot boot boot.ci
#' @importFrom stats sd
#'
#' @returns list containing:
#' \describe{
#'  \item{`MLE`}{the maximum likelihood estimate of the circular dispersion}
#'  \item{`sde`}{standard error of MLE}
#'  \item{`CI`}{lower and upper limit of the confidence interval of MLE}
#' }
#'
#' @seealso [circular_dispersion()]
#'
#' @export
#'
#' @examples
#' data("nuvel1")
#' PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
#' sa.por <- PoR_shmax(san_andreas, PoR, "right")
#' circular_dispersion(sa.por$azi.PoR, y = 135, w = 1 / san_andreas$unc)
#' circular_dispersion_boot(sa.por$azi.PoR, y = 135, w = 1 / san_andreas$unc, R = 1000)
circular_dispersion_boot <- function(x, y = NULL, w = NULL, w.y = NULL, R = 1000, conf.level = .95, ...) {
  if (is.null(w)) {
    w <- rep(1, length(x))
  }
  if (is.null(w.y)) {
    w.y <- rep(1, length(x))
  }

  dat <- data.frame(x = x, y = y, w = w, w.y = w.y)
  cdisp <- boot::boot(dat, circular_dispersion_i, R = R, ...)
  ci <- boot::boot.ci(cdisp, conf = conf.level, type = "perc")

  return(
    list(
      "MLE" = mean(cdisp$t),
      "sde" = stats::sd(cdisp$t),
      "CI" = c(ci$percent[4], ci$percent[5])
    )
  )
}





#' Second Central Momentum
#'
#' Measures the skewness (a measure of the asymmetry of the probability
#' distribution) and the kurtosis (measure of the "tailedness" of the probability
#' distribution). Standardized versions are the skewness and kurtosis normalized
#' by the mean resultant length, Mardia 1972).
#'
#' @inheritParams circular_mean
#'
#' @return list containing
#' \describe{
#' \item{`skewness`}{second central sine momentum, i.e. the skewness}
#' \item{`std_skewness`}{standardized skewness}
#' \item{`kurtosis`}{second central cosine momentum, i.e. the kurtosis}
#' \item{`std_kurtosis`}{standardized kurtosis}
#' }
#' @export
#'
#' @details
#' Negative values of skewness indicate skewed data in counterclockwise
#' direction.
#'
#' Large kurtosis values indicate tailed, values close to `0` indicate packed
#' data.
#'
#'
#'
#' @examples
#' data("nuvel1")
#' PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
#' sa.por <- PoR_shmax(san_andreas, PoR, "right")
#' second_central_moment(sa.por$azi.PoR)
#' second_central_moment(sa.por$azi.PoR, w = 1 / san_andreas$unc)
second_central_moment <- function(x, w = NULL, axial = TRUE, na.rm = FALSE) {
  if (axial) x <- ax2dir(x)
  if (is.null(w)) {
    w <- rep(1, times = length(x))
  }

  data <- cbind(x = x, w = w)
  if (na.rm) {
    data <- data[stats::complete.cases(data), ] # remove NA values
  }

  data <- data[order(data[, "x"]), ]
  x <- data[, "x"]
  w <- data[, "w"]

  n <- length(x)
  Z <- sum(w)

  x_mean <- circular_mean(x, w, axial = FALSE, na.rm = FALSE)

  dev <- x - x_mean

  sin2_dev <- w * sind(2 * dev)
  cos_2dev <- w * cosd(2 * dev)

  b <- sum(sin2_dev) / Z
  a <- sum(cos_2dev) / Z

  R <- mean_resultant_length(x)

  s <- b / (1 - R)^(3 / 2)
  k <- (a - R^4) / (1 - R)^2

  list("skewness" = b, "std_skewness" = s, "kurtosis" = a, "std_kurtosis" = k)
}


#' Circular Sample Median and Deviation
#'
#' Sample median direction for a vector of circular data
#'
#' @param x numeric vector. Values in degrees.
#' @param na.rm logical value indicating whether \code{NA} values in \code{x}
#' should be stripped before the computation proceeds.
#' @param axial logical. Whether the data are axial, i.e. pi-periodical
#' (`TRUE`, the default) or directional, i.e. \eqn{2 \pi}-periodical (`FALSE`).
#'
#' @references
#' N.I. Fisher (1993) Statistical Analysis of Circular Data, Cambridge University Press.
#'
#' @return numeric
#' @importFrom circular circular meandeviation median.circular

#' @name sample_median
#'
#' @examples
#' x <- rvm(n = 100, mean = 0, kappa = 1)
#' circular_sample_median(x)
#' circular_sample_median_deviation(x)
#'
#' data("san_andreas")
#' circular_sample_median(san_andreas$azi)
#' circular_sample_median_deviation(san_andreas$azi)
NULL

#' @rdname sample_median
#' @export
circular_sample_median <- function(x, axial = TRUE, na.rm = TRUE) {
  if (axial) x <- ax2dir(x)
  if (na.rm) x <- na.omit(x)

  x_circular <- circular::circular(deg2rad(x))
  median <- circular::median.circular(x_circular) |>
    as.numeric() |>
    rad2deg()
  if (axial) median <- dir2ax(median)
  median
}

#' @rdname sample_median
#' @export
circular_sample_median_deviation <- function(x, axial = TRUE, na.rm = TRUE) {
  if (axial) x <- ax2dir(x)
  if (na.rm) x <- na.omit(x)

  x_circular <- circular::circular(deg2rad(x))
  md <- circular::meandeviation(x_circular) |>
    as.numeric() |>
    rad2deg()
  if (axial) md <- dir2ax(md)
  md
}

#' Circular Mode
#'
#' Angle of maximum density of a specified von Mises distribution
#'
#' @param x numeric vector. Values in degrees.
#' @param axial logical. Whether the data are axial, i.e. pi-periodical
#' (`TRUE`, the default) or directional, i.e. \eqn{2 \pi}-periodical (`FALSE`).#' @param kappa
#' @param kappa von Mises distribution concentration parameter
#' @param n the number of equally spaced points at which the density is to be estimated.
#'
#' @return numeric
#' @export
#'
#' @examples
#' x <- rvm(10, 0, 100) %% 180
#' circular_mode(x, kappa = 2)
circular_mode <- function(x, kappa, axial = TRUE, n = 512){
  density <- circular_density(x, kappa = kappa, n = n, axial = axial)

  f <- as.numeric(axial) + 1

  angles <- c(1:n)/n*360/f
  angles[which.max(density)]
}


#' Circular Summary Statistics
#'
#' Circular mean, standard deviation, variance, quasi-quantiles, mode,
#' 95% confidence angle, standardized skewness and kurtosis
#'
#' @inheritParams circular_mean
#'
#' @return named vector
#' @export
#'
#' @seealso [circular_mean()], [circular_sd()], [circular_var()],
#' [circular_quantiles()], [confidence_angle()], [second_central_moment()],
#' [circular_mode()]
#' @examples
#' data("nuvel1")
#' PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
#' sa.por <- PoR_shmax(san_andreas, PoR, "right")
#' circular_summary(sa.por$azi.PoR)
#' circular_summary(sa.por$azi.PoR, w = 1 / san_andreas$unc)
circular_summary <- function(x, w = NULL, axial = TRUE, na.rm = FALSE) {
  if (is.null(w)) {
    w <- rep(1, times = length(x))
  }

  data <- cbind(x = x, w = w)
  if (na.rm) {
    data <- data[stats::complete.cases(data), ] # remove NA values
  }

  data <- data[order(data[, "x"]), ]
  x <- data[, "x"]
  w <- data[, "w"]

  n <- length(x)

  x_mean <- circular_mean(x, w, axial, F)
  x_sd <- circular_sd(x, w, axial, F)
  x_var <- circular_var(x, w, axial, F)
  x_CI <- confidence_interval_fisher(x, conf.level = 0.95, w = w, axial = axial, na.rm = F, quiet = TRUE)$conf.angle
  x_quant <- circular_quantiles(x, w, axial, F)
  x_median <- circular_sample_median(x, axial, F)
  x_mode <- circular_mode(x, kappa = 2, axial = axial)
  x_sk <- second_central_moment(x, w, axial, F)
  x_R <- mean_resultant_length(ax2dir(x), w = w, F)

  setNames(
    c(n, x_mean, x_sd, x_var, x_quant[1], x_quant[2], x_quant[3], x_median, x_mode, x_CI, x_sk$std_skewness, x_sk$std_kurtosis, x_R),
    c("n", "mean", "sd", "var", "25%", "quasi-median", "75%", "median", 'mode', "95%CI", "skewness", "kurtosis", "R")
  )
}
