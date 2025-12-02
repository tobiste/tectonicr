#' Mean Cosine and Sine
#'
#' @param x angles in degrees
#' @param w weightings
#' @param na.rm logical
#'
#' @return named two element vector
#'
#' @keywords internal
#' @examples
#' \dontrun{
#' set.seed(1)
#' x <- rvm(100, 0, 5)
#' mean_SC(x)
#' }
mean_SC <- function(x, w = NULL, na.rm = TRUE) {
  stopifnot(any(is.numeric(x)), is.logical(na.rm))

  if (is.null(w)) w <- rep(1, times = length(x))

  if (isTRUE(na.rm)) {
    keep <- !is.na(x) & !is.na(w)
    x <- x[keep]
    w <- w[keep]
  }

  x <- deg2rad(x)

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
#' @importFrom stats runif setNames
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
#' N.I. Fisher (1993) Statistical Analysis of Circular Data, Cambridge University Press.
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
#' set.seed(1)
#' x <- rvm(10, 0, 100) %% 180
#' unc <- stats::runif(100, 0, 10)
#' w <- weighting(unc)
#' circular_mean(x, w)
#' circular_var(x, w)
#' circular_sd(x, w)
#' circular_median(x, w)
#' circular_quantiles(x, w)
#' circular_IQR(x, w)
#'
#' data("san_andreas")
#' w2 <- weighting(san_andreas$unc)
#' circular_mean(san_andreas$azi)
#' circular_mean(san_andreas$azi, w2)
#' circular_median(san_andreas$azi)
#' circular_median(san_andreas$azi, w2)
#' circular_quantiles(san_andreas$azi)
#' circular_quantiles(san_andreas$azi, w2)
#' circular_var(san_andreas$azi)
#' circular_var(san_andreas$azi, w2)
#'
#' data("nuvel1")
#' PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
#' sa.por <- PoR_shmax(san_andreas, PoR, "right")
#' circular_mean(sa.por$azi.PoR, w2)
#' circular_median(sa.por$azi.PoR, w2)
#' circular_var(sa.por$azi.PoR, w2)
#' circular_quantiles(sa.por$azi.PoR, w2)
#' @name circle_stats
NULL

#' @rdname circle_stats
#' @export
circular_mean <- function(x, w = NULL, axial = TRUE, na.rm = TRUE) {
  f <- if (isTRUE(axial)) 2 else 1

  x <- x * f
  m <- mean_SC(x, w, na.rm)
  meanx_rad <- atan2(m["S"], m["C"]) / f
  meanx_deg <- rad2deg(meanx_rad) %% (360 / f)
  unname(meanx_deg)
}
#' @rdname circle_stats
#' @export
circular_var <- function(x, w = NULL, axial = TRUE, na.rm = TRUE) {
  if (isTRUE(axial)) x <- ax2dir(x)

  R <- mean_resultant_length(x = x, w = w, na.rm = na.rm)
  1 - R
}

#' @keywords internal
var_to_sd <- function(v) {
  s <- sqrt(-2 * log(1 - v))
  rad2deg(s)
}

#' @keywords internal
sd_to_var <- function(s) {
  s_rad <- deg2rad(s)
  1 - exp(-s_rad^2 / 2)
}



#' @rdname circle_stats
#' @export
circular_sd <- function(x, w = NULL, axial = TRUE, na.rm = TRUE) {
  f <- if (isTRUE(axial)) 2 else 1
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

  f <- if (isTRUE(axial)) 2 else 1
  mod <- 360 / f
  x <- deg2rad(x * f) %% (2 * pi)

  #remove NA
  if (isTRUE(na.rm)) {
    keep <- !is.na(x) & !is.na(w)
    x <- x[keep]
    w <- w[keep]
  }

  # order x
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]

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
  f <- if (isTRUE(axial)) 2 else 1
  mod <- 360 / f
  x <- deg2rad(f * x) %% (2 * pi)

  if (is.null(w)) {
    w <- rep(1, times = length(x))
  }

  # remove NA
  if (isTRUE(na.rm)) {
    keep <- !is.na(x) & !is.na(w)
    x <- x[keep]
    w <- w[keep]
  }

  # order x
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]

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
    return(rep(NA, 3))
  }
}

#' @rdname circle_stats
#' @export
circular_IQR <- function(x, w = NULL, axial = TRUE, na.rm = TRUE) {
  quantiles <- circular_quantiles(x, w, axial, na.rm)
  res <- deviation_norm(quantiles[3], quantiles[1])
  unname(res)
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
#' @param axial logical. Whether the data are axial, i.e. pi-periodical
#' (`TRUE`, the default) or directional, i.e. \eqn{2 \pi}-periodical (`FALSE`).
#' @param na.rm logical. Whether \code{NA} values in \code{x}
#' should be stripped before the computation proceeds.
#'
#' @details
#' Circular dispersion is a measure for the spread of data like the variance.
#' Dispersion measures the spread about a given angles, whereas
#' the variance measures the spread about the mean (Mardia and Jupp, 1999). When
#' `y = NULL` the dispersion is identical to the variance.
#'
#' Circular standard deviation in [circular_sd2()] is the transformed dispersion
#' instead of the variance as for [circular_sd()].
#'
#' @references Mardia, K.V. (1972). Statistics of Directional Data: Probability
#' and Mathematical Statistics. London: Academic Press.
#'
#' Mardia, K.V., and Jupp, P.E (1999). Directional Statistics,
#' Wiley Series in Probability and Statistics. John Wiley & Sons, Inc.,
#' Hoboken, NJ, USA. \doi{10.1002/9780470316979}
#'
#' @returns `circular_distance` returns a numeric vector of positive numbers,
#' `circular_dispersion` and [circular_sd2()] return a positive number.
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
#' circular_dispersion(sa.por$azi.PoR, y = 135, w = weighting(san_andreas$unc))
#' circular_sd2(sa.por$azi.PoR, y = 135, w = weighting(san_andreas$unc))
NULL

#' @rdname dispersion
#' @export
circular_distance <- function(x, y, axial = TRUE, na.rm = TRUE) {
  f <- if (isTRUE(axial)) 2 else 1
  nx <- length(x)
  ny <- length(y)

  stopifnot(ny == 1 || ny == nx)

  if (ny == 1) {
    y <- rep(y, nx)
  }
  if (nx == 1) {
    x <- rep(x, ny)
  }

  if (nx > 1) {
    if (isTRUE(na.rm)) {
      keep <- !is.na(x) & !is.na(y)
      x <- x[keep]
      y <- y[keep]
    }
  }

  diff <- x - y
  (1 - cosd(f * diff)) / 2
}

#' @rdname dispersion
#' @export
circular_dispersion <- function(x, y = NULL, w = NULL, w.y = NULL, axial = TRUE, na.rm = TRUE) {
  n <- length(x)
  if (is.null(y)) {
    circular_var(x, w, axial, na.rm)
  } else {
    stopifnot(length(y) == 1 | length(y) == n)

    if (is.null(w)) {
      w <- rep(1, times = n)
    }
    if (is.null(w.y)) {
      w.y <- rep(1, times = n)
    }
    if (length(y) == 1) {
      y <- rep(y, times = n)
    }

    # remove NA
    if (na.rm) {
      keep <- !is.na(x) & !is.na(w) & !is.na(y) & !is.na(w.y)
      x <- x[keep]
      w <- w[keep]
      y <- y[keep]
      w.y <- w.y[keep]
    }

    w <- w * w.y
    Z <- sum(w)

    cdists <- circular_distance(x, y, axial, na.rm = FALSE)

    # norm <- !axial
    # md <- ifelse(norm, 2, 1)
    # sum(w * cdists) / (Z * md)
    sum(w * cdists) / Z
  }
}

#' @rdname dispersion
#' @export
circular_sd2 <- function(x, y, w = NULL, axial = TRUE, na.rm = TRUE) {
  f <- if (isTRUE(axial)) 2 else 1

  D <- circular_dispersion(x, y, w, axial, na.rm)
  R <- 1 - D
  sd <- sqrt(-2 * log(R))
  rad2deg(sd / f)
}

#' Sample circular dispersion
#'
#' Alternative versions of variance, dispersion a distance
#' (Mardia and Jupp, 1999; pp. 19-20).
#' These alternative dispersion has a minimum at the sample median.
#'
#' @param x,y vectors of numeric values in degrees. `length(y)` is either
#' `1` or `length(x)`
#' @param w,w.y (optional) Weights. A vector of positive numbers and of the same
#' length as \code{x}. `w.y` is the (optional) weight of `y`.
#' @param axial logical. Whether the data are axial, i.e. pi-periodical
#' (`TRUE`, the default) or directional, i.e. \eqn{2 \pi}-periodical (`FALSE`).
#' @param na.rm logical. Whether \code{NA} values in \code{x}
#' should be stripped before the computation proceeds.
#'
#' @references
#' N.I. Fisher (1993) Statistical Analysis of Circular Data, Cambridge
#' University Press.
#'
#' Mardia, K.V., and Jupp, P.E (1999). Directional Statistics,
#' Wiley Series in Probability and Statistics. John Wiley & Sons, Inc.,
#' Hoboken, NJ, USA. \doi{10.1002/9780470316979}
#'
#' @name sample_dispersion
#'
#' @examples
#' a <- c(0, 2, 359, 6, 354)
#' sample_circular_distance(a, 10) # distance to single value
#'
#' b <- a + 90
#' sample_circular_distance(a, b) # distance to multiple values
#'
#' data("nuvel1")
#' PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
#' sa.por <- PoR_shmax(san_andreas, PoR, "right")
#' sample_circular_variance(sa.por$azi.PoR)
#' sample_circular_dispersion(sa.por$azi.PoR, y = 135)
#' sample_circular_dispersion(sa.por$azi.PoR, y = 135, w = weighting(san_andreas$unc))
NULL

#' @rdname sample_dispersion
#' @export
sample_circular_variance <- function(x, w = NULL, axial = TRUE) {
  # after Fisher 1996
  if (isTRUE(axial)) x <- ax2dir(x)
  Rbar2 <- mean_resultant_length(2 * x, w = w)
  Rbar <- mean_resultant_length(x, w = w)
  (1 - Rbar2) / (2 * Rbar^2)
}

#' @rdname sample_dispersion
#' @export
sample_circular_distance <- function(x, y, axial = TRUE, na.rm = TRUE) {
  f <- if (isTRUE(axial)) 2 else 1

  stopifnot(length(y) == 1 | length(y) == length(x))
  if (length(y) == 1) {
    y <- rep(y, length(x))
  }

  if (length(x) > 1) {
    # remove NA
    if (isTRUE(na.rm)) {
      keep <- !is.na(x) & !is.na(y)
      x <- x[keep]
      y <- y[keep]
    }
  }

  diff <- x - y
  # min(c(diff, 360 - (diff)))
  (180 - abs(180 - abs(diff))) / f
}

#' @rdname sample_dispersion
#' @export
sample_circular_dispersion <- function(x, y = NULL, w = NULL, w.y = NULL, axial = TRUE, na.rm = TRUE) {
  n <- length(x)

  if (is.null(y)) {
    circular_var(x, w, axial, na.rm)
  } else {
    stopifnot(length(y) == 1 | length(y) == n)


    if (is.null(w)) {
      w <- rep(1, times = n)
    }
    if (is.null(w.y)) {
      w.y <- rep(1, times = n)
    }
    if (length(y) == 1) {
      y <- rep(y, times = n)
    }

    if (isTRUE(na.rm)) {
      keep <- !is.na(x) & !is.na(w) & !is.na(y) & !is.na(w.y)
      x <- x[keep]
      w <- w[keep]
      y <- y[keep]
      w.y <- w.y[keep]
    }

    w <- w * w.y
    Z <- sum(w)

    # md <- ifelse(norm, 2, 1)

    cdists <- sample_circular_distance(x, y, axial, na.rm = FALSE)
    sum(w * cdists) / Z
  }
}

#' Circular Mean Difference
#'
#' The circular mean difference is based on the sample circular distance
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
#' @seealso [sample_circular_distance()]
#'
#' @examples
#' data("san_andreas")
#' circular_mean_difference(san_andreas$azi)
#' circular_mean_difference(san_andreas$azi, weighting(san_andreas$unc))
#'
#' circular_mean_difference_alt(san_andreas$azi)
#' circular_mean_difference_alt(san_andreas$azi, weighting(san_andreas$unc))
#' @name circle_mean_diff
NULL

#' @rdname circle_mean_diff
#' @export
circular_mean_difference <- function(x, w = NULL, axial = TRUE, na.rm = TRUE) {
  if (isTRUE(axial)) x <- ax2dir(x)
  Rbar2 <- mean_resultant_length(2 * x, w = w)
  1 - Rbar2
}

#' @rdname circle_mean_diff
#' @export
circular_mean_difference_alt <- function(x, w = NULL, axial = TRUE, na.rm = TRUE) {
  f <- 1
  if (isTRUE(axial)) {
    x <- ax2dir(x)
    f <- 2
  }
  w <- if (is.null(w)) {
    rep(1, times = length(x))
  } else {
    unname(w)
  }

  if (isTRUE(na.rm)) {
    keep <- !is.na(x) & !is.na(w)
    x <- x[keep]
    w <- w[keep]
  }

  Z <- sum(w)
  n <- length(x)

  # d <- matrix(nrow = n, ncol = n)
  # for (j in seq_along(x)) {
  #   for (i in seq_along(x)) {
  #     diff <- x[i] - x[j]
  #     cdists <- (180 - abs(180 - abs(diff))) / f
  #     d[i, j] <- (w[i] * w[j]) * sum(cdists) # why sum? cdist has only one element
  #   }
  # }

  diffmat <- outer(x, x, function(a, b) {
    a - b
  })
  cdists <- (180 - abs(180 - abs(diffmat))) / f
  wmat <- outer(w, w)
  # d <- wmat * sum(cdists)
  d <- wmat * cdists


  (sum(d) / Z^2)
}


#' Circular Range
#'
#' Length of the smallest arc which contains all the observations.
#' The circular range is based on the sample circular distance.
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
#' @seealso [sample_circular_distance()]
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
  f <- if (isTRUE(axial)) 2 else 1
  mod <- 360 / f

  if (isTRUE(na.rm)) x <- x[!is.na(x)]
  x <- (x * f) %% 360
  x <- sort(x)
  n <- length(x)

  # t <- numeric(n)
  # for (i in 1:(n - 1)) {
  #   t[i] <- x[i + 1] - x[i]
  # }
  t <- vapply(1:(n - 1), function(i) {
    x[i + 1] - x[i]
  }, numeric(1))
  t[n] <- 360 - x[n] - x[1]

  w <- 360 - max(t)
  w / f
}



cdist2angle <- function(x, axial = TRUE) {
  f <- if (isTRUE(axial)) 2 else 1
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
#'
#' # San Andreas example:
#' data("nuvel1")
#' por <- subset(nuvel1, nuvel1$plate.rot == "na") # North America relative to Pacific plate
#' data("san_andreas")
#' d <- PoR_distance(san_andreas, por)
#' prd_err(d)
prd_err <- function(dist_PoR, sigma_PoR = 1) {
  x <- 2 * sind(sigma_PoR)^2
  y <- 1 + cosd(dist_PoR)
  acos_beta <- sqrt(1 - x / (sind(dist_PoR)^2) * y)
  acosd(acos_beta) / 2
}

#' @keywords internal
#' @importFrom stats qnorm
z_score <- function(conf.level) {
  stats::qnorm(1 - (1 - conf.level) / 2)
}

#' Standard Error of Mean Direction of Circular Data
#'
#' Measure of the chance variation expected from sample to sample in estimates
#' of the mean direction (after Mardia 1972).
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
#' * Batschelet, E. (1971). Recent statistical methods for orientation data.
#' "Animal Orientation, Symposium 1970 on Wallops Island". Amer. Inst. Biol.
#' Sciences, Washington.
#' * Mardia, K.V. (1972). Statistics of Directional Data: Probability and
#' Mathematical Statistics. London: Academic Press.
#' * N.I. Fisher (1993) Statistical Analysis of Circular Data, Cambridge University Press.
#' * Davis (1986) Statistics and data analysis in geology. 2nd ed., John Wiley & Sons.
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
#' circular_sd_error(sa.por$azi.PoR, w = weighting(san_andreas$unc))
circular_sd_error <- function(x, w = NULL, axial = TRUE, na.rm = TRUE) {
  f <- if (isTRUE(axial)) 2 else 1

  if (is.null(w)) {
    w <- rep(1, times = length(x))
  }

  if (isTRUE(na.rm)) {
    keep <- !is.na(x) & !is.na(w)
    x <- x[keep]
    w <- w[keep]
  }

  n <- length(x)
  # n <- sum(w)


  x <- (x * f) %% 360
  kappa <- est.kappa(x, w = w, axial = FALSE)
  R <- mean_resultant_length(x, w = w, na.rm = FALSE)

  1 / sqrt(n * R * kappa)
}

#' Confidence Interval around the Mean Direction of Circular Data after Batschelet (1971)
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
#' * Batschelet, E. (1971). Recent statistical methods for orientation data.
#' "Animal Orientation, Symposium 1970 on Wallops Island". Amer. Inst. Biol.
#' Sciences, Washington.
#' * Mardia, K.V. (1972). Statistics of Directional Data: Probability and
#' Mathematical Statistics. London: Academic Press. (p. 146)
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
#' confidence_angle(sa.por$azi.PoR, w = weighting(san_andreas$unc))
#' confidence_interval(sa.por$azi.PoR, w = weighting(san_andreas$unc))
#' @name confidence
NULL

#' @rdname confidence
#' @export
confidence_angle <- function(x, conf.level = .95, w = NULL, axial = TRUE, na.rm = TRUE) {
  f <- if (isTRUE(axial)) 2 else 1

  Z_alpha <- z_score(conf.level)
  sde <- circular_sd_error(x, w, axial, na.rm)

  temp <- Z_alpha * sde
  if (temp > 1) temp <- 1 # I don't understand yet why sometimes sde > 1/Z_alpha (which makes asin undefined). Hence I set this term to 1 to make it work. Not ideal though...
  asind(temp) * f
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
#' confidence_interval_fisher(sa.por$azi.PoR, w = weighting(san_andreas$unc))
#' confidence_interval_fisher(sa.por$azi.PoR, w = weighting(san_andreas$unc), boot = TRUE)
confidence_interval_fisher <- function(x, conf.level = 0.95, w = NULL, axial = TRUE, na.rm = TRUE, boot = FALSE, R = 1000L, quiet = FALSE) {
  n <- ifelse(na.rm, length(x[!is.na(x)]), length(x))

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

    temp <- z_score(conf.level) * sde
    if (temp > 1) temp <- 1 # I don't understand yet why sometimes sde > 1/Z_alpha (which makes asin undefined). Hence I set this term to 1 to make it work. Not ideal though...
    conf.angle <- asind(temp)

    mu <- circular_mean(x = x, w = w, axial = axial, na.rm = na.rm)
    conf.interval <- c(mu - conf.angle, mu + conf.angle)
  }
  if (isFALSE(quiet)) message(print_message)
  list(mu = mu, conf.angle = conf.angle, conf.interval = conf.interval)
}

#' @keywords internal
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
#' circular_dispersion(sa.por$azi.PoR, y = 135, w = weighting(san_andreas$unc))
#' circular_dispersion_boot(sa.por$azi.PoR, y = 135, w = weighting(san_andreas$unc), R = 1000)
circular_dispersion_boot <- function(x, y = NULL, w = NULL, w.y = NULL, R = 1000, conf.level = .95, ...) {
  n <- length(x)
  if (is.null(w)) {
    w <- rep(1, n)
  }
  if (is.null(w.y)) {
    w.y <- rep(1, n)
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
#' by the mean resultant length (Mardia 1972).
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
#' second_central_moment(sa.por$azi.PoR, w = weighting(san_andreas$unc))
second_central_moment <- function(x, w = NULL, axial = TRUE, na.rm = FALSE) {
  if (isTRUE(axial)) x <- ax2dir(x)
  if (is.null(w)) {
    w <- rep(1, times = length(x))
  }

  #remove NA
  if (isTRUE(na.rm)) {
    keep <- !is.na(x) & !is.na(w)
    x <- x[keep]
    w <- w[keep]
  }

  # order x
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]

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


#' Sample Circular Median and Deviation
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
#' set.seed(1)
#' x <- rvm(n = 100, mean = 0, kappa = 10)
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
  if (isTRUE(axial)) x <- ax2dir(x)
  if (isTRUE(na.rm)) x <- x[!is.na(x)]

  x_circular <- circular::circular(deg2rad(x))
  median <- circular::median.circular(x_circular) |>
    as.numeric() |>
    rad2deg()
  if (isTRUE(axial)) median <- dir2ax(median)
  median
}

#' @rdname sample_median
#' @export
circular_sample_median_deviation <- function(x, axial = TRUE, na.rm = TRUE) {
  if (isTRUE(axial)) x <- ax2dir(x)
  if (isTRUE(na.rm)) x <- x[!is.na(x)]

  x_circular <- circular::circular(deg2rad(x))
  md <- circular::meandeviation(x_circular) |>
    as.numeric() |>
    rad2deg()
  if (isTRUE(axial)) md <- dir2ax(md)
  md
}

#' Circular Mode
#'
#' MLE angle (maximum density) using a von Mises distribution kernel with
#' specified concentration.
#'
#' @param x numeric vector. Values in degrees.
#' @param axial logical. Whether the data are axial, i.e. pi-periodical
#' (`TRUE`, the default) or directional, i.e. \eqn{2 \pi}-periodical (`FALSE`).#' @param kappa
#' @param kappa von Mises distribution concentration parameter. Will be
#' estimated using [est.kappa()] if not provided.
#' @param n the number of equally spaced points at which the density is to be estimated.
#'
#' @return numeric
#' @export
#'
#' @references
#' N.I. Fisher (1993) Statistical Analysis of Circular Data, Cambridge University Press.
#'
#' @examples
#' set.seed(1)
#' x <- rvm(10, 0, 100)
#' circular_mode(x, kappa = est.kappa(x))
circular_mode <- function(x, kappa = NULL, axial = TRUE, n = 512) {
  if (is.null(kappa)) kappa <- est.kappa(x, axial = axial)
  dns <- circular_density(x, kappa = kappa, n = n, axial = axial)

  # f <- if (axial) 2 else 1
  # angles <- (c(1:n) / n) * 360 / f
  angles <- seq(0, 360, length.out = n)
  angles[which.max(dns)]
}


#' Circular Summary Statistics
#'
#' Circular mean, standard deviation, variance, quasi-quantiles, mode,
#' 95% confidence angle, standardized skewness and kurtosis
#'
#' @inheritParams circular_mean
#' @param mode logical. Whether the circular mode should be calculated or not.
#' @param kappa  numeric. von Mises distribution concentration parameter used
#' for the circular mode. Will be estimated using [est.kappa()] if not provided.
#' @param fisher.CI logical. Whether Fisher's or the default Mardia/Batchelet's
#' confidence interval should be calculated.
#' @param conf.level numeric. Level of confidence: \eqn{(1 - \alpha \%)/100}.
#' (`0.95` by default).
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
#' circular_summary(sa.por$azi.PoR, w = weighting(san_andreas$unc))
circular_summary <- function(x, w = NULL, axial = TRUE, mode = FALSE, kappa = NULL, fisher.CI = FALSE, conf.level = .95, na.rm = FALSE) {
  if (is.null(w)) {
    w <- rep(1, times = length(x))
  }

  #remove NA
  if (na.rm) {
    keep <- !is.na(x) & !is.na(w)
    x <- x[keep]
    w <- w[keep]
  }

  # order x
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]

  # n <- length(x)

  # confidence interval
  ci_fun <- ifelse(fisher.CI, confidence_interval_fisher, confidence_interval)
  x_CI <- ci_fun(x, conf.level = conf.level, w = w, axial = axial, na.rm = FALSE)

  # circular quantiles
  x_quant <- circular_quantiles(x, w, axial, FALSE)

  # central moments
  x_sk <- second_central_moment(x, w, axial, FALSE)


  res <- c(
    n = length(x),
    mean = circular_mean(x, w, axial, FALSE),
    sd = circular_sd(x, w, axial, FALSE),
    var = circular_var(x, w, axial, FALSE),
    x_quant[1],
    "quasi-median" = unname(x_quant[2]),
    x_quant[3],
    "median" = circular_sample_median(x, axial, FALSE),
    # mode = circular_mode(x, kappa = kappa, axial = axial),
    "CI" = x_CI$conf.angle,
    skewness = x_sk$std_skewness,
    kurtosis = x_sk$std_kurtosi,
    R = mean_resultant_length(ax2dir(x), w = w, FALSE)
  )

  if (mode) {
    if (is.null(kappa)) kappa <- est.kappa(x, w = w, axial = axial)
    mode <- circular_mode(x, kappa = kappa, axial = axial)
    append(res, c("mode" = mode), after = 8)
  } else {
    res
  }
}


#' Orientation tensor
#'
#' @inheritParams circular_mean
#' @param norm logical. Whether the tensor should be normalized.
#'
#' @note \deqn{E = x \cdot x^{T}}
#'
#' @return 2x2 matrix
#' @export
#' @seealso [ot_eigen2d()]
#'
#' @examples
#' test <- rvm(100, mean = 0, k = 10)
#' ortensor2d(test)
ortensor2d <- function(x, w = NULL, norm = FALSE){
  if (is.null(w)) w <- rep(1, times = length(x))

  keep <- !is.na(x) & !is.na(w)
  x <- x[keep]
  w <- w[keep]

  x <- deg2rad(x)
  Z <- if(isTRUE(norm)) 1 else sum(w)

  x <- x[!is.na(x)]
  v <- cbind(w * cos(x), w * sin(x))

  1/Z * (t(v) %*% v)
}

#' Decomposition of orientation tensor
#'
#' Eigenvector decomposition of the orientations.
#'
#' @inheritParams circular_mean
#' @param scale logical. Whether the Eigenvalues should be scaled so they sum up to 1.
#'
#' @return list of Eigenvalues and the angles corresponding to the Eigenvectors.
#' @export
#' @seealso [ortensor2d()]
#'
#' @details
#' Eigenvalues can be interpreted as the fraction of the data explained by the
#' orientation of the associated Eigenvector.
#'
#' @examples
#' test <- rvm(100, mean = 0, k = 10)
#' ot_eigen2d(test, axial = FALSE)
#'
#' data("nuvel1")
#' PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
#' sa.por <- PoR_shmax(san_andreas, PoR, "right")
#' sa_eig <- ot_eigen2d(sa.por$azi.PoR, w = weighting(san_andreas$unc), scale = TRUE)
#' print(sa_eig)
#'
#' rose(sa.por$azi.PoR, muci = FALSE)
#' rose_line(sa_eig$vectors, col = c('red', 'green'),
#'   radius = sa_eig$values, lwd = 2)
#' graphics::legend("topright",
#'   legend = round(sa_eig$values, 2),
#'   col = c('red', 'green'), lty = 1)
ot_eigen2d <- function(x, w = NULL, axial = TRUE, scale = FALSE){
  f <- if (isTRUE(axial)) 2 else 1

  ot <- ortensor2d(f * x, w)
  eig <- eigen(ot)

  ev <- t(eig$vectors)
  ev1 <- atand(ev[1, 2] / ev[1, 1]) / f
  eig$vectors <- c(ev1, ev1+90) %% (360 / f)

  if(isTRUE(scale)) eig$values <- eig$values / sum(eig$values)

  eig
}
