mean_SC <- function(x, w, na.rm) {
  stopifnot(any(is.numeric(x)), is.logical(na.rm))

  w <- if (is.null(w)) {
    rep(1, times = length(x))
  } else {
    as.numeric(w)
  }

  data <- cbind(x = x, w = w)
  if (na.rm) {
    data <- data[stats::complete.cases(data), ] # remove NA values
  }

  x <- deg2rad(data[, "x"])
  w <- data[, "w"]

  Z <- sum(w)

  sin2 <- w * sin(x)
  cos2 <- w * cos(x)
  sumsin2 <- sum(sin2)
  sumcos2 <- sum(cos2)
  meansin2 <- sumsin2 / Z
  meancos2 <- sumcos2 / Z
  cbind(C = meancos2, S = meansin2)
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
  R <- sqrt(m[, "C"]^2 + m[, "S"]^2)
  abs(as.numeric(R))
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
#' * Mardia, K.V. (1972). Statistics of Directional Data: Probability and
#' Mathematical Statistics. London: Academic Press.
#' * Ziegler, M. O.; Heidbach O. (2019). Manual of the Matlab Script
#' Stress2Grid v1.1. *WSM Technical Report* 19-02,
#' GFZ German Research Centre for Geosciences. \doi{10.2312/wsm.2019.002}
#' * Heidbach, O., Tingay, M., Barth, A., Reinecker, J., Kurfess, D., & Mueller,
#' B. (2010). Global crustal stress pattern based on the World Stress Map
#' database release 2008. *Tectonophysics* **482**, 3<U+2013>15,
#' \doi{10.1016/j.tecto.2009.07.023}
#'
#' @examples
#' x <- rvm(10, 0, 100) %% 180
#' unc <- stats::runif(100, 0, 10)
#' circular_mean(x, 1 / unc)
#' circular_var(x, 1 / unc)
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
#'
#' data("nuvel1")
#' PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
#' sa.por <- PoR_shmax(san_andreas, PoR, "right")
#' circular_mean(sa.por$azi.PoR, 1 / san_andreas$unc)
#' circular_median(sa.por$azi.PoR, 1 / san_andreas$unc)
#' circular_var(sa.por$azi.PoR, 1 / san_andreas$unc)
#' circular_quantiles(sa.por$azi.PoR, 1 / san_andreas$unc)
#' @name circle_stats
NULL

#' @rdname circle_stats
#' @export
circular_mean <- function(x, w = NULL, axial = TRUE, na.rm = TRUE) {
  f <- ifelse(axial, 2, 1)
  mod <- 360 / f
  x <- (x * f) %% 360
  m <- mean_SC(x, w, na.rm)
  meanx_rad <- atan2(m[, "S"], m[, "C"]) / f
  meanx_deg <- rad2deg(meanx_rad + 2 * pi) %% mod
  as.numeric(meanx_deg)
}
#' @rdname circle_stats
#' @export
circular_var <- function(x, w = NULL, axial = TRUE, na.rm = TRUE) {
  f <- ifelse(axial, 2, 1)
  mod <- 360 / f
  x <- (x * f) %% 360

  R <- mean_resultant_length(x = x, w = w, na.rm = na.rm)
  1 - R
}

#' @rdname circle_stats
#' @export
circular_sd <- function(x, w = NULL, axial = TRUE, na.rm = TRUE) {
  f <- ifelse(axial, 2, 1)
  mod <- 360 / f
  x <- (x * f) %% 360

  R <- mean_resultant_length(x = x, w = w, na.rm = na.rm)
  sd <- sqrt(-2 * log(R)) # / f
  rad2deg(sd + 2 * pi) %% mod
}


#' @rdname circle_stats
#' @export
circular_median <- function(x, w = NULL, axial = TRUE, na.rm = TRUE) {
  if (is.null(w)) {
    w <- rep(1, times = length(x))
  }

  f <- ifelse(axial, 2, 1)
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
  # med <- circular_median(x, w, axial, na.rm)
  f <- ifelse(axial, 2, 1)
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
  deviation_norm(as.numeric(quantiles[3]), as.numeric(quantiles[1]))
}

#' Circular distance and dispersion
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
#' @references Mardia, K.V. (1972). Statistics of Directional Data: Probability
#' and Mathematical Statistics. London: Academic Press.
#'
#' @importFrom stats complete.cases
#'
#' @returns `circular_distance`returns a numeric vector of positive numbers,
#' `circular_dispersion`returns a positive number.
#'
#' @note
#' If `from` is `NULL`, than the circular variance is returned.
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
  f <- ifelse(axial, 2, 1)

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


cdist2angle <- function(x, axial = TRUE) {
  f <- if (axial) {
    2
  } else {
    1
  }
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
#' Measure of the chance variation expected from sample to sample in estimates of the mean direction.
#' The approximated standard error of the mean direction is computed by the mean
#' resultant length and the MLE concentration parameter \eqn{\kappa}.
#'
#' @inheritParams circular_mean
#'
#' @returns numeric
#'
#' @seealso [mean_resultant_length()], [circular_mean()]
#'
#' @references Davis (1986) Statistics and data analysis in geology. 2nd ed., John Wiley & Sons.
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
    mod <- 90
  } else {
    f <- 1
    mod <- 180
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
  Z_alpha <- z_score(conf.level)
  sde <- circular_sd_error(x, w, axial, na.rm)
  asind(Z_alpha * sde)
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


circular_dispersion_i <- function(x, id, ...) {
  circular_dispersion(x$x[id], y = x$mean[id], w = x$w[id], w.y = x$w.y[id], ...)
}


#' Bootstrapped estimates for circular dispersion
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
  f <- ifelse(axial, 2, 1)
  # mod <- 360 / f
  x <- (x * f) %% 360

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


#' Circular Summary statistics
#'
#' Circular mean, standard deviation, variance, quasi-quantiles, 95% confidence
#' angle, standardized skewness and kurtosis
#'
#' @inheritParams circular_mean
#'
#' @return named vector
#' @export
#'
#' @seealso [circular_mean()], [circular_sd()], [circular_var()],
#' [circular_quantiles()], [confidence_angle()], [second_central_moment()]
#' @examples
#' data("nuvel1")
#' PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
#' sa.por <- PoR_shmax(san_andreas, PoR, "right")
#' circular_summary(sa.por$azi.PoR)
#' circular_summary(sa.por$azi.PoR, w = 1 / san_andreas$unc)
circular_summary <- function(x, w = NULL, axial = TRUE, na.rm = FALSE) {
  f <- ifelse(axial, 2, 1)
  mod <- 360 / f
  x <- (x * f) %% 360

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

  x_mean <- (circular_mean(x, w, F, F) / f) %% mod
  x_sd <- circular_sd(x, w, F, F)
  x_var <- circular_var(x, w, F, F)
  x_CI <- confidence_angle(x, 0.95, w, F, F)
  x_quant <- (circular_quantiles(x, w, F, F) / f) %% mod
  x_sk <- second_central_moment(x, w, F, F)

  setNames(
    c(n, x_mean, x_sd, x_var, x_quant[1], x_quant[2], x_quant[3], x_CI, x_sk$std_skewness, x_sk$std_kurtosis),
    c("n", "mean", "sd", "var", "25%", "median", "75%", "95%CI", "skewness", "kurtosis")
  )
}
