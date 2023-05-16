mean_SC <- function(x, w, na.rm) {
  stopifnot(any(is.numeric(x)), is.logical(na.rm))

  if (is.null(w)) {
    w <- rep(1, times = length(x))
  } else {
    w <- as.numeric(w)
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


#' Mean resultant length
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
#' @export
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

#' @title Summary statistics of directional data
#'
#' @description Calculate the (weighted median) and standard deviation
#' of orientation data.
#'
#' @param x numeric vector. Values in degrees, for which the
#' mean, median or standard deviation are required.
#' @param w (optional) Weights. A vector of positive numbers, of the same length as
#' \code{x}.
#' @param na.rm logical value indicating whether \code{NA} values in \code{x}
#' should be stripped before the computation proceeds.
#' @param axial logical. Whether the data are axial, i.e. pi-periodical
#' (`TRUE`, the default) or directional, i.e. \eqn{2 \pi}-periodical (`FALSE`).
#' @importFrom stats complete.cases
#' @return numeric vector
#' @note Weighting may be the reciprocal of the data uncertainties.
#'
#' Weightings have no effect on median and quantiles if
#' `length(x) %% 2 != 1` and `length(x) %% 4 == 0`, respectively.
#' @references
#' * Mardia, K.V. (1972). Statistics of Directional Data: Probability and
#' Mathematical Statistics. London: Academic Press.
#' * Ziegler, M. O.; Heidbach O. (2019). Manual of the Matlab Script
#' Stress2Grid v1.1. *WSM Technical Report* 19-02,
#' GFZ German Research Centre for Geosciences. \doi{10.2312/wsm.2019.002}
#' * Heidbach, O., Tingay, M., Barth, A., Reinecker, J., Kurfess, D., & Mueller,
#' B. (2010). Global crustal stress pattern based on the World Stress Map
#' database release 2008. *Tectonophysics* **482**, 3–15,
#' \doi{10.1016/j.tecto.2009.07.023}
#' @examples
#' x <- c(175, 179, 0, 2, 4)
#' unc <- c(5, 1, 0.1, 2, 4)
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
#' ep <- subset(nuvel1, nuvel1$plate.rot == "na")
#' sa.por <- PoR_shmax(san_andreas, ep, "right")
#' circular_mean(sa.por$azi.PoR, 1 / san_andreas$unc)
#' circular_median(sa.por$azi.PoR, 1 / san_andreas$unc)
#' circular_var(sa.por$azi.PoR, 1 / san_andreas$unc)
#' circular_quantiles(sa.por$azi.PoR, 1 / san_andreas$unc)
#' @name circle_stats
NULL

#' @rdname circle_stats
#' @export
circular_mean <- function(x, w = NULL, axial = TRUE, na.rm = TRUE) {
  if (axial) {
    f <- 2
    mod <- 180
  } else {
    f <- 1
    mod <- 360
  }
  x <- (x * f) %% 360
  m <- mean_SC(x, w, na.rm)
  meanx_rad <- atan2(m[, "S"], m[, "C"]) / f
  meanx_deg <- rad2deg(meanx_rad + 2 * pi) %% mod
  as.numeric(meanx_deg)
}
#' @rdname circle_stats
#' @export
circular_var <- function(x, w = NULL, axial = TRUE, na.rm = TRUE) {
  if (axial) {
    f <- 2
    mod <- 180
  } else {
    f <- 1
    mod <- 360
  }
  x <- (x * f) %% 360

  R <- mean_resultant_length(x = x, w = w, na.rm = na.rm)
  1 - R
}

#' @rdname circle_stats
#' @export
circular_sd <- function(x, w = NULL, axial = TRUE, na.rm = TRUE) {
  if (axial) {
    f <- 2
    mod <- 180
  } else {
    f <- 1
    mod <- 360
  }
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

  if (axial) {
    f <- 2
    mod <- 180
  } else {
    f <- 1
    mod <- 360
  }
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

  if (axial) {
    f <- 2
    mod <- 180
  } else {
    f <- 1
    mod <- 360
  }
  x <- deg2rad(f * x) %% (2 * pi)

  if (is.null(w)) {
    w <- rep(1, times = length(x))
  }

  data <- cbind(x = x, w = w)
  if (na.rm) {
    data <- data[stats::complete.cases(data), ] # remove NA values
  }
  data <- data[order(data[, "x"]), ]

  # x_first <- data[1, "x"]
  # x_last <- data[length(data[, 1]), "x"]

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

    quantiles <- c(lq, med, uq) / f
    names(quantiles) <- c("25%", "50%", "75%")
    return(quantiles %% mod)
  } else {
    message("x needs more than 3 values")
    return(NULL)
  }
}

#' @rdname circle_stats
#' @export
circular_IQR <- function(x, w = NULL, axial = TRUE, na.rm = TRUE) {
  quantiles <- circular_quantiles(x, w, axial, na.rm)
  deviation_norm(as.numeric(quantiles[3] - quantiles[1]))
}

#' Error of model's prediction
#'
#' The maximum error in the model's predicted azimuth given the Pole of
#' rotations uncertainty and distance of the data point to the pole.
#'
#' @param dist_PoR Distance to Euler pole (great circle distance, in degree)
#' @param sigma_PoR uncertainty of the position of the Pole of rotation
#' (in degree).
#' @references Ramsay, J.A. Folding and fracturing of rocks. McGraw-Hill, New York, 1967.
#' @returns numeric vector. The maximum error for azimuths prediction (in degree)
#' @seealso  [PoR_shmax()] and [model_shmax()] for the model's prediction, and
#' [orthodrome()] for great circle distances.
#' @export
#' @examples
#' prd_err(67, 1)
prd_err <- function(dist_PoR, sigma_PoR = 1) {
  x <- 2 * sind(sigma_PoR)^2
  y <- 1 + cosd(dist_PoR)
  acos_beta <- sqrt(1 - x / (sind(dist_PoR)^2) * y)
  acosd(acos_beta) / 2
}



#' Apply Rolling Functions using circular statistics
#'
#' A generic function for applying a function to rolling margins of an array.
#'
#' @inheritParams circular_mean
#' @param width numeric vector or list. In the simplest case this is an integer
#' specifying the window width (in numbers of observations) which is aligned to
#' the original sample according to the `align` argument. Alternatively, width
#' can be a list regarded as offsets compared to the current time.
#' @param FUN the function to be applied
#' @param by.column logical. If `TRUE`, FUN is applied to each column separately.
#' @param fill a three-component vector or list (recycled otherwise) providing
#' filling values at the left/within/to the right of the data range. See the
#' fill argument of [zoo::na.fill()] for details
#' @param partial logical or numeric. If `FALSE` then `FUN` is only
#' applied when all indexes of the rolling window are within the observed time
#' range. If `TRUE` (default), then the subset of indexes that are in range
#' are passed to `FUN`. A numeric argument to partial can be used to determine
#' the minimal window size for partial computations. See below for more details.
#' @inheritDotParams zoo::rollapply -data
#' @returns numeric vector  with the results of the rolling function.
#' @note If the rolling statistics are applied to values that are a function of
#' distance it is recommended to sort the values first.
#' @importFrom zoo rollapply
#' @returns numeric vector
#' @export
#' @examples
#' data("plates")
#' plate_boundary <- subset(plates, plates$pair == "na-pa")
#' data("san_andreas")
#' ep <- subset(nuvel1, nuvel1$plate.rot == "na")
#' distance <- distance_from_pb(
#'   x = san_andreas,
#'   euler = ep,
#'   pb = plate_boundary,
#'   tangential = TRUE
#' )
#' dat <- san_andreas[order(distance), ]
#' roll_circstats(dat$azi, w = 1 / dat$unc, circular_mean, width = 51)
roll_circstats <- function(x, w = NULL,
                           FUN,
                           axial = TRUE, na.rm = TRUE,
                           width, by.column = FALSE,
                           partial = TRUE,
                           fill = NA,
                           ...) {
  FUN <- match.fun(FUN)

  zoo::rollapply(
    cbind(x, w),
    width = width,
    FUN = function(x) {
      FUN(x[, 1], w = x[, 2], axial, na.rm)
    },
    by.column = by.column,
    partial = partial,
    fill = fill,
    ...
  )
}

#' Apply Rolling Functions using circular statistics
#'
#' A generic function for applying a function to rolling margins of an array.
#'
#' @inheritParams norm_chisq
#' @param width numeric vector or list. In the simplest case this is an integer
#' specifying the window width (in numbers of observations) which is aligned to
#' the original sample according to the `align` argument. Alternatively, width
#' can be a list regarded as offsets compared to the current time.
#' @param by.column logical. If `TRUE`, FUN is applied to each column separately.
#' @param fill a three-component vector or list (recycled otherwise) providing
#' filling values at the left/within/to the right of the data range. See the
#' fill argument of [zoo::na.fill()] for details
#' @param partial logical or numeric. If `FALSE` then `FUN` is only
#' applied when all indexes of the rolling window are within the observed time
#' range. If `TRUE` (default), then the subset of indexes that are in range
#' are passed to `FUN`. A numeric argument to partial can be used to determine
#' the minimal window size for partial computations. See below for more details.
#' @inheritDotParams zoo::rollapply -data
#' @returns numeric vector  with the results of the rolling function.
#' @note If the rolling statistics are applied to values that are a function of
#' distance it is recommended to sort the values first.
#' @importFrom zoo rollapply
#' @export
#' @returns numeric vector
#' @examples
#' data("plates")
#' plate_boundary <- subset(plates, plates$pair == "na-pa")
#' data("san_andreas")
#' ep <- subset(nuvel1, nuvel1$plate.rot == "na")
#' distance <- distance_from_pb(
#'   x = san_andreas,
#'   euler = ep,
#'   pb = plate_boundary,
#'   tangential = TRUE
#' )
#' dat <- san_andreas[order(distance), ]
#' dat.PoR <- PoR_shmax(san_andreas, ep, "right")
#' roll_normchisq(dat.PoR$azi.PoR, 135, dat$unc, width = 51)
roll_normchisq <- function(obs, prd, unc = NULL,
                           width, by.column = FALSE,
                           partial = TRUE,
                           fill = NA,
                           ...) {
  zoo::rollapply(
    cbind(obs, prd, unc),
    width = width,
    FUN = function(x) {
      norm_chisq(x[, 1], x[, 2], x[, 3])
    },
    by.column = by.column,
    partial = partial,
    fill = fill,
    ...
  )
}


z_score <- function(conf.level) {
  stats::qnorm(1 - (1 - conf.level) / 2)
}

#' Standard error of mean direction
#'
#' Measure of the chance variation expected from sample to sample in estimates of the mean direction.
#' The approximated standard error of the mean direction is computed by the mean
#' resultant length and the MLE concentration parameter \eqn{\kappa}.
#'
#' @inheritParams circular_mean
#' @returns Angle in degrees
#' @seealso [mean_resultant_length()], [circular_mean()]
#' @references Davis (1986) Statistics and data analysis in geology. 2nd ed., John Wiley & Sons.
#' @export
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
#' ep <- subset(nuvel1, nuvel1$plate.rot == "na")
#' sa.por <- PoR_shmax(san_andreas, ep, "right")
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

  # n <- length(x)
  n <- sum(w)

  kappa <- est.kappa(x, w = w, axial = axial, na.rm = na.rm)

  x <- (x * f) %% 360
  R <- mean_resultant_length(x, w = w, na.rm = na.rm)

  sde <- 1 / sqrt(n * R * kappa)
  rad2deg(sde + 2 * pi) %% mod
}

#' Confidence interval around the mean direction
#'
#' Probabilistic limit on the location of the true or population mean direction,
#' assuming that the estimation errors are normally distributed.
#'
#' @inheritParams circular_mean
#' @param conf.level Level of confidence: \eqn{(1 - \alpha \%)/100}.
#' (`0.95` by default).
#' @returns Angle in degrees
#' @seealso [mean_resultant_length()], [circular_sd_error()]
#' @references Davis (1986) Statistics and data analysis in geology. 2nd ed., John Wiley & Sons.
#'
#' Jammalamadaka, S. Rao and Sengupta, A. (2001). Topics in Circular Statistics, Sections 3.3.3 and 3.4.1, World Scientific Press, Singapore.
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
#' ep <- subset(nuvel1, nuvel1$plate.rot == "na")
#' sa.por <- PoR_shmax(san_andreas, ep, "right")
#' confidence_angle(sa.por$azi.PoR, w = 1 / san_andreas$unc)
#' confidence_interval(sa.por$azi.PoR, w = 1 / san_andreas$unc)
#' @name confidence
NULL

#' @rdname confidence
#' @export
confidence_angle <- function(x, conf.level = .95, w = NULL, axial = TRUE, na.rm = TRUE) {
  (circular_sd_error(x, w, axial, na.rm) * z_score(conf.level)) %% 180
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

# Tests ####
nchisq_eq <- function(obs, prd, unc) {
  # if (is.na(obs)) {
  #   x <- y <- NA
  # } else {
  if (is.na(unc) || unc == 0) unc <- 1 # uncertainty cannot be 0
  w <- deviation_norm(obs - prd)
  x <- (w / unc)^2
  y <- (90 / unc)^2
  # }
  return(c(x, y))
}

#' Normalized Chi-Squared Test
#'
#' A quantitative comparison between the predicted and observed directions of
#' \eqn{\sigma_{Hmax}}{SHmax} is obtained by the calculation of the average
#' azimuth and by a normalized \eqn{\chi^2}{chi-squared} test.
#'
#' @references Wdowinski, S., 1998, A theory of intraplate
#'   tectonics. *Journal of Geophysical Research: Solid Earth*, **103**,
#'   5037-5059, doi: 10.1029/97JB03390.
#' @param prd Numeric vector containing the modeled azimuths of
#' \eqn{\sigma_{Hmax}}{SHmax}, i.e.
#' the return object from \code{model_shmax()}
#' @param obs Numeric vector containing the observed azimuth of
#' \eqn{\sigma_{Hmax}}{SHmax},
#' same length as \code{prd}
#' @param unc Uncertainty of observed \eqn{\sigma_{Hmax}}{SHmax}, either a
#' numeric vector or a number
#' @return Numeric vector
#' @details
#' The normalized \eqn{\chi^2}{chi-squared} test is
#' \deqn{ {Norm} \chi^2_i =
#'  = \frac{
#'    \sum^M_{i = 1} \left( \frac{\alpha_i - \alpha_{{predict}}}{\sigma_i}
#'    \right) ^2}
#'    {\sum^M_{i = 1} \left( \frac{90}{\sigma_i} \right) ^2 }}{
#'    (sum( ((obs-prd)/unc)^2 ) / sum( (90/unc)^2 )
#'    }
#' The value of the chi-squared test statistic is a number between 0 and 1
#' indicating the quality of the predicted \eqn{\sigma_{Hmax}}{SHmax}
#' directions. Low values
#' (\eqn{\le 0.15}) indicate good agreement,
#' high values (\eqn{> 0.7}) indicate a systematic misfit between predicted and
#' observed \eqn{\sigma_{Hmax}}{SHmax} directions.
#' @importFrom stats complete.cases
#' @export
#' @examples
#' data("nuvel1")
#' euler <- subset(nuvel1, nuvel1$plate.rot == "na") # North America relative to
#' # Pacific plate
#' data(san_andreas)
#' point <- data.frame(lat = 45, lon = 20)
#' prd <- model_shmax(point, euler)
#' norm_chisq(obs = c(50, 40, 42), prd$sc, unc = c(10, NA, 5))
#'
#' data(san_andreas)
#' prd2 <- PoR_shmax(san_andreas, euler, type = "right")
#' norm_chisq(obs = prd2$azi.PoR, 135, unc = san_andreas$unc)
norm_chisq <- function(obs, prd, unc) {
  N <- length(obs)
  if (length(prd) == 1) {
    prd <- rep(prd, N)
  }

  if (length(unc) == 1) {
    unc <- rep(unc, N)
  }
  stopifnot(length(prd) == N, length(unc) == N)

  x <- cbind(obs, prd, unc)
  x_compl <- matrix(
    x[stats::complete.cases(x[, 1]) & stats::complete.cases(x[, 2]), ],
    ncol = 3
  ) # remove NA values
  # stopifnot(length(x) > 0)

  xy <- mapply(
    FUN = nchisq_eq,
    obs = x_compl[, 1], prd = x_compl[, 2], unc = x_compl[, 3]
  )
  sum(xy[1, ], na.rm = TRUE) / sum(xy[2, ], na.rm = TRUE)
}


#' Rayleigh Test of Circular Uniformity
#'
#' Performs a Rayleigh test of uniformity (or randomness), assessing the
#' significance of the mean resultant length.
#' The alternative hypothesis is an unimodal distribution with unknown mean
#' direction and unknown mean resultant length if `mu` is `NULL`.
#' If `mu` is specified the alternative hypothesis is a unimodal distribution with a
#' specified mean direction and unknown mean resultant length.
#'
#' @param x numeric vector. Values in degrees
#' @param axial logical. Whether the data are axial, i.e. \eqn{\pi}-periodical
#' (`TRUE`, the default) or directional, i.e. \eqn{2 \pi}-periodical (`FALSE`).
#' @param mu (optional) The specified or known mean direction (in degrees) in alternative hypothesis
#' @details
#' If `statistic > p.value`, the null hypothesis is rejected,
#' i.e. the length of the mean resultant differs significantly from zero.
#' If not, randomness (uniform distribution) cannot be excluded.
#'
#' @note Although the Rayleigh test is consistent against (non-uniform)
#' von Mises alternatives, it is not consistent against alternatives with `p = 0`
#' (in particular, distributions with antipodal symmetry, i.e. axial data).
#' Tests of non-uniformity which are consistent against all alternatives
#' include Kuiper’s test ([kuiper_test()]) and Watson’s \eqn{U^2} test
#' ([watson_test()]).
#' @returns a list with the components:
#' \describe{
#'  \item{`statistic`}{mean resultant length}
#'  \item{`p.value`}{significance level of the test statistic}
#'  \item{`p.value2`}{modified significance level (Cordeiro and Ferrari, 1991)}
#' }
#' @references
#' Mardia and Jupp (2000). Directional Statistics. John Wiley and Sons.
#'
#' Wilkie (1983): Rayleigh Test for Randomness of Circular Data. Appl. Statist. 32, No. 3, pp. 311-312
#'
#' Jammalamadaka, S. Rao and Sengupta, A. (2001). Topics in Circular Statistics, Sections 3.3.3 and 3.4.1, World Scientific Press, Singapore.
#' @seealso [mean_resultant_length()], [circular_mean()], [norm_chisq()], [kuiper_test()], [watson_test()]
#' @export
#' @examples
#' # Example data from Mardia and Jupp (2001), pp. 93
#' pidgeon_homing <- c(55, 60, 65, 95, 100, 110, 260, 275, 285, 295)
#' rayleigh_test(pidgeon_homing, axial = FALSE)
#'
#' # Example data from Davis (1986), pp. 316
#' finland_stria <- c(
#'   23, 27, 53, 58, 64, 83, 85, 88, 93, 99, 100, 105, 113,
#'   113, 114, 117, 121, 123, 125, 126, 126, 126, 127, 127, 128, 128, 129, 132,
#'   132, 132, 134, 135, 137, 144, 145, 145, 146, 153, 155, 155, 155, 157, 163,
#'   165, 171, 172, 179, 181, 186, 190, 212
#' )
#' rayleigh_test(finland_stria, axial = FALSE)
#' rayleigh_test(finland_stria, mu = 105, axial = FALSE)
#'
#' # Example data from Mardia and Jupp (2001), pp. 99
#' atomic_weight <- c(
#'   rep(0, 12), rep(3.6, 1), rep(36, 6), rep(72, 1),
#'   rep(108, 2), rep(169.2, 1), rep(324, 1)
#' )
#' rayleigh_test(atomic_weight, 0, axial = FALSE)
#'
#' data(san_andreas)
#' rayleigh_test(san_andreas$azi)
#' data("nuvel1")
#' ep <- subset(nuvel1, nuvel1$plate.rot == "na")
#' sa.por <- PoR_shmax(san_andreas, ep, "right")
#' rayleigh_test(sa.por$azi.PoR, mu = 135)
rayleigh_test <- function(x, mu = NULL, axial = TRUE) {
  if (axial) {
    f <- 2
  } else {
    f <- 1
  }
  x <- (na.omit(x) * f) %% 360
  n <- length(x)

  if (is.null(mu)) {
    R <- mean_resultant_length(x, na.rm = FALSE)
    S <- 2 * n * R^2
    S2 <- (1 - 1 / (2 * n)) * S + (n * R^4) / 2
    # if(n <= 10){
    #  p.value <- p_value3(R, n)
    # } else  {
    p.value <- rayleigh_p_value1(S / 2, n)
    # }
    p.value2 <- rayleigh_p_value1(S2 / 2, n)

    result <- list(
      statistic = R,
      p.value = p.value,
      p.value2 = p.value2
    )
    if (R > p.value2) {
      message("Reject Null Hypothesis\n")
    } else {
      message("Do Not Reject Null Hypothesis\n")
    }
  } else {
    mu <- (f * mu) %% 360
    C <- (sum(cosd(x - mu))) / n
    s <- sqrt(2 * n) * C
    p.value <- rayleigh_p_value2(s, n)

    result <- list(
      statistic = C,
      p.value = p.value
    )
    if (C > p.value) {
      message("Reject Null Hypothesis\n")
    } else {
      message("Do Not Reject Null Hypothesis\n")
    }
  }

  return(result)
}

rayleigh_p_value1 <- function(K, n) {
  # Pearson. 1906; Greenwood and Durand, 1955
  P <- exp(-K)
  if (n < 50) {
    temp <- 1 +
      (2 * K - K^2) / (4 * n) -
      (24 * K - 132 * K^2 + 76 * K^3 - 9 * K^4) / (288 * n^2)
  } else {
    temp <- 1
  }
  P * temp
  min(max(P * temp, 0), 1)
}

# rayleigh_p_value3 <- function(R, n) {
#   # Wilkie 1983
#   Rn <- R * n
#   temp <- sqrt(1 + 4 * n + 4 * (n^2 - Rn^2)) - (1 + 2 * n)
#   round(exp(temp), 3)
# }

rayleigh_p_value2 <- function(K, n) {
  # Greenwood and Durand, 1957
  pK <- stats::pnorm(K) # distribution function of standard normal distribution
  fK <- stats::dnorm(K) # density function of standard normal distribution
  P <- 1 - pK + fK * (
    (3 * K - K^3) / (16 * n) +
      (15 * K + 305 * K^3 - 125 * K^5 + 9 * K^7) / (4608 * n^2)
  )
  min(max(P, 0), 1)
}

#' Normalized Rayleigh goodness-of-fit test
#'
#' Performs a Rayleigh test of uniformity using observations weighted by
#' their uncertainties.
#'
#' @param x numeric vector. Values in degrees
#' @param unc numeric. The standard deviations of `x`. If `NULL`, the non-weighted
#' Rayleigh test is performed.
#' @param prd The specified or known mean direction (in degrees) in alternative hypothesis
#' @param axial logical. Whether the data are axial, i.e. \eqn{\pi}-periodical
#' (`TRUE`, the default) or directional, i.e. \eqn{2 \pi}-periodical (`FALSE`).
#'
#' @details
#' The Null hypothesis is that the distributions of the observations `x` is an
#' unimodal distribution with a specified mean direction (`prd`).
#' If `statistic > p.value`, the null hypothesis is rejected.
#' If not, randomness (uniform distribution) cannot be excluded.
#' @returns a list with the components:
#' \describe{
#'  \item{`statistic`}{Test statistic}
#'  \item{`p.value`}{significance level of the test statistic}
#' }
#' @seealso [rayleigh_test()]
#' @export
#'
#' @examples
#' data("cpm_models")
#' data(san_andreas)
#' ep <- equivalent_rotation(cpm_models %>% filter(model == "NNR-MORVEL56"), "na", "pa")
#' sa.por <- PoR_shmax(san_andreas, ep, "right")
#' data("iceland")
#' ep.ice <- equivalent_rotation(cpm_models %>% filter(model == "NNR-MORVEL56"), "eu", "na")
#' ice.por <- PoR_shmax(iceland, ep.ice, "out")
#' data("tibet")
#' ep.tib <- equivalent_rotation(cpm_models %>% filter(model == "NNR-MORVEL56"), "eu", "in")
#' tibet.por <- PoR_shmax(tibet, ep.tib, "in")
#'
#' norm_rayleigh_test(tibet.por$azi.PoR, unc = tibet$unc, prd = 90)
#' norm_rayleigh_test(ice.por$azi.PoR, unc = iceland$unc, prd = 0)
#' norm_rayleigh_test(sa.por$azi.PoR, unc = san_andreas$unc, prd = 135)
norm_rayleigh_test <- function(x, unc, prd = NULL, axial = TRUE) {
  if (is.null(unc)) {
    rayleigh_test(x, mu = prd, axial = axial)
  } else {
    if (axial) {
      f <- 2
    } else {
      f <- 1
    }

    data <- cbind(x = x, unc = unc)
    data <- data[stats::complete.cases(data), ] # remove NA values

    x <- deg2rad(data[, "x"])
    prd <- deg2rad(prd)
    unc <- data[, "unc"]
    w <- 1 / unc
    Z <- sum(w)

    a <- w * (1 - cos(f * (x - prd)))
    b <- 1 - (w / 90)

    Csq <- sum(a * a) / sum(b * b)
    C <- sqrt(Csq)

    s <- sqrt(2 * Z) * C
    p.value <- rayleigh_p_value2(s, Z)

    result <- list(
      statistic = C,
      p.value = p.value
    )
    if (C > p.value) {
      message("Reject Null Hypothesis\n")
    } else {
      message("Do Not Reject Null Hypothesis\n")
    }
    return(result)
  }
}

#' Kuiper Test of Circular Uniformity
#'
#' Kuiper test for circular random distribution.
#'
#' @param x numeric vector containing the circular data which are expressed in degrees
#' @param alpha Significance level of the test. Valid levels are `0.01`, `0.05`, and `0.1`.
#' This argument may be omitted (`NULL`, the default), in which case, a range for the p-value will be returned.
#' @param axial logical. Whether the data are axial, i.e. \eqn{\pi}-periodical
#' (`TRUE`, the default) or circular, i.e. \eqn{2 \pi}-periodical (`FALSE`).
#' @returns list containing the test statistic `statistic` and the significance
#' level `p.value`.
#' @details
#' Kuiper's test statistic is a rotation-invariant Kolmogorov-type test statistic.
#' The critical values of a modified Kuiper's test statistic are used according
#' to the tabulation given in Stephens (1970).
#'
#' If `statistic > p.value`, the null hypothesis is rejected.
#' If not, randomness (uniform distribution) cannot be excluded.
#' @export
#' @examples
#' # Example data from Mardia and Jupp (2001), pp. 93
#' pidgeon_homing <- c(55, 60, 65, 95, 100, 110, 260, 275, 285, 295)
#' kuiper_test(pidgeon_homing, alpha = .05)
#'
#' data(san_andreas)
#' data("nuvel1")
#' ep <- subset(nuvel1, nuvel1$plate.rot == "na")
#' sa.por <- PoR_shmax(san_andreas, ep, "right")
#' kuiper_test(sa.por$azi.PoR, alpha = .05)
kuiper_test <- function(x, alpha = 0, axial = TRUE) {
  if (!any(c(0, 0.01, 0.025, 0.05, 0.1, 0.15) == alpha)) {
    stop("'alpha' must be one of the following values: 0, 0.01, 0.025, 0.05, 0.1, 0.15")
  }
  kuiper.crits <- cbind(
    c(0.15, 0.1, 0.05, 0.025, 0.01),
    c(1.537, 1.62, 1.747, 1.862, 2.001)
  )
  if (axial) {
    f <- 2
  } else {
    f <- 1
  }
  x <- (na.omit(x) * f) %% 360
  u <- sort(deg2rad(x) %% (2 * pi)) / (2 * pi)
  n <- length(x)
  i <- 1:n
  D.P <- max(i / n - u)
  D.M <- max(u - (i - 1) / n)
  sqrt_n <- sqrt(n)
  V <- D.P + D.M
  V <- V * (sqrt_n + 0.155 + 0.24 / sqrt_n)

  if (alpha == 0) {
    if (V < 1.537) {
      p.value <- "P-value > 0.15"
    } else if (V < 1.62) {
      p.value <- "0.10 < P-value < 0.15"
    } else if (V < 1.747) {
      p.value <- "0.05 < P-value < 0.10"
    } else if (V < 1.862) {
      p.value <- "0.025 < P-value < 0.05"
    } else if (V < 2.001) {
      p.value <- "0.01 < P-value < 0.025"
    } else {
      p.value <- "P-value < 0.01"
    }
  } else {
    p.value <- kuiper.crits[(1:5)[alpha == c(kuiper.crits[, 1])], 2]

    if (V > p.value) {
      message("Reject Null Hypothesis\n")
    } else {
      message("Do Not Reject Null Hypothesis\n")
    }
  }
  return(
    list(
      statistic = V,
      p.value = p.value
    )
  )
}

#' Watson's \eqn{U^2} Test of Circular Uniformity
#'
#' Watson's test for circular random distribution.
#'
#' @param x numeric vector. Values in degrees
#' @param alpha Significance level of the test. Valid levels are `0.01`, `0.05`, and `0.1`.
#' This argument may be omitted (`NULL`, the default), in which case, a range for the p-value will be returned.
#' @param axial logical. Whether the data are axial, i.e. \eqn{\pi}-periodical
#' (`TRUE`, the default) or circular, i.e. \eqn{2 \pi}-periodical (`FALSE`).
#' @param dist Distribution to test for. The default, `"uniform"`, is the
#' uniform distribution. `"vonmises"` tests the von Mises distribution.
#' @returns list containing the test statistic `statistic` and the significance
#' level `p.value`.
#' @details
#' Watson's test statistic is a rotation-invariant Cramer - von Mises test statistic.
#' If `statistic > p.value`, the null hypothesis is rejected.
#' If not, randomness (uniform distribution) cannot be excluded.
#' @references Mardia and Jupp (2000). Directional Statistics. John Wiley and Sons.
#' @export
#' @examples
#' # Example data from Mardia and Jupp (2001), pp. 93
#' pidgeon_homing <- c(55, 60, 65, 95, 100, 110, 260, 275, 285, 295)
#' watson_test(pidgeon_homing, alpha = .05)
#'
#' data(san_andreas)
#' data("nuvel1")
#' ep <- subset(nuvel1, nuvel1$plate.rot == "na")
#' sa.por <- PoR_shmax(san_andreas, ep, "right")
#' watson_test(sa.por$azi.PoR, alpha = .05)
#' watson_test(sa.por$azi.PoR, alpha = .05, dist = "vonmises")
watson_test <- function(x, alpha = 0, dist = c("uniform", "vonmises"), axial = TRUE) {
  if (!any(c(0, 0.01, 0.025, 0.05, 0.1) == alpha)) {
    stop("'alpha' must be one of the following values: 0, 0.01, 0.025, 0.05, 0.1")
  }

  dist <- match.arg(dist)
  x <- na.omit(x)
  n <- length(x)

  if (dist == "uniform") {
    if (axial) {
      f <- 2
    } else {
      f <- 1
    }
    x <- (x * f) %% 360

    # U2 Statistic:
    u <- sort(deg2rad(x)) / (2 * pi)
    u.bar <- mean(u)
    i <- 1:n
    u2 <- sum((u - u.bar - (i - .5) / n + .5)^2) + 1 / (12 * n)
    statistic <- (u2 - 0.1 / n + 0.1 / (n^2)) * (1 + 0.8 / n)

    # P-value:
    crits <- c(99, 0.267, 0.221, 0.187, 0.152)
    if (n < 8) {
      p.value <- NA
      warning("Total Sample Size < 8:  Results are not valid")
    }
    if (alpha == 0) {
      if (statistic > 0.267) {
        p.value <- "P-value < 0.01"
      } else if (statistic > 0.221) {
        p.value <- "0.01 < P-value < 0.025"
      } else if (statistic > 0.187) {
        p.value <- "0.025 < P-value < 0.05"
      } else if (statistic > 0.152) {
        p.value <- "0.05 < P-value < 0.10"
      } else {
        p.value <- "P-value > 0.10"
      }
    } else {
      index <- (1:5)[alpha == c(0, 0.01, 0.025, 0.05, 0.1)]
      p.value <- crits[index]

      if (statistic > p.value) {
        message("Reject Null Hypothesis\n")
      } else {
        message("Do Not Reject Null Hypothesis\n")
      }
    }
  } else {
    u2.crits <- cbind(
      c(0, 0.5, 1, 1.5, 2, 4, 100),
      c(0.052, 0.056, 0.066, 0.077, 0.084, 0.093, 0.096),
      c(0.061, 0.066, 0.079, 0.092, 0.101, 0.113, 0.117),
      c(0.081, 0.09, 0.11, 0.128, 0.142, 0.158, 0.164)
    )

    mu <- circular_mean(x, axial = axial, na.rm = FALSE)
    kappa.mle <- est.kappa(x, axial = axial)
    x <- deg2rad(x - mu)
    x <- matrix(x, ncol = 1)
    z <- apply(x, 1, pvm, 0, kappa.mle)
    z <- sort(z)
    z.bar <- mean(z)
    i <- c(1:n)
    sum.terms <- (z - (2 * i - 1) / (2 * n))^2
    statistic <- sum(sum.terms) - n * (z.bar - 0.5)^2 + 1 / (12 * n)
    if (kappa.mle < 0.25) {
      row <- 1
    } else if (kappa.mle < 0.75) {
      row <- 2
    } else if (kappa.mle < 1.25) {
      row <- 3
    } else if (kappa.mle < 1.75) {
      row <- 4
    } else if (kappa.mle < 3) {
      row <- 5
    } else if (kappa.mle < 5) {
      row <- 6
    } else {
      row <- 7
    }
    if (alpha != 0) {
      if (alpha == 0.1) {
        col <- 2
      } else if (alpha == 0.05) {
        col <- 3
      } else if (alpha == 0.01) {
        col <- 4
      } else {
        stop("Invalid input for alpha", "\n", "\n")
      }
      p.value <- u2.crits[row, col]
      if (statistic > p.value) {
        message("Reject Null Hypothesis\n")
      } else {
        message("Do Not Reject Null Hypothesis\n")
      }
    } else {
      if (statistic < u2.crits[row, 2]) {
        p.value <- "P-value > 0.10"
      } else if ((statistic >= u2.crits[row, 2]) && (statistic < u2.crits[row, 3])) {
        p.value <- "0.05 < P-value > 0.10"
      } else if ((statistic >= u2.crits[row, 3]) && (statistic < u2.crits[row, 4])) {
        p.value <- "0.01 < P-value > 0.05"
      } else {
        p.value <- "P-value < 0.01"
      }
    }
  }
  list(
    statistic = statistic,
    p.value = p.value
  )
}


# Distribution ####
pvm.mu0 <- function(theta, kappa, acc) {
  flag <- TRUE
  p <- 1
  sum <- 0
  while (flag) {
    term <- (besselI(x = kappa, nu = p, expon.scaled = FALSE) *
      sin(p * theta)) / p
    sum <- sum + term
    p <- p + 1
    if (abs(term) < acc) {
      flag <- FALSE
    }
  }
  theta / (2 * pi) + sum / (pi * besselI(
    x = kappa, nu = 0,
    expon.scaled = FALSE
  ))
}

pvm <- function(theta, mu, kappa, acc = 1e-20) {
  theta <- theta %% (2 * pi)
  mu <- mu %% (2 * pi)

  if (mu == 0) {
    pvm.mu0(theta, kappa, acc)
  } else {
    if (theta <= mu) {
      upper <- (theta - mu) %% (2 * pi)
      if (upper == 0) {
        upper <- 2 * pi
      }
      lower <- (-mu) %% (2 * pi)
      pvm.mu0(upper, kappa, acc) - pvm.mu0(lower, kappa, acc)
    } else {
      upper <- theta - mu
      lower <- mu %% (2 * pi)
      pvm.mu0(upper, kappa, acc) + pvm.mu0(lower, kappa, acc)
    }
  }
}

A1inv <- function(x) {
  ifelse(0 <= x & x < 0.53, 2 * x + x^3 + (5 * x^5) / 6,
    ifelse(x < 0.85, -0.4 + 1.39 * x + 0.43 / (1 - x), 1 / (x^3 - 4 * x^2 + 3 * x))
  )
}

est.kappa <- function(x, bias = FALSE, ...) {
  x <- na.omit(x)
  mean.dir <- circular_mean(x, ...)
  kappa <- abs(A1inv(mean.default(cosd(x - mean.dir))))
  if (bias) {
    kappa.ml <- kappa
    n <- length(x)
    if (kappa.ml < 2) {
      kappa <- max(kappa.ml - 2 * (n * kappa.ml)^-1, 0)
    }
    if (kappa.ml >= 2) {
      kappa <- ((n - 1)^3 * kappa.ml) / (n^3 + n)
    }
  }
  kappa
}
