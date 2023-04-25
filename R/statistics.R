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

mean_resultant <- function(x, w, na.rm) {
  m <- mean_SC(x, w, na.rm)
  R <- sqrt(m[, "C"]^2 + m[, "S"]^2)
  as.numeric(R)
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
#' (TRUE, the default) or circular, i.e. 2pi-periodical (FALSE).
#' @importFrom stats complete.cases
#' @return numeric vector
#' @note Weighting may be the reciprocal of the data uncertainties.
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
  x <- (x %% mod) * f
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
  x <- (x %% mod) * f

  R <- mean_resultant(x = x, w = w, na.rm = na.rm)
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
  x <- (x %% mod) * f

  R <- mean_resultant(x = x, w = w, na.rm = na.rm)
  sd <- sqrt(-2 * log(R)) # / f
  rad2deg(sd + 2 * pi) %% mod
}


#' @rdname circle_stats
#' @export
circular_median <- function(x, w = NULL, axial = TRUE, na.rm = TRUE) {
  meanx <- circular_mean(x, w, axial = TRUE, na.rm)

  if (meanx <= 25 || meanx >= 155) {
    sub <- 90
    x <- x + sub
  } else {
    sub <- 0
  }

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
  x <- (x %% mod)
  data <- cbind(x = x, w = w)


  if (na.rm) {
    data <- data[stats::complete.cases(data), ] # remove NA values
  }

  data <- data[order(data[, "x"]), ]

  x <- f * deg2rad(data[, "x"])
  w <- data[, "w"]

  n <- length(x)

  if (n %% 2 != 0) { # if odd
    m <- (n - 1) / 2
    sumsin2 <- sin(x[m + 1])
    sumcos2 <- cos(x[m + 1])
  } else { # if even
    m <- n / 2
    sumsin2 <- (w[m] * sin(x[m]) + w[m + 1] * sin(x[m + 1])) / (w[m] + w[m + 1])
    sumcos2 <- (w[m] * cos(x[m]) + w[m + 1] * cos(x[m + 1])) / (w[m] + w[m + 1])
  }
  ((atan2d(sumsin2, sumcos2) / f) - sub) %% mod
}



#' @rdname circle_stats
#' @export
circular_quantiles <- function(x, w = NULL, axial = TRUE, na.rm = TRUE) {
  if (circular_mean(x, w, axial = TRUE, na.rm) < 10) {
    sub <- 90
    x <- x + sub
  } else {
    sub <- 0
  }

  med <- circular_median(x, w, axial, na.rm)

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
  x <- x %% mod

  data <- cbind(x = x, w = w)

  if (na.rm) {
    data <- data[stats::complete.cases(data), ] # remove NA values
  }

  data <- data[order(data[, "x"]), ]

  x_first <- data[1, "x"]
  x_last <- data[length(data[, 1]), "x"]

  x <- f * deg2rad(data[, "x"])
  w <- data[, "w"]
  n <- length(x)

  if (n > 3) {
    if (n %% 4 == 0) {
      m <- n / 4
      sum.sin.lq <- sind(x[m + 1])
      sum.cos.lq <- cosd(x[m + 1])

      sum.sin.uq <- sind(x[3 * m + 1])
      sum.cos.uq <- cosd(x[3 * m + 1])

      Zu <- Zl <- 1
    } else if (n %% 4 == 1) {
      m <- (n - 1) / 4

      sum.sin.lq <- 3 * w[m] * sind(x[m]) + w[m + 1] * sind(x[m + 1])
      sum.cos.lq <- 3 * w[m] * cosd(x[m]) + w[m + 1] * cosd(x[m + 1])

      sum.sin.uq <- 3 * w[3 * m] * sind(x[3 * m]) + w[3 * m + 1] *
        sind(x[3 * m + 1])
      sum.cos.uq <- 3 * w[3 * m] * cosd(x[3 * m]) + w[3 * m + 1] *
        cosd(x[3 * m + 1])

      Zl <- w[m] + w[m + 1]
      Zu <- w[3 * m] + w[3 * m + 1]
    } else if (n %% 4 == 2) {
      m <- (n - 2) / 4

      sum.sin.lq <- w[m] * sind(x[m]) + w[m + 1] * sind(x[m + 1])
      sum.cos.lq <- w[m] * cosd(x[m]) + w[m + 1] * cosd(x[m + 1])

      sum.sin.uq <- w[3 * m] * sind(x[3 * m]) + w[3 * m + 1] *
        sind(x[3 * m + 1])
      sum.cos.uq <- w[3 * m] * cosd(x[3 * m]) + w[3 * m + 1] *
        cosd(x[3 * m + 1])

      Zl <- w[m] + w[m + 1]
      Zu <- w[3 * m] + w[3 * m + 1]
    } else { # if (n %% 4 == 3) {
      m <- (n - 2) / 4

      sum.sin.lq <- w[m] * sind(x[m]) + 3 * w[m + 1] * sind(x[m + 1])
      sum.cos.lq <- w[m] * cosd(x[m]) + 3 * w[m + 1] * cosd(x[m + 1])

      sum.sin.uq <- w[3 * m] * sind(x[3 * m]) + 3 * w[3 * m + 1] *
        sind(x[3 * m + 1])
      sum.cos.uq <- w[3 * m] * cosd(x[3 * m]) + 3 * w[3 * m + 1] *
        cosd(x[3 * m + 1])

      Zl <- w[m] + w[m + 1]
      Zu <- w[3 * m] + w[3 * m + 1]
    }
    mean.sin.lq <- sum.sin.lq / Zl
    mean.cos.lq <- sum.cos.lq / Zl

    mean.sin.uq <- sum.sin.uq / Zu
    mean.cos.uq <- sum.cos.uq / Zu

    lq <- atan2d(mean.sin.lq, mean.cos.lq) / f
    uq <- atan2d(mean.sin.uq, mean.cos.uq) / f

    quantiles <- c(
      x_first, rad2deg(lq), med, rad2deg(uq), x_last
    )
    names(quantiles) <- c("0%", "25%", "50%", "75%", "100%")
    return(quantiles - sub)
  } else {
    message("x needs more than 3 values")
    return(NULL)
  }
}

#' @rdname circle_stats
#' @export
circular_IQR <- function(x, w = NULL, axial = TRUE, na.rm = TRUE) {
  quantiles <- circular_quantiles(x, w, axial, na.rm)
  deviation_norm(as.numeric(quantiles[4] - quantiles[2]))
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
