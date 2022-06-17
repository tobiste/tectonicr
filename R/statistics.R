#' Normalized Chi-Squared Test
#'
#' A quantitative comparison between the predicted and observed directions of
#' \eqn{\sigma_{Hmax}}{SHmax} is obtained by the calculation of the average
#' azimuth and by a normalized \eqn{\chi^2}{chi-squared} test.
#'
#' @references Wdowinski, S., 1998, A theory of intraplate
#'   tectonics. *Journal of Geophysical Research: Solid Earth*, v. 103, p.
#'   5037-5059, doi: 10.1029/97JB03390.
#' @inheritParams misfit_shmax
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
#' @export
#' @examples
#' data("nuvel1")
#' euler <- subset(nuvel1, nuvel1$plate.rot == "na") # North America relative to
#' # Pacific plate
#' point <- data.frame(lat = 45, lon = 20)
#' prd <- model_shmax(point, euler)
#' norm_chisq(obs = c(50, 40, 42), prd$sc, unc = c(10, NA, 5))
norm_chisq <- function(obs, prd, unc) {
  stopifnot(is.numeric(obs) & is.numeric(prd) & is.numeric(unc))

  if (length(prd) == 1) {
    prd <- rep(prd, length(obs))
  }

  if (length(unc) == 1) {
    unc <- rep(unc, length(obs))
  }

  if (anyNA(obs)) {
    x <- subset(data.frame(
      obs = obs, prd = prd, unc = unc
    ), !is.na(obs) | !is.na(prd))
    obs <- x[, 1]
    prd <- x[, 2]
    unc <- x[, 3]
    message("NA values have been removed")
  }

  stopifnot(
    length(prd) == length(obs) &
      length(unc) == length(obs) &
      length(unc) == length(prd)
  )

  w <- c()
  x <- c()
  y <- c()
  for (i in seq_along(obs)) {
    if (is.na(obs[i])) {
      x[i] <- NA
      y[i] <- NA
    } else {
      if (!is.na(unc[i]) & unc[i] == 0) {
        unc[i] <- 0.01
      } # uncertainty cannot be 0
      w[i] <- deviation_norm(obs[i] - prd[i])
      x[i] <- (w[i] / unc[i])^2
      y[i] <- (90 / unc[i])^2
    }
  }
  sum(x, na.rm = TRUE) / sum(y, na.rm = TRUE)
}


#' @title Median and statistics on Pi-periodic Data
#'
#' @description Calculate the mean, median, quartile, interquartile range,
#' variance, deviation, and error of orientation data.
#'
#' @param x Numeric vector in degrees.
#' @param quiet logical. If false, a warning message is printed if there are NA
#' values.
#'
#' @return Numeric vector
#'
#' @details Quasi median on the circle, quasi quartiles on a circle, quasi
#' interquartile range on a circle.
#'
#' @source [median()], [quantile()], and [IQR()] are the
#' equivalents for non-periodic data.
#'
#' @references
#' * Ratanaruamkarn, S., Niewiadomska-Bugaj, M., Wang, J.-C. (2009).
#' A New Estimator of a Circular Median. *Communications in Statistics -
#' Simulation and Computation*, 38(6), 1269--1291.
#' \doi{10.1080/03610910902899950}.
#' * Reiter, K., Heidbach, O., Schmitt, D., Haug, K., Ziegler, M., & Moeck, I.
#' (2014). A revised crustal stress orientation database for Canada.
#' *Tectonophysics*, 636, 111-124. \doi{10.1016/j.tecto.2014.08.006}
#'
#' @importFrom stats median
#' @examples
#' x <- c(0, 45, 55, 40 + 180, 50 + 180)
#' circular_mean(x)
#' circular_quasi_median(x)
#' circular_quasi_quantile(x)
#' circular_quasi_IQR(x)
#' circular_var(x)
#' circular_mean_deviation(x)
#' circular_median_deviation(x)
#' circular_mean_error(x)
#'
#' data("san_andreas")
#' circular_quasi_median(san_andreas$azi)
#' @name circle_median
NULL

#' @rdname circle_median
#' @export
circular_quasi_median <- function(x, quiet = TRUE) {
  stopifnot(any(is.numeric(x)))

  if (NA %in% x & quiet) {
    message("NA values have been dropped")
  }

  x <- deg2rad(x) %% pi
  x <- sort(x[!is.na(x)])
  n <- length(x)

  if (n %% 2 != 0) { # if odd
    m <- (n - 1) / 2
    # atand(
    #   sind(x[m+1]) / cosd(x[m+1])
    # ) %% 180
    ss <- sin(x[m + 1])
    cs <- cos(x[m + 1])
  } else { # if even
    m <- n / 2
    # atand(
    #   (sind(x[m]) + sind(x[m + 1])) /
    #     (cosd(x[m]) + cosd(x[m + 1]))
    # ) %% 180
    ss <- sin(x[m]) + sin(x[m + 1])
    cs <- cos(x[m]) + cos(x[m + 1])
  }
  atan2d_spec(ss, cs) %% 180
}

#' @rdname circle_median
#' @export
circular_mean <- function(x) {
  # stopifnot(any(is.numeric(x)))
  #
  # data <- data.frame(x)
  # if (quiet) {
  #   data <- subset(data, !is.na(x))
  # }
  #
  # x <- deg2rad(data$x)
  #
  # s.m <- sum(sin(2 * x))
  # c.m <- sum(cos(2 * x))
  #
  # (atan2d_spec(s.m, c.m) / 2) %% 180
  circular_weighted_mean(x, w = 1, na.rm = TRUE)
}

#' @rdname circle_median
#' @export
circular_quasi_quantile <- function(x, quiet = TRUE) {
  stopifnot(any(is.numeric(x)))

  if (NA %in% x & quiet) {
    message("NA values have been dropped")
  }
  x <- x %% 180
  x <- sort(x[!is.na(x)])
  n <- length(x)

  if (n > 3) {
    # ms <- 1:n

    med <- circular_quasi_median(x)

    if (n %% 4 == 0) {
      m <- n / 4
      lq <- atand(
        sind(x[m + 1]) / cosd(x[m + 1])
      )
      uq <- atand(
        sind(x[3 * m + 1]) / cosd(x[3 * m + 1])
      )
    }
    if (n %% 4 == 1) {
      m <- (n - 1) / 4
      lq <- atand(
        (3 * sind(x[m]) + sind(x[m + 1])) /
          (3 * cosd(x[m]) + cosd(x[m + 1]))
      )
      uq <- atand(
        (3 * sind(x[3 * m]) + sind(x[3 * m + 1])) /
          (3 * cosd(x[3 * m]) + cosd(x[3 * m + 1]))
      )
    }
    if (n %% 4 == 2) {
      m <- (n - 2) / 4
      lq <- atand((sind(x[m]) + sind(x[m + 1])) /
        (cosd(x[m]) + cosd(x[m + 1])))
      uq <- atand((sind(x[3 * m]) + sind(x[3 * m + 1])) /
        (cosd(x[3 * m]) + cosd(x[3 * m + 1])))
    }
    if (n %% 4 == 3) {
      m <- (n - 2) / 4
      lq <- atand((sind(x[m]) + 3 * sind(x[m + 1])) /
        (cosd(x[m]) + 3 * cosd(x[m + 1])))
      uq <- atand((sind(x[3 * m]) +
        3 * sind(x[3 * m + 1])) /
        (cosd(x[3 * m]) +
          3 * cosd(x[3 * m + 1])))
    }


    quantiles <- c(x[1], lq, med, uq, x[length(x)])
    names(quantiles) <- c("0%", "25%", "50%", "75%", "100%")
    return(as.numeric(quantiles))
  } else {
    message("x needs more than 3 values")
    return(NULL)
  }
}

#' @rdname circle_median
#' @export
circular_quasi_IQR <- function(x, quiet = TRUE) {
  quantiles <- circular_quasi_quantile(x)
  deviation_norm(quantiles[4] - quantiles[2])
}

#' @rdname circle_median
#' @export
circular_var <- function(x, quiet = TRUE) {
  stopifnot(any(is.numeric(x)))

  if (NA %in% x & quiet) {
    message("NA values have been dropped")
  }
  x <- x %% 180
  n <- length(x)

  for (i in 1:n) {
    cs <- cosd(x[i])
    ss <- sind(x[i])
  }
  R <- sqrt(sum(cs)^2 + sum(ss)^2)
  1 - R / n
}

#' @rdname circle_median
#' @export
circular_mean_deviation <- function(x, quiet = TRUE) {
  stopifnot(any(is.numeric(x)))

  if (NA %in% x & quiet) {
    message("NA values have been dropped\n")
  }
  x <- x %% 180
  n <- length(x)

  for (i in 1:n) {
    k <- abs(
      180 - abs(x[i] - circular_quasi_median(x))
    )
  }
  180 - (1 / n * sum(k))
}

#' @rdname circle_median
#' @export
circular_median_deviation <- function(x, quiet = TRUE) {
  stopifnot(any(is.numeric(x)))

  if (NA %in% x & quiet) {
    message("NA values have been dropped\n")
  }
  x <- x %% 180

  for (i in seq_along(x)) {
    k <- 180 - abs(180 - abs(x[i] - circular_quasi_median(x)))
  }
  stats::median(k)
}

#' @rdname circle_median
#' @export
circular_mean_error <- function(x, quiet = TRUE) {
  stopifnot(any(is.numeric(x)))

  if (NA %in% x & quiet) {
    message("NA values have been dropped\n")
  }
  x <- x %% 180
  n <- length(x)

  for (i in 1:n) {
    k <- abs(180 - abs(x[i] - circular_quasi_median(x)))
  }
  180 - (1 / n * sum(k))
}

#' @title Weighted Mean and Statistics on Pi-periodic Data
#'
#' @description Calculate the weighted median and standard deviation
#' of orientation data. Weighting is based on the reciprocal of the data
#' uncertainties.
#'
#' @param x Data values. A vector of numeric values in degrees, for which the
#' mean, median or standard deviation are required.
#' @param w Weights. A vector of positive numbers, of the same length as
#' \code{x}.
#' @param na.rm logical value indicating whether \code{NA} values in \code{x}
#' should be stripped before the computation proceeds.
#' @references
#' * Mardia, K.V. (1972). Statistics of Directional Data: Probability and
#' Mathematical Statistics. London: Academic Press.
#' * Ziegler, M. O.; Heidbach O. (2019). Manual of the Matlab Script
#' Stress2Grid v1.1. (WSM Technical Report; 19-02),
#' GFZ German Research Centre for Geosciences. \doi{10.2312/wsm.2019.002}
#' * Heidbach, O., Tingay, M., Barth, A., Reinecker, J., Kurfeß, D., & Müller,
#' B. (2010). Global crustal stress pattern based on the World Stress Map
#' database release 2008. *Tectonophysics* 482, 3–-15,
#' \doi{10.1016/j.tecto.2009.07.023}
#' @examples
#' x <- c(175, 179, 0, 2, 4) + 90
#' unc <- c(5, 1, 0.1, 2, 4)
#' circular_weighted_mean(x, 1 / unc)
#' circular_weighted_var(x, 1 / unc)
#' circular_weighted_sd(x, 1 / unc)
#' circular_weighted_median(x, 1 / unc)
#' circular_weighted_quantiles(x, 1 / unc)
#' circular_weighted_IQR(x, 1 / unc)
#'
#' data("san_andreas")
#' circular_weighted_mean(san_andreas$azi, 1 / san_andreas$unc)
#' @name weighted_circle_stats
NULL

#' @rdname weighted_circle_stats
#' @export
circular_weighted_mean <- function(x, w = NULL, na.rm = TRUE) {
  stopifnot(any(is.numeric(x)))

  if (is.null(w)) {
    w <- rep(1, times = length(x))
  } else {
    w <- as.double(w)
  }
  na.rm <- as.logical(na.rm)
  if (is.na(na.rm)) na.rm <- FALSE

  data <- data.frame(x, w)
  if (na.rm) {
    data <- subset(data, !is.na(x) & !is.na(w))
  }

  x <- deg2rad(data$x) %% pi
  w <- data$w

  Z <- sum(w)

  sin2 <- w * sin(2 * x)
  cos2 <- w * cos(2 * x)
  sumsin2 <- sum(sin2)
  sumcos2 <- sum(cos2)
  meansin2 <- sumsin2 / Z
  meancos2 <- sumcos2 / Z

  meanx_rad <- atan2(meansin2, meancos2) / 2
  rad2deg(meanx_rad) %% 180
}
#' @rdname weighted_circle_stats
#' @export
circular_weighted_var <- function(x, w = NULL, na.rm = TRUE) {
  circular_weighted_sd(x, w, na.rm)^2
}


#' @rdname weighted_circle_stats
#' @export
circular_weighted_sd <- function(x, w = NULL, na.rm = TRUE) {
  stopifnot(any(is.numeric(x)))

  if (is.null(w)) {
    w <- rep(1, times = length(x))
  } else {
    w <- as.double(w)
  }
  na.rm <- as.logical(na.rm)
  if (is.na(na.rm)) na.rm <- FALSE


  data <- data.frame(x, w)
  if (na.rm) {
    data <- subset(data, !is.na(x) & !is.na(w))
  }

  x <- deg2rad(data$x) %% pi
  w <- data$w

  Z <- sum(w)

  sin2 <- w * sin(2 * x)
  cos2 <- w * cos(2 * x)
  sumsin2 <- sum(sin2)
  sumcos2 <- sum(cos2)
  meansin2 <- sumsin2 / Z
  meancos2 <- sumcos2 / Z
  meanR <- sqrt(meancos2^2 + meansin2^2)

  sd <- sqrt(-2 * log(meanR)) / 2
  rad2deg(sd)
}


#' @rdname weighted_circle_stats
#' @export
circular_weighted_median <- function(x, w = NULL, na.rm = TRUE) {
  stopifnot(any(is.numeric(x)))

  if (is.null(w)) {
    w <- rep(1, times = length(x))
  } else {
    w <- as.double(w)
  }
  na.rm <- as.logical(na.rm)
  if (is.na(na.rm)) na.rm <- FALSE

  data <- data.frame(x = x %% 180, w)

  if (na.rm) {
    data <- subset(data, !is.na(x) & !is.na(w))
  }

  data <- dplyr::arrange(data, x)

  x <- deg2rad(data$x)
  w <- data$w

  n <- length(x)

  if (n %% 2 != 0) { # if odd
    m <- (n - 1) / 2
    # atand(
    #   sind(x[m+1]) / cosd(x[m+1])
    # ) %% 180

    sumsin2 <- sin(x[m + 1])
    sumcos2 <- cos(x[m + 1])
  } else { # if even
    m <- n / 2
    # atand(
    #   (sind(x[m]) + sind(x[m + 1])) /
    #     (cosd(x[m]) + cosd(x[m + 1]))
    # ) %% 180

    sumsin2 <- (w[m] * sin(x[m]) + w[m + 1] * sin(x[m + 1])) / (w[m] + w[m + 1])
    sumcos2 <- (w[m] * cos(x[m]) + w[m + 1] * cos(x[m + 1])) / (w[m] + w[m + 1])
  }
  atan2d(sumsin2, sumcos2) %% 180
}



#' @rdname weighted_circle_stats
#' @export
circular_weighted_quantiles <- function(x, w = NULL, na.rm = TRUE) {
  stopifnot(any(is.numeric(x)))

  if (is.null(w)) {
    w <- rep(1, times = length(x))
  } else {
    w <- as.double(w)
  }
  na.rm <- as.logical(na.rm)
  if (is.na(na.rm)) na.rm <- FALSE

  data <- data.frame(x = x %% 180, w)

  if (na.rm) {
    data <- subset(data, !is.na(x) & !is.na(w))
  }

  data <- dplyr::arrange(data, x)

  med <- circular_weighted_median(data$x, data$w)

  x <- deg2rad(data$x)
  w <- data$w

  n <- length(x)

  if (n > 3) {
    if (n %% 4 == 0) {
      m <- n / 4
      sum.sin.lq <- sind(x[m + 1])
      sum.cos.lq <- cosd(x[m + 1])

      sum.sin.uq <- sind(x[3 * m + 1])
      sum.cos.uq <- cosd(x[3 * m + 1])

      Zu <- Zl <- 1
    }
    if (n %% 4 == 1) {
      m <- (n - 1) / 4

      sum.sin.lq <- 3 * w[m] * sind(x[m]) + w[m + 1] * sind(x[m + 1])
      sum.cos.lq <- 3 * w[m] * cosd(x[m]) + w[m + 1] * cosd(x[m + 1])

      sum.sin.uq <- 3 * w[3 * m] * sind(x[3 * m]) + w[3 * m + 1] *
        sind(x[3 * m + 1])
      sum.cos.uq <- 3 * w[3 * m] * cosd(x[3 * m]) + w[3 * m + 1] *
        cosd(x[3 * m + 1])

      Zl <- w[m] + w[m + 1]
      Zu <- w[3 * m] + w[3 * m + 1]
    }
    if (n %% 4 == 2) {
      m <- (n - 2) / 4

      sum.sin.lq <- w[m] * sind(x[m]) + w[m + 1] * sind(x[m + 1])
      sum.cos.lq <- w[m] * cosd(x[m]) + w[m + 1] * cosd(x[m + 1])

      sum.sin.uq <- w[3 * m] * sind(x[3 * m]) + w[3 * m + 1] *
        sind(x[3 * m + 1])
      sum.cos.uq <- w[3 * m] * cosd(x[3 * m]) + w[3 * m + 1] *
        cosd(x[3 * m + 1])

      Zl <- w[m] + w[m + 1]
      Zu <- w[3 * m] + w[3 * m + 1]
    }
    if (n %% 4 == 3) {
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

    lq <- atan2d(mean.sin.lq, mean.cos.lq)
    uq <- atan2d(mean.sin.uq, mean.cos.uq)

    quantiles <- c(
      rad2deg(x[1]), rad2deg(lq), med, rad2deg(uq),
      rad2deg(x[length(x)])
    )
    names(quantiles) <- c("0%", "25%", "50%", "75%", "100%")
    return(as.numeric(quantiles))
  } else {
    message("x needs more than 3 values")
    return(NULL)
  }
}

#' @rdname weighted_circle_stats
#' @export
circular_weighted_IQR <- function(x, w = NULL, na.rm = TRUE) {
  quantiles <- circular_weighted_quantiles(x, w)
  deviation_norm(as.numeric(quantiles[4] - quantiles[2]))
}
