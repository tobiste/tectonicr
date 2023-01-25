nchisq_eq <- function(obs, prd, unc) {
  if (is.na(obs)) {
    x <- NA
    y <- NA
  } else {
    if (!is.na(unc) && unc == 0) {
      unc <- 1
    } # uncertainty cannot be 0
    w <- obs - prd
    x <- (w / unc)^2
    y <- (90 / unc)^2
  }
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
#' @inheritParams deviation_shmax
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
#' @importFrom magrittr %>%
#' @importFrom tidyr drop_na
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
  # stopifnot(is.numeric(obs), is.numeric(prd), is.numeric(unc))
  if (length(prd) == 1) {
    prd <- rep(prd, length(obs))
  }

  if (length(unc) == 1) {
    unc <- rep(unc, length(obs))
  }

  # if (anyNA(obs)) {
  x <- data.frame(
    obs = obs, prd = prd, unc = unc
  ) %>%
    tidyr::drop_na(obs, prd)
  # obs <- x[, 1]
  # prd <- x[, 2]
  # unc <- x[, 3]
  # message("NA values have been removed")
  # }

  xy <- mapply(FUN = nchisq_eq, obs = x[, 1], prd = x[, 2], unc = x[, 3])
  sum(xy[1, ], na.rm = TRUE) / sum(xy[2, ], na.rm = TRUE)
}

mean_SC <- function(x, w, na.rm) {
  stopifnot(any(is.numeric(x)), is.logical(na.rm))

  if (is.null(w)) {
    w <- rep(1, times = length(x))
  } else {
    w <- as.numeric(w)
  }

  data <- data.frame(x, w)
  if (na.rm) {
    data <- tidyr::drop_na(data)
  }

  x <- deg2rad(data$x)
  w <- data$w

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
#' @param x Data values. A vector of numeric values in degrees, for which the
#' mean, median or standard deviation are required.
#' @param w (optional) Weights. A vector of positive numbers, of the same length as
#' \code{x}.
#' @param na.rm logical value indicating whether \code{NA} values in \code{x}
#' should be stripped before the computation proceeds.
#' @param axial logical. Whether the data are axial, i.e. pi-periodical
#' (TRUE, the default) or circular, i.e. 2pi-periodical (FALSE).
#' @importFrom dplyr arrange
#' @importFrom tidyr drop_na
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
#' circular_mean(san_andreas$azi, 1 / san_andreas$unc)
#' circular_var(san_andreas$azi, 1 / san_andreas$unc)
#' circular_quantiles(san_andreas$azi, 1 / san_andreas$unc)
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
  } else {
    w <- as.numeric(w)
  }

  if (axial) {
    f <- 2
    mod <- 180
  } else {
    f <- 1
    mod <- 360
  }
  x <- (x %% mod)
  data <- data.frame(x = x, w)


  if (na.rm) {
    data <- tidyr::drop_na(data)
  }

  data <- dplyr::arrange(data, x)

  x <- f * deg2rad(data$x)
  w <- data$w

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
  } else {
    w <- as.numeric(w)
  }

  if (axial) {
    f <- 2
    mod <- 180
  } else {
    f <- 1
    mod <- 360
  }
  x <- x %% mod

  data <- data.frame(x = x, w)

  if (na.rm) {
    data <- tidyr::drop_na(data)
  }

  data <- dplyr::arrange(data, x)

  x_first <- data$x[1]
  x_last <- data$x[length(data$x)]

  x <- f * deg2rad(data$x)
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
