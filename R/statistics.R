#' Normalized Chi-Squared Test
#'
#' A quantitative comparison between the predicted and observed directions of
#' \eqn{\sigma_\text{Hmax}}{SHmax} is obtained by the calculation of the average
#' azimuth and by a normalized \eqn{\chi^2}{chi-squared} test.
#'
#' @references Wdowinski, S., 1998, A theory of intraplate
#'   tectonics. *Journal of Geophysical Research: Solid Earth*, v. 103, p.
#'   5037--5059, \doi{10.1029/97JB03390}.
#' @inheritParams misfit_shmax
#' @param unc Uncertainty of observed \eqn{\sigma_\text{Hmax}}{SHmax}, either a
#' numeric vector or a number
#' @return Numeric vector
#' @details
#' The normalized \eqn{\chi^2}{chi-squared} test is
#' \deqn{ \text{Norm} \chi^2_i =
#'  = \frac{
#'    \sum^M_{i = 1} \left( \frac{\alpha_i - \alpha_{\text{predict}}}{\sigma_i} \right) ^2}
#'    {\sum^M_{i = 1} \left( \frac{90}{\sigma_i} \right) ^2 }}{
#'    (sum( ((obs-prd)/unc)^2 ) / sum( (90/unc)^2 )
#'    }
#' The value of the chi-squared test statistic is a number between 0 and 1
#' indicating the quality of the predicted \eqn{\sigma_\text{Hmax}}{SHmax}
#' directions. Low values
#' (\eqn{\le 0.15}) indicate good agreement,
#' high values (\eqn{> 0.7}) indicate a systematic misfit between predicted and
#' observed \eqn{\sigma_\text{Hmax}}{SHmax} directions.
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
#' @description Calculate the mean, median, quartile, interquartile range, variance,
#' deviation, and error of orientation data.
#'
#' @param x Numeric vector (NA values will be removed).
#' @param quiet logical. If false, a warning message is printed if there are NA
#' values.
#'
#' @return Numeric vector
#'
#' @details Quasi median on the circle, quasi quartiles on a circle, quasi
#' interquartile range on a circle.
#'
#' @importFrom stats median
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
#' @examples
#' x <- c(0, 45, 55, 40 + 180, 50 + 180)
#' circular_quasi_mean(x)
#' circular_quasi_median(x)
#' circular_quasi_quartile(x)
#' circular_quasi_IQR(x)
#' circular_var(x)
#' circular_mean_deviation(x)
#' circular_median_deviation(x)
#' circular_mean_error(x)
#' @name circle_median
NULL

#' @rdname circle_median
#' @export
circular_quasi_median <- function(x, quiet = TRUE) {
  stopifnot(any(is.numeric(x)))

  if (NA %in% x & quiet) {
    message("NA values have been dropped\n")
  }

  x <- deg2rad(sort(x[!is.na(x)])) %% pi
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
circular_quasi_mean <- function(x, quiet = TRUE) {
  stopifnot(any(is.numeric(x)))

  data <- data.frame(x)
  if (quiet) {
    data <- subset(data, !is.na(x))
  }

  x <- deg2rad(data$x)

  s.m <- sum(sin(2 * x))
  c.m <- sum(cos(2 * x))

  (atan2d_spec(s.m, c.m) / 2) %% 180
}

#' @rdname circle_median
#' @export
circular_quasi_quartile <- function(x, quiet = TRUE) {
  stopifnot(any(is.numeric(x)))

  if (NA %in% x & quiet) {
    message("NA values have been dropped\n")
  }
  x <- sort(x[!is.na(x)]) %% 180
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
    return(quantiles)
  } else {
    message("x needs more than 3 values\n")
    return(NULL)
  }
}

#' @rdname circle_median
#' @export
circular_quasi_IQR <- function(x, quiet = TRUE) {
  stopifnot(any(is.numeric(x)))

  if (NA %in% x & quiet) {
    message("NA values have been dropped\n")
  }

  quantiles <- circular_quasi_quartile(x)
  deviation_norm(as.numeric(quantiles[4] - quantiles[2]))
}

#' @rdname circle_median
#' @export
circular_var <- function(x, quiet = TRUE) {
  stopifnot(any(is.numeric(x)))

  if (NA %in% x & quiet) {
    message("NA values have been dropped\n")
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
#' @param x Numeric vector
#' @param unc Numeric vector or number of the uncertainties of x.
#' @param na.rm  logical value indicating whether \code{NA} values in x should
#' be stripped before the computation proceeds.
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
#' x <- c(175, 179, 2, 4)
#' unc <- c(10, 20, 10, 20)
#' circular_weighted_mean(x, unc)
#' circular_weighted_median(x, unc)
#' circular_weighted_sd(x, unc)
#' @name weighted_circle_stats
NULL

#' @rdname weighted_circle_stats
#' @export
circular_weighted_mean <- function(x, unc, na.rm = TRUE) {
  stopifnot(any(is.numeric(x)) & any(is.numeric(unc)))

  data <- data.frame(x, unc)
  if (na.rm) {
    data <- subset(data, !is.na(x) & !is.na(unc))
  }

  x <- deg2rad(data$x)
  unc <- deg2rad(data$unc)

  w <- 1 / unc
  Z <- sum(w)

  s.m <- 1 / Z * sum(w * sin(2 * x))
  c.m <- 1 / Z * sum(w * cos(2 * x))

  (atan2d_spec(s.m, c.m) / 2) %% 180
}

#' @rdname weighted_circle_stats
#' @export
circular_weighted_median <- function(x, unc, na.rm = TRUE) {
  stopifnot(any(is.numeric(x)) & any(is.numeric(unc)))

  data <- data.frame(x, unc)

  if (na.rm) {
    data <- subset(data, !is.na(x) & !is.na(unc))
  }

  data <- data[order(x), ]

  x <- deg2rad(data$x)
  unc <- deg2rad(data$unc)

  n <- length(x)
  w <- 1 / unc

  if (n %% 2 != 0) { # if odd
    m <- (n - 1) / 2
    # atand(
    #   sind(x[m+1]) / cosd(x[m+1])
    # ) %% 180

    s.m <- w[m + 1] * sin(2 * x[m + 1])
    c.m <- w[m + 1] * cos(2 * x[m + 1])
  } else { # if even
    m <- n / 2
    # atand(
    #   (sind(x[m]) + sind(x[m + 1])) /
    #     (cosd(x[m]) + cosd(x[m + 1]))
    # ) %% 180

    s.m <- (w[m] * sin(2 * x[m]) + w[m + 1] * sin(2 * x[m + 1]))
    c.m <- (w[m] * cos(2 * x[m]) + w[m + 1] * cos(2 * x[m + 1]))
  }
  atan2d_spec(s.m, c.m) %% 180
}


#' @rdname weighted_circle_stats
#' @export
circular_weighted_sd <- function(x, unc, na.rm = TRUE) {
  stopifnot(any(is.numeric(x)) & any(is.numeric(unc)))

  data <- data.frame(x, unc)
  if (na.rm) {
    data <- subset(data, !is.na(x) & !is.na(unc))
  }

  x <- deg2rad(data$x)
  unc <- deg2rad(data$unc)

  w <- 1 / unc
  Z <- sum(w)

  s.m <- 1 / Z * sum(w * sin(2 * x))
  c.m <- 1 / Z * sum(w * cos(2 * x))
  R <- sqrt(c.m^2 + s.m^2)

  sd <- sqrt((-2 * log(R)) / 2)
}
