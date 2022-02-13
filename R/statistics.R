#' Normalized Chi-Square Test
#'
#' A quantitative comparison between the predicted and observed directions of
#' \eqn{\sigma_\text{Hmax}}{SHmax} is obtained by the calculation of the average azimuth and by a
#' normalized \eqn{\chi^2}{chi-square} test.
#'
#' @references Wdowinski, S., 1998, A theory of intraplate
#'   tectonics: Journal of Geophysical Research: Solid Earth, v. 103, p.
#'   5037-5059, \doi{10.1029/97JB03390}.
#' @inheritParams misfit_shmax
#' @param unc Uncertainty of observed \eqn{\sigma_\text{Hmax}}{SHmax}, either numeric vector of length of
#' \code{obs} or a number
#' @return Numeric vector
#' @details
#' The normalized \eqn{\chi^2}{chi-square} test is
#' \deqn{ \text{Norm} \chi^2_i =
#'  = \frac{
#'    \sum^M_{i = 1} \left( \frac{\alpha_i - \alpha_{\text{predict}}}{\sigma_i} \right) ^2}
#'    {\sum^M_{i = 1} \left( \frac{90}{\sigma_i} \right) ^2 }}
#' The test result are values are between 0 and 1 indicating the quality of
#' the predicted \eqn{\sigma_\text{Hmax}}{SHmax} directions. Low values (\eqn{\le 0.15}) indicate good agreement,
#' high values (\eqn{> 0.7}) indicate a systematic misfit between predicted and
#' observed \eqn{\sigma_\text{Hmax}}{SHmax} directions.
#' @export
#' @examples
#' data("nuvel1")
#' euler <- subset(nuvel1, nuvel1$ID == "na") # North America relative to
#' # Pacific plate
#' point <- data.frame(lat = 45, lon = 20)
#' prd <- model_shmax(point, euler)
#' norm_chi2(obs = 90, prd$sc, unc = 10)
norm_chi2 <- function(obs, prd, unc) {
  if (length(prd) != length(obs)) {
    stop("Observed and predicted values must have the same length!\n")
  }

  if (length(unc) != 1 & length(unc) != length(obs)) {
    stop("Uncertainties must be either a numeric value or a numeric value with
         length of observed values!\n")
  }

  if (length(unc) == 1) {
    unc <- rep(unc, length(obs))
  }

  w <- c()
  x <- c()
  y <- c()
  for (i in seq_along(obs)) {
    if (is.na(prd[i]) | is.na(obs[i])) {
      x[i] <- NA
      y[i] <- NA
    } else {
      w[i] <- deviation_norm(prd[i] - (obs[i] + 180) %% 180)
      x[i] <- (w[i] / unc[i])^2
      y[i] <- (90 / unc[i])^2
    }
  }
  nchi2 <- sum(x, na.rm = TRUE) / sum(y, na.rm = TRUE)

  return(nchi2)
}


#' @title Median and Quartile on Pi-periodic Data
#'
#' @description Calculate the median, quartile, and interquartile range of
#' orientation data.
#'
#' @param x Numeric vector
#'
#' @return Numeric vector
#'
#' @details Quasi median on the circle, quasi quartiles on a circle, quasi interquartile range on a circle.
#'
#' @source [stats::median()], [stats::stats()], and [stats::IQR()] are the
#' equivalents for non-periodic data.
#'
#' @references
#' * Ratanaruamkarn, S., Niewiadomska-Bugaj, M., Wang, J.-C. (2009).
#' A New Estimator of a Circular Median. Communications in Statistics -
#' Simulation and Computation, 38(6), 1269-1291.
#' \doi{10.1080/03610910902899950}.
#'
#' * Reiter, K., Heidbach, O., Schmitt, D., Haug, K., Ziegler, M.,
#' Moeck, I. (2014). A revised crustal stress orientation database for Canada.
#' Tectonophysics, 636, 111-124. \doi{10.1016/j.tecto.2014.08.006}.
#' @examples
#' x <- c(0, 45, 55, 40 + 180, 50 + 180)
#' circular_quasi_median(x)
#' circular_quasi_quartile(x)
#' circular_quasi_interquartile_range(x)
#' @name circle_median
NULL

#' @rdname circle_median
#' @export
circular_quasi_median <- function(x) {
  stopifnot(any(is.numeric(x)))

  if (NA %in% x) {
    message("NA values have been dropped\n")
  }
  x <- sort(x[!is.na(x)])
  n <- length(x)

  if (n %% 2 != 0) {
    m <- (n - 1) / 2
    qmed <- atand(
      sind(x[m + 1]) / cosd(x[m + 1])
    )
  } else {
    m <- n / 2
    qmed <- atand(
      (sind(x[m]) + sind(x[m + 1])) /
        (cosd(x[m]) + cosd(x[m + 1]))
    )
  }
  return(qmed)
}

#' @rdname circle_median
#' @export
circular_quasi_quartile <- function(x) {
  stopifnot(any(is.numeric(x)))

  if (NA %in% x) {
    message("NA values have been dropped\n")
  }
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
    return(quantiles)
  } else {
    message("x needs more than 3 values\n")
    return(NULL)
  }
}

#' @rdname circle_median
#' @export
circular_quasi_interquartile_range <- function(x) {
  stopifnot(any(is.numeric(x)))

  if (NA %in% x) {
    message("NA values have been dropped\n")
  }
  x <- sort(x[!is.na(x)])

  quantiles <- circular_quasi_quartile(x)
  qiroc <- as.numeric(quantiles[4] - quantiles[2])
  return(qiroc)
}
