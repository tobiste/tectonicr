#' Apply Rolling Functions using Circular Statistics
#'
#' A generic function for applying a function to rolling margins of an array.
#'
#' @inheritParams circular_mean
#' @param width integer specifying the window width (in numbers of observations)
#' which is aligned to the original sample according to the `align` argument.
#' If `NULL`, an optimal width is calculated.
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
                           width = NULL, by.column = FALSE,
                           partial = TRUE,
                           fill = NA,
                           ...) {
  FUN <- match.fun(FUN)

  if (is.null(width)) {
    width <- optimal_rollwidth(x)
  }

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

#' Apply Rolling Functions using Circular Statistical Tests for Uniformity
#'
#' A generic function for applying a function to rolling margins of an array.
#'
#' @inheritParams norm_chisq
#' @param x numeric. Directions in dgerees
#' @inheritParams circular_dispersion
#' @param conf.level Level of confidence: \eqn{(1 - \alpha \%)/100}.
#' (`0.95` by default).
#' @param R The number of bootstrap replicates.
#' @param width integer specifying the window width (in numbers of observations)
#' which is aligned to the original sample according to the `align` argument.
#' If `NULL`, an optimal width is estimated.
#' @param by.column logical. If `TRUE`, FUN is applied to each column separately.
#' @param fill a three-component vector or list (recycled otherwise) providing
#' filling values at the left/within/to the right of the data range. See the
#' fill argument of [zoo::na.fill()] for details
#' @param partial logical or numeric. If `FALSE` then `FUN` is only
#' applied when all indexes of the rolling window are within the observed time
#' range. If `TRUE` (default), then the subset of indexes that are in range
#' are passed to `FUN`. A numeric argument to partial can be used to determine
#' the minimal window size for partial computations. See below for more details.
#' @param ... optional arguments passed to [zoo::rollapply()]
#' @returns numeric vector with the test statistic of the rolling test.
#' `roll_dispersion_CI` returns a 2-column matrix with the lower and the upper confidence limits
#' @note If the rolling functions are applied to values that are a function of
#' distance it is recommended to sort the values first.
#' @importFrom zoo rollapply
#' @name rolling_test
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
#' roll_normchisq(dat.PoR$azi.PoR, 135, dat$unc)
#' roll_rayleigh(dat.PoR$azi.PoR, prd = 135, unc = dat$unc)
#' roll_dispersion(dat.PoR$azi.PoR, mean = 135, w = 1 / dat$unc)
#' roll_confidence(dat.PoR$azi.PoR, w = 1 / dat$unc)
#' roll_dispersion_CI(dat.PoR$azi.PoR, mean = 135, w = 1 / dat$unc, R = 10)
NULL

#' @rdname rolling_test
#' @export
roll_normchisq <- function(obs, prd, unc = NULL,
                           width = NULL, by.column = FALSE,
                           partial = TRUE,
                           fill = NA,
                           ...) {
  if (is.null(width)) {
    width <- optimal_rollwidth(obs)
  }

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

#' @rdname rolling_test
#' @export
roll_rayleigh <- function(obs, prd, unc = NULL,
                          width = NULL, by.column = FALSE,
                          partial = TRUE,
                          fill = NA,
                          ...) {
  if (is.null(width)) {
    width <- optimal_rollwidth(obs)
  }

  zoo::rollapply(
    cbind(obs, prd, unc),
    width = width,
    FUN = function(x) {
      suppressMessages(weighted_rayleigh(x[, 1], x[, 2], x[, 3]))$statistic
    },
    by.column = by.column,
    partial = partial,
    fill = fill,
    ...
  )
}

#' @rdname rolling_test
#' @export
roll_dispersion <- function(x, mean, w = NULL,
                            width = NULL, by.column = FALSE,
                            partial = TRUE,
                            fill = NA,
                            ...) {
  if (is.null(width)) {
    width <- optimal_rollwidth(x)
  }

  zoo::rollapply(
    cbind(x, mean, w),
    width = width,
    FUN = function(x) {
      suppressMessages(circular_dispersion(x[, 1], x[, 2], x[, 3], norm = TRUE))
    },
    by.column = by.column,
    partial = partial,
    fill = fill,
    ...
  )
}

#' @rdname rolling_test
#' @export
roll_confidence <- function(x, conf.level = .95, w = NULL, axial = TRUE,
                            width = NULL, by.column = FALSE, partial = TRUE,
                            fill = NA,
                            ...) {
  if (is.null(width)) {
    width <- optimal_rollwidth(x)
  }

  zoo::rollapply(
    cbind(x, rep(conf.level, length(x)), w),
    width = width,
    FUN = function(x) {
      confidence_angle(x[, 1], x[1, 2], x[, 3])
    },
    by.column = by.column,
    partial = partial,
    fill = fill,
    ...
  )
}

#' @rdname rolling_test
#' @export
roll_dispersion_CI <- function(x, mean, w = NULL, R, conf.level = .95,
                               width = NULL, by.column = FALSE, partial = TRUE, fill = NA, ...) {
  if (is.null(width)) {
    width <- optimal_rollwidth(x)
  }

  zoo::rollapply(
    cbind(x, mean, w),
    width = width,
    FUN = function(x) {
      suppressMessages(circular_dispersion_MLE(x[, 1], x[, 2], x[, 3], R = R, conf.level = conf.level, ...)$CI)
    },
    by.column = by.column,
    partial = partial,
    fill = fill,
    ...
  )
}


optimal_rollwidth <- function(x) {
  round((2 * circular_IQR(x) / length(x)^(1 / 3)) * 2 * pi)
}