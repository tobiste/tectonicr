#' Apply Rolling Functions using Circular Statistics
#'
#' A generic function for applying a function to rolling margins of an array.
#'
#' @inheritParams circular_mean
#' @param width integer specifying the window width (in numbers of observations)
#' which is aligned to the original sample according to the `align` argument.
#' If `NULL`, an optimal width is calculated.
#' @param FUN the function to be applied
#' @param by.column logical. If `TRUE`, `FUN` is applied to each column
#' separately.
#' @param fill a three-component vector or list (recycled otherwise) providing
#' filling values at the left/within/to the right of the data range. See the
#' fill argument of [zoo::na.fill()] for details
#' @param partial logical or numeric. If `FALSE` then `FUN` is only
#' applied when all indexes of the rolling window are within the observed time
#' range. If `TRUE` (default), then the subset of indexes that are in range
#' are passed to `FUN`. A numeric argument to partial can be used to determine
#' the minimal window size for partial computations. See below for more details.
#' @param ... optional arguments passed to [zoo::rollapply()]
#'
#' @returns numeric vector  with the results of the rolling function.
#'
#' @note If the rolling statistics are applied to values that are a function of
#' distance it is recommended to sort the values first.
#'
#' @importFrom zoo rollapply
#'
#' @export
#'
#' @examples
#' data("plates")
#' plate_boundary <- subset(plates, plates$pair == "na-pa")
#' data("san_andreas")
#' PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
#' distance <- distance_from_pb(
#'   x = san_andreas,
#'   PoR = PoR,
#'   pb = plate_boundary,
#'   tangential = TRUE
#' )
#' dat <- san_andreas[order(distance), ]
#' roll_circstats(dat$azi, w = 1 / dat$unc, circular_mean, width = 51) |> head()
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
#' @param x,y numeric. Directions in degrees
#' @param w,w.y (optional) Weights of `x` and `y`, respectively. A vector of
#' positive numbers and of the same length as `x`.
#' @inheritParams circular_dispersion
#' @param conf.level Level of confidence: \eqn{(1 - \alpha \%)/100}.
#' (`0.95` by default).
#' @param R The number of bootstrap replicates.
#' @param width integer specifying the window width (in numbers of observations)
#' which is aligned to the original sample according to the `align` argument.
#' If `NULL`, an optimal width is estimated.
#' @param by.column logical. If `TRUE`, FUN is applied to each column
#' separately.
#' @param fill a three-component vector or list (recycled otherwise) providing
#' filling values at the left/within/to the right of the data range. See the
#' fill argument of [zoo::na.fill()] for details
#' @param partial logical or numeric. If `FALSE` then `FUN` is only
#' applied when all indexes of the rolling window are within the observed time
#' range. If `TRUE` (default), then the subset of indexes that are in range
#' are passed to `FUN`. A numeric argument to partial can be used to determine
#' the minimal window size for partial computations. See below for more details.
#' @param ... optional arguments passed to [zoo::rollapply()]
#'
#' @returns numeric vector with the test statistic of the rolling test.
#' `roll_dispersion_CI` returns a 2-column matrix with the lower and the upper
#' confidence limits
#'
#' @note If the rolling functions are applied to values that are a function of
#' distance it is recommended to sort the values first.
#'
#' @importFrom zoo rollapply
#'
#' @name rolling_test
#'
#' @examples
#' data("plates")
#' plate_boundary <- subset(plates, plates$pair == "na-pa")
#' data("san_andreas")
#' PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
#' distance <- distance_from_pb(
#'   x = san_andreas,
#'   PoR = PoR,
#'   pb = plate_boundary,
#'   tangential = TRUE
#' )
#' dat <- san_andreas[order(distance), ]
#' dat.PoR <- PoR_shmax(san_andreas, PoR, "right")
#' roll_normchisq(dat.PoR$azi.PoR, 135, dat$unc) |> head()
#' roll_rayleigh(dat.PoR$azi.PoR, prd = 135, unc = dat$unc) |> head()
#' roll_dispersion(dat.PoR$azi.PoR, y = 135, w = 1 / dat$unc) |> head()
#' roll_confidence(dat.PoR$azi.PoR, w = 1 / dat$unc) |> head()
#' \donttest{
#' roll_dispersion_CI(dat.PoR$azi.PoR, y = 135, w = 1 / dat$unc, R = 10) |> head()
#' }
NULL

#' @rdname rolling_test
#' @export
roll_normchisq <- function(obs, prd, unc = NULL,
                           width = NULL, by.column = FALSE,
                           partial = TRUE,
                           fill = NA,
                           ...) {
  if (is.null(width)) width <- optimal_rollwidth(obs)

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
  if (is.null(width)) width <- optimal_rollwidth(obs)

  zoo::rollapply(
    cbind(obs, prd, 1 / unc),
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
roll_dispersion <- function(x, y, w = NULL, w.y = NULL,
                            width = NULL, by.column = FALSE,
                            partial = TRUE,
                            fill = NA,
                            ...) {
  if (is.null(width)) width <- optimal_rollwidth(x)
  if (is.null(w)) w <- rep(1, length(x))
  if (is.null(w.y)) w.y <- rep(1, length(x))

  zoo::rollapply(
    cbind(x, y, w, w.y),
    width = width,
    FUN = function(x) {
      suppressMessages(
        circular_dispersion(x[, 1], x[, 2], x[, 3], x[, 4])
      )
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
  if (is.null(width)) width <- optimal_rollwidth(x)

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
roll_dispersion_CI <- function(x, y, w = NULL, w.y = NULL, R, conf.level = .95,
                               width = NULL, by.column = FALSE, partial = TRUE, fill = NA, ...) {
  if (is.null(width)) width <- optimal_rollwidth(x)
  if (is.null(w)) w <- rep(1, length(x))
  if (is.null(w.y)) w.y <- rep(1, length(x))

  zoo::rollapply(
    cbind(x, y, w, w.y),
    width = width,
    FUN = function(x) {
      suppressMessages(circular_dispersion_boot(x[, 1], x[, 2], x[, 3], x[, 4], R = R, conf.level = conf.level, ...)$CI)
    },
    by.column = by.column,
    partial = partial,
    fill = fill,
    ...
  )
}

#' @rdname rolling_test
#' @export
roll_dispersion_sde <- function(x, y, w = NULL, w.y = NULL, R, conf.level = .95,
                                width = NULL, by.column = FALSE, partial = TRUE, fill = NA, ...) {
  if (is.null(width)) width <- optimal_rollwidth(x)
  if (is.null(w)) w <- rep(1, length(x))
  if (is.null(w.y)) w.y <- rep(1, length(x))

  zoo::rollapply(
    cbind(x, y, w, w.y),
    width = width,
    FUN = function(x) {
      suppressMessages(
        circular_dispersion_boot(
          x[, 1], x[, 2], x[, 3], x[, 4],
          R = R, conf.level = conf.level, ...
        )$sde
      )
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

#' Distance Binned Summary Statistics
#'
#' Circular summary statistics over intervals of distances.
#'
#' @param azi numeric. Azimuth values in degrees.
#' @param prd (optional) numeric. A predicted orientation in degrees.
#' @param prd.error (optional) numeric. The uncertainty of the predicted
#' orientation in degrees.
#' @param distance  numeric. the independent variable along the values in `azi`
#' are sorted, e.g. the plate boundary distances
#' @param n.breaks numeric. number (greater than or equal to 2) giving the
#' number of equal-sized intervals into which `distance` is to be cut.
#' Default is 10. Will be ignored if `width.breaks` is specified.
#' @param width.breaks numeric. The width of the intervals into which `distance`
#' is to be cut.
#' @param unc (optional) Uncertainties of `azi` (in degrees) acting as inverse
#' weighting factors for statistics.
#' @param kappa numeric. Concentration parameter applied for the circular mode.
#' @param R integer. Number of bootstrap iterates for estimating the error of
#' the dispersion.
#' @param conf.level The level of confidence for confidence interval and
#' bootstrapped standard error of dispersion.
#' @param ... optional arguments passed to [ggplot2::cut_number()] and [ggplot2::cut_width()]
#'
#' @return tibble containing the `n` values for `azi`in each bin, min/median/max
#' distance of the bin, and the summary statistics for `azi`.
#' If `prd` is specified, the normal Chi-squared statistic, dispersion and its
#' standard error are returned as well.
#'
#' @importFrom dplyr tibble mutate summarise select
#' @importFrom ggplot2 cut_width cut_number
#'
#' @seealso [circular_summary()], [circular_dispersion()], and [circular_dispersion_boot()]
#'
#' @export
#'
#' @examples
#' data("plates")
#' plate_boundary <- subset(plates, plates$pair == "na-pa")
#' data("san_andreas")
#' PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
#' san_andreas$distance <- distance_from_pb(
#'   x = san_andreas,
#'   PoR = PoR,
#'   pb = plate_boundary,
#'   tangential = TRUE
#' )
#' dat <- san_andreas |> cbind(PoR_shmax(san_andreas, PoR, "right"))
#'
#' distance_binned_stats(dat$azi.PoR,
#'   distance = dat$distance, width.breaks = 1,
#'   unc = dat$unc, prd = 135
#' ) |> head()
distance_binned_stats <- function(azi, distance, n.breaks = 10, width.breaks = NULL,
                                  unc = NULL, prd = NULL, prd.error = NULL,
                                  kappa = 2, R = 1000, conf.level = 0.95, ...) {
  # lifecycle::signal_stage("experimental", "distance_binned_stats()")
  stopifnot(length(azi) == length(distance))

  # binding of to be defined variables
  n <- w <- w.prd <- nchisq <- dispersion <- dispersion_sde <- NULL

  if (!is.null(width.breaks)) {
    bins <- ggplot2::cut_width(distance, width.breaks, ...)
  } else {
    bins <- ggplot2::cut_number(distance, n.breaks, ...)
  }
  no_prd <- is.null(prd)

  res <- dplyr::tibble(
    azi = azi,
    unc = unc,
    prd = prd,
    prderr = prd.error,
    distance = distance,
    bins = bins
  ) |>
    dplyr::mutate(
      # bins = cut(distance, n = n.breaks, ...),
      # bins = cut_fun(distance, cut_param, ...),
      w = ifelse(is.null(unc), 1, 1 / unc),
      w.prd = ifelse(is.null(prd.error), 1, 1 / prd.error),
      prd = ifelse(is.null(prd), NA, prd)
    ) |>
    dplyr::group_by(bins) |>
    dplyr::summarise(
      n = length(azi),
      distance_min = min(distance),
      distance_median = stats::median(distance, na.rm = TRUE),
      distance_max = max(distance),
      mean = ifelse(n > 2, circular_mean(azi, w = w), NA),
      sd = ifelse(n > 2, circular_sd(azi, w = w), NA),
      var = ifelse(n > 2, circular_var(azi, w = w), NA),
      lq = ifelse(n > 3, circular_quantiles(azi, w = w)[1], NA),
      quasimedian = ifelse(n > 2, circular_median(azi, w = w), NA),
      uq = ifelse(n > 3, circular_quantiles(azi, w = w)[3], NA),
      median = circular_sample_median(azi),
      mode = ifelse(n > 2, circular_mode(azi, kappa = kappa, kernel = "vonmises"), NA),
      CI = ifelse(n > 2, confidence_angle(azi, w = w, conf.level = conf.level), NA),
      skewness = ifelse(n > 2, second_central_moment(azi, w = w)$skewness, NA),
      kurtosis = ifelse(n > 2, second_central_moment(azi, w = w)$kurtosis, NA),
      nchisq = norm_chisq(azi, prd, w),
      dispersion = ifelse(n > 2, circular_dispersion(azi, prd, w = w, w.y = w.prd), NA),
      dispersion_sde = ifelse(n > 2, circular_dispersion_boot(azi, prd,
        w = w, w.y = w.prd,
        conf.level = conf.level, R = R
      )$sde, NA)
    )

  if (no_prd) {
    dplyr::select(res, -c(nchisq, dispersion, dispersion_sde))
  } else {
    res
  }
}
