#' @title Rose Diagram
#'
#' @description Plots a rose diagram (rose of directions), the analogue of a
#' histogram or density plot for angular data.
#'
#' @param x Data to be plotted. A numeric vector containing angles.
#' @param binwidth The width of the bins.
#' @param bins number of arcs to partition the circle width.
#' Overridden by `binwidth`.
#' @param axial Logical. Whether data are uniaxial (`axial=FALSE`)
#' or biaxial (`TRUE`, the default).
#' @param clockwise Logical. Whether angles increase in the
#' clockwise direction (`clockwise=TRUE`, the default) or anti-clockwise,
#' counter-clockwise direction (`FALSE`).
#' @param unit The unit in which the angles are expressed.
#' `"degrees"` (the default), or `"radians"`.
#' @param main,sub Character string specifying the title and subtitle of the
#' plot. If `sub = NULL`, it will show the bin width.
#' @param ... Additional arguments passed to [spatstat.core::rose()].
#' @note If `bins` and `binwidth` are `NULL`, an optimal bin width will be
#' calculated using:
#' \deqn{ \frac{2 IQR(x)}{n^{\frac{1}{3}}}
#' }
#' with n being the length of `x`.
#' @return A window (class `"owin"`) containing the plotted region.
#' @importFrom spatstat.core rose
#' @importFrom graphics hist title
#' @importFrom stats na.omit
#' @export
#' @examples
#' x <- runif(100, 60, 210)
#' rose(x)
#'
#' data("san_andreas")
#' rose(san_andreas$azi, col = "grey", axial = TRUE)
rose <- function(x, binwidth = NULL, bins = NULL, axial = TRUE, clockwise = TRUE, unit = c("degree", "radian"), main = "N", sub, ...) {
  x <- as.vector(x %% 360)


  if (!is.null(bins) & is.null(binwidth)) {
    bins <- round(bins)
    stopifnot(bins > 0)
    binwidth <- 360 / bins # bin width
  } else if (is.null(bins) & is.null(binwidth)) {
    # bins <- length(x)
    binwidth <- 2 * circular_quasi_IQR(x) / length(stats::na.omit(x))^(1 / 3)
  }
  stopifnot(binwidth > 0)
  breaks <- seq(0, 360, binwidth)
  if (!(360 %in% breaks)) {
    breaks <- c(breaks, 360)
  }

  if (axial) {
    x2 <- (x + 180) %% 360 # add data to the other side of the circle
    x <- graphics::hist(x = c(x, x2), plot = FALSE, breaks = breaks)
  }
  spatstat.core::rose(
    x,
    breaks = breaks,
    clockwise = clockwise, start = "N", unit = unit, main = main, xlab = NULL, ...
  )

  if (missing(sub)) sub <- paste0("Bin width: ", binwidth)
  graphics::title(sub = sub, ylab = NULL)
}
