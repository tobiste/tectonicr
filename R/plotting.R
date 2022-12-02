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
    binwidth <- 2 * circular_IQR(x) / length(stats::na.omit(x))^(1 / 3)
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

#' Plotting the \eqn{\sigma_{Hmax}}{SHmax} azimuth
#'
#' Creates a set of plots including
#' the azimuth as a function of the distance to the plate boundary,
#' the Norm Chi-squared as a function of the distance to the plate boundary,
#' and a rose diagram of the frequency distribution of the azimuths.
#'
#' @param azi numeric. Azimuth of \eqn{\sigma_{Hmax}}{SHmax}
#' @param distance numeric. Distance to plate boundary
#' @param prd numeric. the predicted direction of \eqn{\sigma_{Hmax}}{SHmax}
#' @param unc numeric. Uncertainty of observed \eqn{\sigma_{Hmax}}{SHmax}, either a
#' numeric vector or a number
#' @param regime character vector. The stress
#' regime (following the classification of the World Stress Map)
#' @param k integer. window width (in number of observations) for rolling
#' statistics. Has to be an odd number.
#' @param ... optional arguments to `zoo::rollapply()`, [rose()], or `plot()`
#' @importFrom zoo rollapply rollmedian
#' @seealso [PoR_shmax()], [distance_from_pb()], [circular_median()], [circular_IQR()], [norm_chisq()]
#' @export
#' @examples
#' data("nuvel1")
#' na_pa <- subset(nuvel1, nuvel1$plate.rot == "na")
#'
#' data("plates")
#' plate_boundary <- subset(plates, plates$pair == "na-pa")
#'
#' data("san_andreas")
#' res <- PoR_shmax(san_andreas, na_pa, "right")
#' d <- distance_from_pb(san_andreas, na_pa, plate_boundary, tangential=TRUE)
#' PoR_plot(res$azi.PoR, d, res$prd, san_andreas$unc, san_andreas$regime)
PoR_plot <- function(azi, distance, prd, unc, regime, k = 51, ...){
  stopifnot(k >= 3, k %% 2 == 1)
  if(missing(regime)){
    regime <- rep(NA, length(azi))
  }
  nchisq_i <- NULL
  regime <- ifelse(is.na(regime), "U", regime)
  t <- data.frame(azi, distance, prd, unc, regime = factor(regime, levels = c("U", "N", "NS", "S", "TS", "T"))) %>%
    arrange(distance) %>%
    mutate(nchisq_i =  (deviation_norm(azi - prd)/unc)^2 / (90/unc)^2,
           azi.rmedian = zoo::rollapply(azi, width = k, FUN = circular_median, align = "center", fill = NA, ...),
           azi.riqr = zoo::rollapply(azi, width = k, FUN = circular_IQR, align = "center", fill = NA, ...),
           nchi2.rmedian = zoo::rollmedian(nchisq_i, k = k, align = "center", fill = NA, ...),
           nchi2.rmad = zoo::rollapply(nchisq_i, width = k, FUN = stats::mad, align = "center", fill = NA, na.rm = TRUE, ...)
    )

  nchisq <- norm_chisq(azi, prd, unc)
  azi.PoR.median <- circular_median(azi, 1 / unc)
  azi.PoR.IQR <- circular_IQR(azi, 1 / unc)

  subtitle <- paste0(
    "N: ", length(azi),
    " | Median azimuth: ", round(azi.PoR.median, 1), "\u00B0 \u00B1 ", round(azi.PoR.IQR / 2, 1),
    "\u00B0 | Norm \u03C7\u00B2:", round(nchisq, 2)
  )

  grDevices::palette(c("grey60","#D55E00", "#E69F00", "#009E73", "#56B4E9", "#0072B2"))

  plot(0, type="n",
       xlab = "Distance from plate boundary", ylab = "Azimuth wrt. EP (\u00B0)",
       sub = subtitle,
       xlim = range(distance),
       ylim = c(0, 180), yaxp = c(0, 180, 8), ...
  )
  graphics::arrows(y0=t$azi-t$unc, x0=t$distance, y1=t$azi+t$unc, x1=t$distance, code=0, lwd=.25, col = t$regime)
  graphics::points(azi ~ distance, data = t, col = t$regime)

  graphics::lines(azi.rmedian-azi.PoR.IQR/2 ~ distance, data = t, type = "S", col = "#85112A7D", lty =3)
  graphics::lines(azi.rmedian+azi.PoR.IQR/2 ~ distance, data = t, type = "S", col = "#85112A7D", lty =3)
  graphics::lines(azi.rmedian ~ distance, data = t, type = "S", col = "#85112AFF")
  graphics::abline(h = unique(prd), col = 'black', lty = 2)
  graphics::legend("bottomright", inset = .05, cex = .5, legend=c("N", "NS", "S", "TS", "T", "U"), title = "Stress regime", fill=c("#D55E00", "#E69F00", "#009E73", "#56B4E9", "#0072B2", "grey60" ))

  grDevices::dev.new()
  plot(nchisq_i ~ distance, data = t, col = t$regime,
       xlab = "Distance from plate boundary", ylab = expression(Norm ~ chi[i]^2),
       #sub = subtitle,
       xlim = range(distance),
       ylim = c(0, 1), yaxp = c(0, 1, 4), ...
  )
  graphics::lines(nchi2.rmedian+t$nchi2.rmad ~ distance, data = t, type = "S", col = "#85112A7D", lty=3)
  graphics::lines(nchi2.rmedian-t$nchi2.rmad ~ distance, data = t, type = "S", col = "#85112A7D", lty=3)
  graphics::lines(nchi2.rmedian ~ distance, data = t, type = "S", col = "#85112AFF")
  graphics::abline(h = .15, col = 'black', lty = 2)

  grDevices::dev.new()
  rose(azi, main = "Euler pole", ...)
}
