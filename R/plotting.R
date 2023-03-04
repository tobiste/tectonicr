#' @title Selecting optimal number of bins and width for rose diagrams
#'
#' @param n Integer. number of data
#' @param round Logical. Whether bin width is round to zero digits (`round=TRUE`, the default)
#' or as is (`FALSE`).
#' @param axial Logical. Whether data are uniaxial (`axial=FALSE`)
#' or biaxial (`TRUE`, the default).
#' @param ... Additional arguments passed to [rose_bw()].
#' @name rose_bw
NULL

#' @rdname rose_bw
rose_bins <- function(n, round = TRUE) {
  b <- 2 * n^(1 / 3) # bins
  if (round) {
    round(b)
  } else {
    b
  }
}

#' @rdname rose_bw
rose_binwidth <- function(n, axial = TRUE, ...) {
  if (axial) {
    r <- 180
  } else {
    r <- 360
  }
  r / rose_bins(n, ...)
}



#' @title Rose Diagram
#'
#' @description Plots a rose diagram (rose of directions), the analogue of a
#' histogram or density plot for angular data.
#'
#' @param x Data to be plotted. A numeric vector containing angles.
#' @param weights Optional vector of numeric weights associated with x.
#' @param binwidth The width of the bins.
#' @param bins number of arcs to partition the circle width.
#' Overridden by `binwidth`.
#' @param axial Logical. Whether data are uniaxial (`axial=FALSE`)
#' or biaxial (`TRUE`, the default).
#' @param equal_area Logical. Whether the radii of the bins are proportional to
#' the frequencies (`equal_area=FALSE`, i.e. equal-angle) or proportional to the
#' square-root of the frequencies (`equal_area=TRUE`, the default).
#' @param clockwise Logical. Whether angles increase in the
#' clockwise direction (`clockwise=TRUE`, the default) or anti-clockwise,
#' counter-clockwise direction (`FALSE`).
#' @param unit The unit in which the angles are expressed.
#' `"degrees"` (the default), or `"radians"`.
#' @param round_binwidth Logical. Whether bin width is round to zero digits
#' (`round_binwidth=TRUE`, the default) or as is (`FALSE`).
#' @param main,sub Character string specifying the title and subtitle of the
#' plot. If `sub = NULL`, it will show the bin width.
#' @param at Optional vector of angles at which tick marks should be plotted.
#' Set `at=numeric(0)` to suppress tick marks.
#' @param add_pts logical. Whether a circular dot plot should be added
#' (`FALSE` is the default).
#' @param pts_cex,pts_pch,pts_col Plotting arguments for circular dot plot
#' @param ... Additional arguments passed to [spatstat.explore::rose()].
#' @note If `bins` and `binwidth` are `NULL`, an optimal bin width will be
#' calculated using Scott (1979):
#' \deqn{ \frac{R}{n^{\frac{1}{3}}}
#' }
#' with n being the length of `x`, and the range R being either 180 or 360
#' degree for axial or directional data, respectively.
#' @return A window (class `"owin"`) containing the plotted region.
#' @importFrom spatstat.explore rose
#' @importFrom graphics hist title points
#' @importFrom stats na.omit
#' @export
#' @examples
#' x <- runif(100, 60, 210)
#' rose(x)
#'
#' data("san_andreas")
#' rose(san_andreas$azi, col = "grey", axial = TRUE, stack = TRUE)
#' rose(san_andreas$azi, weights = 1 / san_andreas$unc, col = "grey", axial = TRUE)
rose <- function(x, weights = NULL, binwidth = NULL, bins = NULL, axial = TRUE,
                 equal_area = TRUE, clockwise = TRUE, unit = c("degree", "radian"),
                 round_binwidth = TRUE, main = "N", sub, at = seq(0, 360 - 45, 45),
                 add_pts = FALSE, pts_pch = 1, pts_cex = 1, pts_col = "grey", ...) {
  x <- as.vector(x %% 360)

  if (!is.null(bins) && is.null(binwidth)) {
    bins <- round(bins)
    stopifnot(bins > 0)
    binwidth <- 360 / bins # bin width
  } else if (is.null(bins) && is.null(binwidth)) {
    # binwidth <- 2 * circular_IQR(x) / length(stats::na.omit(x))^(1 / 3)
    binwidth <- rose_binwidth(length(stats::na.omit(x)), axial = axial)
  }

  if (round_binwidth) binwidth <- round(binwidth)

  stopifnot(binwidth > 0)
  breaks <- seq(0, 360, binwidth)
  if (!(360 %in% breaks)) {
    breaks <- c(breaks, 360)
  }

  if (axial) {
    x2 <- (x + 180) %% 360 # add data to the other side of the circle
    # x <- graphics::hist(x = c(x, x2), plot = FALSE, breaks = breaks)
    x <- c(x, x2)
    weights <- c(weights, weights)
  }

  freqs <- graphics::hist(x = x, plot = FALSE, breaks = breaks)

  if (equal_area) {
    freqs$density <- sqrt(freqs$density)
  }

  spatstat.explore::rose(
    freqs,
    weights = weights,
    breaks = breaks,
    clockwise = clockwise, start = "N", unit = unit, main = main, xlab = NULL,
    at = seq(0, 360 - 45, 45),
    ...
  )

  if (add_pts) {
    scale <- 1.1 * max(freqs$density)
    u <- deg2rad(90 - x)
    n <- length(x)
    z <- cos(u) * scale
    y <- sin(u) * scale
    # if (stack == FALSE) {
    graphics::points(z, y, cex = pts_cex, pch = pts_pch, col = pts_col)
    # } else {
    #   bins <- 180/binwidth
    #   bins.count <- c(1:bins)
    #   arc <- (2 * pi) / bins
    #   for (i in 1:bins) {
    #     bins.count[i] <- sum(u <= i * arc & u > (i - 1) * arc)
    #   }
    #   mids <- seq(arc / 2, 2 * pi - pi / bins, length = bins)
    #   index <- pts_cex / dotsep
    #   for (i in 1:bins) {
    #     if (bins.count[i] != 0) {
    #       for (j in 0:(bins.count[i] - 1)) {
    #         r <- 1 + j * index
    #         z <- r * cos(mids[i]) * scale
    #         y <- r * sin(mids[i]) * scale
    #         points(z, y, cex = pts_cex, pch = pts_pch)
    #       }
    #     }
    #   }
    # }
  }

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
#' @importFrom dplyr arrange mutate
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
#' d <- distance_from_pb(san_andreas, na_pa, plate_boundary, tangential = TRUE)
#' PoR_plot(res$azi.PoR, d, res$prd, san_andreas$unc, san_andreas$regime)
PoR_plot <- function(azi, distance, prd, unc = NULL, regime, k = 51, ...) {
  stopifnot(k >= 3, k %% 2 == 1)
  if (missing(regime)) {
    regime <- rep(NA, length(azi))
  }
  nchisq_i <- NULL
  regime <- ifelse(is.na(regime), "U", regime)



  t <- data.frame(azi, distance, prd, unc, regime = factor(regime, levels = c("U", "N", "NS", "S", "TS", "T"))) %>%
    dplyr::arrange(distance) %>%
    dplyr::mutate(
      nchisq_i = (deviation_norm(azi - prd) / unc)^2 / (90 / unc)^2,
    )

  t$roll_mean = zoo::rollapply(
    t %>% select(azi, unc),
    width = k,
    FUN = function(x){circular_mean(x[,"azi"],1/x[,"unc"])},
    by.column = FALSE,
    partial = TRUE,
    align = "center",
    fill = NA,
    ...
    )
  t$roll_sd = zoo::rollapply(
    t %>% select(azi, unc),
    width = k,
    FUN = function(x){circular_sd(x[,"azi"],1/x[,"unc"])},
    by.column = FALSE,
    partial = TRUE,
    align = "center",
    fill = NA,
    ...
    )
  t$roll_nchisq <- zoo::rollapply(
    t %>% select(azi, prd, unc),
    width = k,
    FUN = function(x){norm_chisq(x[,"azi"],x[,"prd"],x[,"unc"])},
    by.column = FALSE,
    partial = TRUE,
    align = "center",
    fill = NA,
    ...
  )

  nchisq <- norm_chisq(azi, prd, unc)
  azi.PoR.mean <- circular_mean(azi, 1 / unc)
  azi.PoR.sd <- circular_sd(azi, 1 / unc)

  subtitle <- paste0(
    "N: ", length(azi),
    " | Mean azimuth: ", round(azi.PoR.mean, 1), "\u00B0 \u00B1 ", round(azi.PoR.sd, 1),
    "\u00B0 | Norm \u03C7\u00B2: ", round(nchisq, 2)
  )

  grDevices::palette(c("grey60", "#D55E00", "#E69F00", "#009E73", "#56B4E9", "#0072B2"))

  plot(0,
    type = "n",
    xlab = "Distance from plate boundary", ylab = "Azimuth wrt. PoR (\u00B0)",
    sub = subtitle,
    xlim = range(distance),
    ylim = c(0, 180), yaxp = c(0, 180, 8), ...
  )
  graphics::arrows(y0 = t$azi - t$unc, x0 = t$distance, y1 = t$azi + t$unc, x1 = t$distance, code = 0, lwd = .25, col = t$regime)
  graphics::points(azi ~ distance, data = t, col = t$regime)

  graphics::lines(roll_mean - roll_sd ~ distance, data = t, type = "S", col = "#85112A7D", lty = 3)
  graphics::lines(roll_mean + roll_sd ~ distance, data = t, type = "S", col = "#85112A7D", lty = 3)
  graphics::lines(roll_mean ~ distance, data = t, type = "S", col = "#85112AFF")
  graphics::abline(h = unique(prd), col = "black", lty = 2)
  graphics::legend("bottomright", inset = .05, cex = .5, legend = c("N", "NS", "S", "TS", "T", "U"), title = "Stress regime", fill = c("#D55E00", "#E69F00", "#009E73", "#56B4E9", "#0072B2", "grey60"))

  grDevices::dev.new()
  plot(nchisq_i ~ distance,
    data = t, col = t$regime,
    xlab = "Distance from plate boundary", ylab = expression(Norm ~ chi^2),
    # sub = subtitle,
    xlim = range(distance),
    ylim = c(0, 1), yaxp = c(0, 1, 4), ...
  )
  graphics::lines(roll_nchisq ~ distance, data = t, type = "S", col = "#85112AFF")
  graphics::abline(h = .15, col = "black", lty = 2)

  grDevices::dev.new()
  rose(azi, weights = 1 / unc, main = "PoR", ...)
}
