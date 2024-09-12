#' Colors for input variables
#'
#' assigns colors to continuous or categorical values for plotting
#'
#' @param x values for color assignment
#' @param n integer. number of colors for continuous colors (i.e.
#' `categorical = FALSE``).
#' @param pal either a named vector specifying the colors for categorical
#' values, or a color function. If `NULL`, default colors are
#' `RColorBrewer::brewer.pal()`
#' (`categorical = TRUE`) and `viridis::viridis()` (`categorical = FALSE`).
#' @param categorical logical.
#' @param na.value color for `NA` values (categorical).
#' @param ... optional arguments passed to palette function
#'
#' @return named color vector
#'
#' @importFrom viridis viridis
#' @importFrom RColorBrewer brewer.pal
#' @importFrom dplyr left_join mutate
#' @export
#'
#' @examples
#' val1 <- c("N", "S", "T", "T", NA)
#' tectonicr.colors(val1, categorical = TRUE)
#' tectonicr.colors(val1, pal = stress_colors(), categorical = TRUE)
#'
#' val2 <- runif(10)
#' tectonicr.colors(val2, n = 5)
tectonicr.colors <- function(x, n = 10, pal = NULL, categorical = FALSE, na.value = "grey", ...) {
  code <- val <- NULL
  if (categorical) {
    # dat <- data.frame(val = x) # |> mutate(val = ifelse(is.na(val), "NA", val))
    if (is.null(pal)) {
      val <- unique(x)
      n <- length(val)
      if (n < 9) {
        pal <- structure(RColorBrewer::brewer.pal(n, name = "Set2"), names = val)
      } else {
        pal <- structure(RColorBrewer::brewer.pal(n, name = "Set3"), names = val)
      }
    }
    cols <- pal[x]
    cols[is.na(cols)] <- na.value
    return(cols)
  } else {
    breaks <- pretty(x, n = n + 1)
    n2 <- length(breaks) - 1
    # order = findInterval(x, sort(x))
    if (is.null(pal)) {
      cols <- viridis::viridis(n = n2, ...) # [order]
    } else {
      cols <- pal(n = n2, ...) # [order]
    }

    ret <- cut(x, breaks = breaks, labels = cols, include.lowest = TRUE) |>
      as.character()
    names(ret) <- cut(x, breaks = breaks, include.lowest = TRUE)
    return(ret)
  }
}


#' Color palette for stress regime
#'
#' @return function
#' @export
#'
#' @examples
#' stress_colors()
stress_colors <- function() {
  structure(
    c("#D55E00", "#E69F00", "#009E73", "#56B4E9", "#0072B2", "grey60"),
    names = c("N", "NS", "S", "TS", "T", "U")
  )
}


#' Plot axes
#'
#' @param x,y coordinates of points
#' @param angle Azimuth in degrees
#' @param radius length of axis
#' @param arrow.code integer. Kind of arrow head. The default is `1`, i.e. no
#' arrow head. See [graphics::arrows()] for details
#' @param arrow.length numeric Length of the edges of the arrow head (in
#' inches). (Ignored if `arrow.code = 1`)
#' @param add logical. add to existing plot?
#' @param ... optional arguments passed to [graphics::arrows()]
#'
#' @returns No return value, called for side effects
#'
#' @export
#'
#' @examples
#' data("san_andreas")
#' axes(san_andreas$lon, san_andreas$lat, san_andreas$azi, add = FALSE)
axes <- function(x, y, angle, radius = .5, arrow.code = 1, arrow.length = 0, add = FALSE, ...) {
  if (!add) plot(x, y, cex = 0, ...)
  graphics::arrows(x, y, x1 = x + radius / 2 * cosd(270 - angle), y1 = y + radius / 2 * sind(270 - angle), length = arrow.length, code = arrow.code, ...)
  graphics::arrows(x, y, x1 = x + radius / 2 * cosd(90 - angle), y1 = y + radius / 2 * sind(90 - angle), length = arrow.length, code = arrow.code, ...)
}

#' Class for Central Position of Spoke Marker
#'
#' position subclass \code{"center_spoke"} to center \code{ggplot::geom_spoke()}
#' marker at its origin
#'
#' @noRd
position_center_spoke <- function() PositionCenterSpoke #

#' @title  Centrically aligned geom_spoke marker
#'
#' @description \code{"position"} subclass "center_spoke" to center
#' \code{ggplot::geom_spoke()} marker at its origin
#'
#' @export
#' @keywords internal
#' @source \url{https://stackoverflow.com/questions/55474143/how-to-center-geom-spoke-around-their-origin}
#'
#' @importFrom ggplot2 ggproto Position
PositionCenterSpoke <- ggplot2::ggproto("PositionCenterSpoke", ggplot2::Position,
  compute_panel = function(self, data, params, scales) {
    data$x <- 2 * data$x - data$xend
    data$y <- 2 * data$y - data$yend
    data$radius <- 2 * data$radius
    data
  }
)

#' Quantile-Quantile Linearised Plot for Circular Distributions
#'
#'
#' Uniformly distributed orientations should yield a straight line through the
#' origin. Systematic departures from linearity will indicate preferred
#' orientation in some manner.
#'
#' @param x numeric. Angles in degrees
#' @param axial Logical. Whether data are uniaxial (`axial=FALSE`)
#' @param xlab,ylab,main plot labels.
#' @param ... graphical parameters
#'
#' @return plot
#' @importFrom graphics plot abline lines points
#' @export
#'
#' @references Borradaile, G. J. (2003). Statistics of earth
#' science data: their distribution in time, space, and orientation (Vol. 351,
#' p. 329). Berlin: Springer.
#'
#' @examples
#' x_vm <- rvm(100, mean = 0, kappa = 2)
#' circular_qqplot(x_vm, pch = 20)
#'
#' x_norm <- rnorm(100, mean = 0, sd = 25)
#' circular_qqplot(x_norm, pch = 20)
#'
#' x_unif <- runif(100, 0, 360)
#' circular_qqplot(x_unif, pch = 20)
circular_qqplot <- function(x, axial = TRUE,
                            xlab = paste("i/(n+1)"),
                            ylab = NULL, main = "Circular Quantile-Quantile Plot", ...) {
  if (axial) {
    f <- 2
  } else {
    f <- 1
  }
  if (is.null(ylab)) {
    ylab <- paste0(deparse1(substitute(x)), "/", 360 / f)
  }

  x <- (x %% (360 / f)) / (360 / f)
  x <- sort(x)
  n <- length(x)
  xin <- seq_along(x) / (n + 1)

  graphics::plot(xin, x,
    type = "p", xlim = c(0, 1), ylim = c(0, 1), asp = 1,
    xlab = xlab, ylab = ylab, main = main,
    sub = bquote("N" == .(n)),
    ...
  )
  graphics::abline(a = 0, b = 1, col = "slategrey")
  graphics::lines(xin, x)
  invisible(xin)
  # graphics::points(xin, x, col = "slategrey")
}

#' von Mises Quantile-Quantile Plot
#'
#' Produces a Q-Q plot of the data against a specified von Mises distribution
#' to graphically assess the goodness of fit of the model.
#'
#' @param x numeric. Angles in degrees
#' @param w numeric. optional weightings for `x` to estimate `mean` and `kappa`.
#' @param axial Logical. Whether data are uniaxial (`axial=FALSE`)
#' @param xlab,ylab,main plot labels.
#' @param mean numeric. Circular mean of the von Mises distribution. If `NULL`,
#' it will be estimated from `x`.
#' @param kappa numeric. Concentration parameter of the von Mises distribution.
#' If `NULL`, it will be estimated from `x`.
#' @param ... graphical parameters
#'
#' @return plot
#' @importFrom stats ecdf
#' @importFrom graphics plot lines
#'
#' @export
#'
#' @examples
#' x_vm <- rvm(100, mean = 0, kappa = 4)
#' vm_qqplot(x_vm, axial = FALSE, pch = 20)
#'
#' x_unif <- runif(100, 0, 360)
#' vm_qqplot(x_unif, axial = FALSE, pch = 20)
vm_qqplot <- function(x, w = NULL, axial = TRUE, mean = NULL, kappa = NULL,
                      xlab = "von Mises quantile function",
                      ylab = "Empirical quantile function",
                      main = "von Mises Q-Q Plot", ...) {
  if (axial) {
    f <- 2
  } else {
    f <- 1
  }

  n <- length(x)

  # if (stretch) {
  #   k <- 4
  # } else {
  #   k <- 2
  # }

  if (is.null(mean)) mean <- circular_mean(x, w = w, axial = axial)
  if (is.null(kappa)) kappa <- est.kappa(x, w = w, axial = axial)

  caption <- bquote(
    bar(alpha) == .(round(mean, 1)) * degree ~ "|" ~ kappa == .(round(kappa, 1))
  )

  xf <- (x * f) %% 360
  xf <- sort(xf)

  edf <- stats::ecdf(xf)


  # z <- sind((xf - mean*f) / k)
  #
  # z_sort <- sort(z) #|> scales::rescale(c(-1, 1))
  #
  #
  # probs <- seq(0, 1, length.out = n)
  # quantiles <- qvm(edf(xf), mean = 0, kappa = kappa, from = 0)
  # quantiles <- qvm(probs, mean = 0, kappa = kappa, from = 0)

  #
  # sin_q <- sind(quantiles / k)

  # graphics::plot(0, 0, type = "n", xlim = c(-1, 1), ylim = c(-1, 1), asp = 1, xlab = xlab, ylab = ylab, main = main, ...)
  # graphics::abline(a = 0, b = 1, col = "slategrey")
  # graphics::points(sin_q, z_sort, pch=20)

  tqf <- qvm(edf(xf), mean * f, kappa, from = 0)
  #
  graphics::plot(tqf / f, xf / f, xlim = c(0, 360 / f), ylim = c(0, 360 / f), asp = 1, xlab = xlab, ylab = ylab, main = main, sub = bquote("N" == .(n)), ...)
  graphics::abline(a = 0, b = 1, col = "slategrey")
  mtext(caption)
  invisible(tqf)
}


# Circular diagram -----------------------------------------------------------------


#' Circular plot
#'
#' @param main Character string specifying the title of the plot.
#' @param at Optional vector of angles at which tick marks should be plotted.
#' Set `at=numeric(0)` to suppress tick marks.
#' @param labels Either a logical value indicating whether to plot labels
#' next to the tick marks, or a vector of labels for the tick marks.
#' @param cborder logical. Border of rose plot.
#' @param ... optional arguments passed to `plot.default()`
#' @importFrom spatstat.geom disc plot.owin
#' @import spatstat.utils
#' @import spatstat.explore
#' @importFrom spatstat.univar whist
#'
#' @note Polar diagram where angles increase clockwise.
#'
#' @export
#'
#' @return none
#'
#' @keywords internal
circular_plot <- function(main = NULL, labels = TRUE,
                          at = seq(0, 360 - 45, 45), cborder = TRUE, ...) {
  unit <- "degree"
  ymax <- 1
  insideclearance <- 0.1
  outsidespace <- if (!is.null(at) && length(at) == 0) {
    0
  } else if (identical(labels, FALSE)) {
    0.1
  } else {
    0.25
  }
  R <- (1 + insideclearance) * ymax
  DD <- spatstat.geom::disc(R)
  Rout <- (1 + outsidespace) * R
  disco <- spatstat.geom::disc(Rout)
  spatstat.utils::dont.complain.about(DD, disco)
  result <- spatstat.utils::do.call.matched(
    spatstat.geom::plot.owin, spatstat.utils::resolve.defaults(list(
      x = quote(disco),
      main = main, type = "n", border = cborder
    ), list(...))
  )
  spatstat.utils::do.call.matched(
    spatstat.geom::plot.owin, spatstat.utils::resolve.defaults(
      list(
        x = quote(DD),
        hatch = FALSE, add = TRUE, border = cborder
      ), list(...)
    ),
    extrargs = spatstat.utils::graphicsPars("owin"), skipargs = "col"
  )
  spatstat.explore::circticks(
    R,
    at = at, unit = unit, start = -90, clockwise = TRUE, labels = labels
  )
  return(invisible(result))
}


rose_histogram <- function(x, ..., main, labels = TRUE, at = NULL,
                           cborder = TRUE, axial = FALSE, add = FALSE) {
  if (missing(main) || is.null(main)) {
    main <- spatstat.utils::short.deparse(substitute(x))
  }

  bks <- x$mids
  bw <- x$binwidth
  y <- x$density

  if (!add) circular_plot(smain = main, labels, at = at, cborder = cborder)

  for (i in seq_along(y)) {
    rose_fan(bks[i], d = bw, radius = y[i], axial = axial, add = TRUE, ...)
  }
}

graphicsAargh <- c(
  "density", "angle", "col", "border",
  "xlim", "ylim", "xlab", "ylab", "axes"
)

rose_freq <- function(x, bins = NULL, ..., weights = NULL, binwidth = NULL,
                      round_binwidth = 0, equal_area = TRUE,
                      main = NULL, axial = TRUE) {
  if (missing(main) || is.null(main)) {
    main <- spatstat.utils::short.deparse(substitute(x))
  }

  stopifnot(is.numeric(x))
  if (!is.null(weights)) {
    spatstat.utils::check.nvector(weights, length(x),
      things = "observations",
      vname = "weights"
    )
  }

  if (axial) {
    if (!is.null(bins) && is.null(binwidth)) {
      stopifnot(bins > 0)
      binwidth <- 180 / bins # bin width
    } else if (is.null(bins) && is.null(binwidth)) {
      binwidth <- rose_binwidth(length(stats::na.omit(x)), axial = axial)
    }

    binwidth <- round(binwidth)
    stopifnot(binwidth > 0)

    binwidth <- symmetric_bw(binwidth)

    x <- x %% 180
    breaks <- seq(0, 180, binwidth) |> add_end(180)
  } else {
    if (!is.null(bins) && is.null(binwidth)) {
      stopifnot(bins > 0)
      binwidth <- 360 / bins # bin width
    } else if (is.null(bins) && is.null(binwidth)) {
      binwidth <- rose_binwidth(length(stats::na.omit(x)), axial = axial)
    }

    binwidth <- round(binwidth)
    stopifnot(binwidth > 0)

    breaks <- seq(0, 360, binwidth) |> add_end(360)
  }

  h <- spatstat.utils::do.call.matched(graphics::hist.default, list(
    x = x, breaks = breaks, ..., plot = FALSE
  ),
  skipargs = graphicsAargh, sieve = TRUE
  )

  result <- numeric()
  result <- h$result
  result$otherargs <- h$otherargs

  freqs <- spatstat.univar::whist(
    x = x, breaks = breaks, weights = weights
  )

  result$count <- freqs
  result$density <- freqs / binwidth

  if (equal_area) {
    result$density <- sqrt(result$density)
  }
  result$xname <- main
  result$density <- result$density / max(result$density)
  result$binwidth <- binwidth

  result
}



#' @title Selecting optimal number of bins and width for rose diagrams
#'
#' @param n Integer. number of data
#' @param round Logical. Whether bin width is round to zero digits
#' (`round=TRUE`, the default)
#' or as is (`FALSE`).
#' @param axial Logical. Whether data are uniaxial (`axial=FALSE`)
#' or biaxial (`TRUE`, the default).
#' @param ... Additional arguments passed to [rose_bw()].
#'
#' @keywords internal
#' @name rose_bw
NULL

#' @rdname rose_bw
rose_bins <- function(n, round = FALSE) {
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
  r / rose_bins(n)
}

is.naturalnumber <- function(x, tol = .Machine$double.eps^0.5) x > tol & abs(x - round(x)) < tol

symmetric_bw <- function(x) {
  div <- 180 / seq(1, 180, 1)
  cond <- is.naturalnumber(div) & div < 180
  allowed <- div[cond]
  target.index <- which(abs(allowed - x) == min(abs(allowed - x)))
  min(allowed[target.index])
}

add_end <- function(x, end) {
  check <- end %in% x
  if (check) {
    x
  } else {
    x <- c(x, end)
  }
}

#' @title Rose Diagram
#'
#' @description Plots a rose diagram (rose of directions), the analogue of a
#' histogram or density plot for angular data.
#'
#' @param x Data to be plotted. A numeric vector containing angles (in degrees).
#' @param weights Optional vector of numeric weights associated with x.
#' @param binwidth The width of the bins (in degrees).
#' @param bins number of arcs to partition the circle width.
#' Overridden by `binwidth`.
#' @param axial Logical. Whether data are uniaxial (`axial=FALSE`)
#' or biaxial (`TRUE`, the default).
#' @param equal_area Logical. Whether the radii of the bins are proportional to
#' the frequencies (`equal_area=FALSE`, i.e. equal-angle) or proportional to the
#' square-root of the frequencies (`equal_area=TRUE`, the default).
#' @param muci logical. Whether the mean and its 95% CI are added to the plot
#' or not.
#' @param round_binwidth integer. Number of decimal places of bin width (0 by
#' default).
#' @param mtext character. String to be drawn at the top margin of the plot
#' (`"N"` by default)
#' @param main,sub Character string specifying the title and subtitle of the
#' plot. If `sub = NULL`, it will show the bin width.
#' @inheritParams circular_plot
#' @param col fill color of bins
#' @param dots logical. Whether a circular dot plot should be added
#' (`FALSE` is the default).
#' @param jitter_factor Add a small amount of noise to the angles' radius that
#' is added to `scale`. Jitter is ignored when `stack==TRUE`).
#' If `0`, no jitter is added (by default); if negative, the points fall into
#' the circle.
#' @param stack logical. Groups and stacks the dots if `TRUE`. Default is `FALSE`.
#' @param dot_cex,dot_pch,dot_col Plotting arguments for circular dot plot
#' @param add logical.
#' @param ... Additional arguments passed to [spatstat.explore::rose()].
#'
#' @note If `bins` and `binwidth` are `NULL`, an optimal bin width will be
#' calculated using Scott (1979):
#' \deqn{ w_b = \frac{R}{n^{\frac{1}{3}}}
#' }
#' with n being the length of `x`, and the range R being either 180 or 360
#' degree for axial or directional data, respectively.
#'
#' If `"axial" == TRUE`, the binwidth is adjusted to guarantee symmetrical fans.
#'
#' @return A window (class `"owin"`) containing the plotted region or a `list`
#' of the calculated frequencies.
#'
#' @importFrom spatstat.explore rose
#' @importFrom spatstat.utils short.deparse
#' @importFrom graphics hist title points mtext
#' @importFrom stats na.omit
#'
#' @export
#'
#' @examples
#' x <- rvm(100, mean = 90, k = 5)
#' rose(x, axial = FALSE, border = TRUE)
#'
#' data("san_andreas") #'
#' rose(san_andreas$azi, main = "equal area")
#' rose(san_andreas$azi, equal_area = FALSE, main = "equal angle")
#'
#' # weighted frequencies:
#' rose(san_andreas$azi, weights = 1 / san_andreas$unc, main = "weighted")
#'
#' # add dots
#' rose(san_andreas$azi, dots = TRUE, main = "dot plot", jitter = .2)
#' rose(san_andreas$azi,
#'   dots = TRUE, stack = TRUE, dot_cex = 0.5, dot_pch = 21,
#'   main = "stacked dot plot"
#' )
rose <- function(x, weights = NULL, binwidth = NULL, bins = NULL, axial = TRUE,
                 equal_area = TRUE, muci = TRUE,
                 round_binwidth = 0, mtext = "N", main = NULL, sub = NULL,
                 at = seq(0, 360 - 45, 45), cborder = TRUE, labels = TRUE,
                 col = "grey", dots = FALSE, dot_pch = 1, dot_cex = 1,
                 dot_col = "slategrey", stack = FALSE, jitter_factor = 0,
                 add = FALSE, ...) {
  if (!add) {
    if (missing(main) || is.null(main)) {
      main <- spatstat.utils::short.deparse(substitute(x))
    }
    circular_plot(main = main, labels = labels, at = at, cborder = cborder)
  }

  if (axial) {
    x <- x %% 180
    x[x >= 180] <- 180 - 2 * .Machine$double.eps
  } else {
    x <- x %% 360
    x[x >= 360] <- 360 - 4 * .Machine$double.eps
  }

  freqs <- rose_freq(
    x,
    bins = bins, ..., weights = weights, binwidth = binwidth,
    round_binwidth = round_binwidth, equal_area = equal_area,
    axial = axial
  )

  rose_histogram(freqs, ...,
    col = col, axial = axial,
    main = main, labels = TRUE, at = at, cborder = TRUE, add = TRUE
  )

  if (dots) {
    plot_points(x,
      axial = axial, stack = stack, cex = dot_cex, pch = dot_pch,
      col = dot_col, jitter_factor = jitter_factor, add = TRUE
    )
  }

  if (is.null(sub)) sub <- paste("Bin width:", freqs$binwidth)
  graphics::title(main = NULL, sub = sub, ylab = NULL)
  graphics::mtext(mtext)

  if (muci) rose_stats(x, weights = weights, axial = axial)
  invisible(freqs)
}


#' Direction Lines and Fans in Circular Diagram
#'
#' @param x angles in degrees
#' @param d width of a fan (in degrees)
#' @param radius of the plotted circle
#' @param axial Logical. Whether `x` are uniaxial (`axial=FALSE`)
#' or biaxial (`TRUE`, the default).
#' @param add logical. Add to existing plot?
#' @param ... optional arguments passed to [graphics::segments()] or
#' [graphics::polygon()]
#'
#' @returns No return value, called for side effects
#'
#' @importFrom graphics segments polygon
#' @name rose_geom
#' @examples
#' angles <- c(0, 10, 45)
#' radius <- c(.7, 1, .2)
#' lwd <- c(2, 1, .75)
#' col <- c(1, 2, 3)
#' rose_line(angles, radius = radius, axial = FALSE, add = FALSE, lwd = lwd, col = col)
NULL

#' @rdname rose_geom
#' @export
rose_line <- function(x, radius = 1, axial = TRUE, add = TRUE, ...) {
  xrad <- deg2rad(90 - x)
  tx <- radius * cos(xrad)
  ty <- radius * sin(xrad)

  if (!add) circular_plot()
  graphics::segments(0, 0, tx, ty, ...)
  if (axial) {
    graphics::segments(0, 0, -tx, -ty, ...)
  }
  invisible()
}

#' @rdname rose_geom
#' @export
rose_fan <- function(x, d, radius = 1, axial = TRUE, add = TRUE, ...) {
  xrad <- deg2rad(x)
  drad <- deg2rad(d) / 2

  eps <- (pi / 128) / 2
  aa <- (pi / 2) - seq(xrad - drad, xrad + drad, by = eps)

  tx <- radius * cos(aa)
  ty <- radius * sin(aa)
  xx <- c(0, tx, 0)
  yy <- c(0, ty, 0)

  if (!add) circular_plot()

  graphics::polygon(x = xx, y = yy, ...)
  if (axial) {
    graphics::polygon(x = -xx, y = -yy, ...)
  }
}

#' Show Average Direction and Spread in Rose Diagram
#'
#' Adds the average direction (and its spread) to an existing rose diagram.
#'
#' @param x Data to be plotted. A numeric vector containing angles (in degrees).
#' @param weights Optional vector of numeric weights associated with x.
#' @param axial Logical. Whether data are uniaxial (`axial=FALSE`)
#' or biaxial (`TRUE`, the default).
#' @param avg character. The average estimate for x. Either the circular mean
#' (`"mean"`, the default), the circular Quasi Median (`"median"`), or the
#' sample median (`"sample_median"`).
#' @param spread character. The measure of spread to be plotted as a fan.
#' Either 95% confidence interval (`"CI"`, the default), Fishers confidence interval (`"fisher"`), the circular
#' standard deviation (`"sd"`), the Quasi interquartile range on the circle
#' (`"IQR"`), or the sampke median deviation (`"mdev"`). `NULL` if no fan should be drawn.
#' @param avg.col color for the average line
#' @param avg.lty line type of the average line
#' @param avg.lwd  line width of the average line
#' @param spread.col color of the spread fan
#' @param spread.border logical. Whether to draw a border of the fan or not.
#' @param spread.lty line type of the spread fan's border
#' @param spread.lwd line width of the spread fan's border
#' @param add logical.
#' @param ... optional arguments to `circular_plot()` if add is `FALSE`.
#' @importFrom ggplot2 alpha
#'
#' @seealso [rose()] for plotting the rose diagram, and
#' [circular_mean()], [circular_median()], [circular_sample_median()],
#' [confidence_interval()], [confidence_interval_fisher()],
#' [circular_sd()], [circular_IQR()], [circular_sample_median_deviation()]
#' for statistical parameters.
#'
#' @returns No return value, called for side effects
#' @export
#'
#' @examples
#' data("san_andreas")
#' rose(san_andreas$azi, weights = 1 / san_andreas$unc, muci = FALSE)
#' rose_stats(san_andreas$azi, weights = 1 / san_andreas$unc, avg = "sample_median", spread = "mdev")
rose_stats <- function(x, weights = NULL, axial = TRUE, avg = c("mean", "median", "sample_median"), spread = c("CI", "fisher", "sd", "IQR", "mdev"),
                       avg.col = "#85112AFF", avg.lty = 2, avg.lwd = 1.5,
                       spread.col = ggplot2::alpha("#85112AFF", .2), spread.border = FALSE, spread.lty = NULL, spread.lwd = NULL, add = TRUE, ...) {
  avg <- match.arg(avg)
  mu <- switch(avg,
    mean = circular_mean(x, weights, axial),
    median = circular_median(x, weights, axial),
    sample_median = circular_sample_median(x, axial)
  )
  # mu_text <- switch(avg,
  #   mean = "Mean: ",
  #   median = "Median: "
  # )

  if (!is.null(spread)) {
    spread <- match.arg(spread)
    ci <- switch(spread,
      CI = confidence_interval(x, w = weights, axial = axial)$conf.angle,
      fisher = confidence_interval_fisher(x, w = weights, axial = axial, quiet = TRUE)$conf.angle,
      sd = circular_sd(x, weights, axial),
      IQR = circular_IQR(x, weights, axial),
      mdev = circular_sample_median_deviation(x, axial)
    )
    rose_fan(mu, ci,
      radius = 1.1, axial = axial, col = spread.col,
      border = spread.border, lty = spread.lty, lwd = spread.lwd,
      add = add, ...
    )
  }

  rose_line(mu,
    radius = 1.1, axial = axial, col = avg.col, lty = avg.lty,
    lwd = avg.lwd, add = TRUE, ...
  )
  invisible(c(mu, ci))
}

# Dot plot ---------------------------------------------------------------------


#' Add Points to a Circular Plot
#'
#' Add points to a plot of circular data points on the current graphics device.
#'
#' @param x Data to be plotted. A numeric vector containing angles (in degrees).
#' @param axial Logical. Whether data are uniaxial (`axial=FALSE`)
#' or biaxial (`TRUE`, the default).
#' @param stack logical: if `TRUE`, points are stacked on the perimeter of the
#' circle.
#' Otherwise, all points are plotted on the perimeter of the circle. Default is
#' `FALSE`.
#' @param binwidth numeric. Bin width (in degrees) for the stacked dot plots.
#' ignored when `stack==FALSE`. Is set to `1` degree by default.
#' @param cex character (or symbol) expansion: a numerical vector. This works as
#' a multiple of `par("cex")`.
#' @param sep constant used to specify the distance between stacked points, if
#' `stack==TRUE` or in the case of more than one dataset. Default is `0.025`;
#' smaller values will create smaller spaces.
#' @param jitter_factor numeric. Adds a small amount of random variation to the
#' location of each points along radius that is added to `scale`. Jitter is
#' ignored when `stack==TRUE`). If `0`, no jitter is added (by default); if
#' negative, the points fall into the circle.
#' @param ... Further graphical parameters may also be supplied as arguments.
#' @param scale radius of plotted circle. Default is `1.1`.
#' Larger values shrink the circle, while smaller values enlarge the circle.
#' @param add logical
#' @inheritParams circular_plot
#'
#' @importFrom graphics points
#'
#' @return A list with information on the plot
#' @export
#'
#' @examples
#' x <- rvm(100, mean = 90, k = 5)
#' plot_points(x, add = FALSE)
#' plot_points(x, jitter_factor = .2, add = FALSE) # jittered plot
#' plot_points(x, stack = TRUE, binwidth = 3, add = FALSE) # stacked plot
plot_points <- function(x, axial = TRUE, stack = FALSE, binwidth = 1, cex = 1, sep = 0.025, jitter_factor = 0, ..., scale = 1.1, add = TRUE,
                        main = NULL, labels = TRUE,
                        at = seq(0, 360 - 45, 45), cborder = TRUE) {
  if (!add) {
    if (missing(main) || is.null(main)) {
      main <- spatstat.utils::short.deparse(substitute(x))
    }
    circular_plot(main = main, labels = labels, at = at, cborder = cborder)
  }

  f <- as.numeric(axial) + 1

  if (!stack) {
    if (axial) {
      x_shift <- (x + 180) %% 360
      x <- c(x, x_shift)
    }
    xr <- deg2rad(x)
    u <- pi / 2 - xr
    n <- length(x)
    r <- scale + runif(n, min(0, jitter_factor), max(0, jitter_factor))
    z <- cos(u) * r
    y <- sin(u) * r
    graphics::points(z, y, cex = cex, ...)
  } else {
    freqs <- rose_freq(x, axial = axial, binwidth = binwidth)
    if (axial) {
      freqs$mids <- freqs$mids %% 180
      freqs$count <- rep(freqs$count, 2)
      freqs$mids <- c(freqs$mids, freqs$mids + 180)
    }

    bins <- f * 180 / freqs$binwidth
    bins.count <- freqs$count
    mids <- deg2rad(-90 - freqs$mids)
    index <- cex * sep
    # index <- cex * freqs$binwidth

    for (i in 1:bins) {
      if (bins.count[i] > 0) {
        for (j in 0:(bins.count[i] - 1)) {
          r <- scale + j * index
          z <- r * cos(mids[i])
          y <- r * sin(mids[i])
          graphics::points(z, y, cex = cex, ...)
        }
      }
    }
  }
}

# Plot density lines on rose ---------------------------------------------------

calc_circular_density <- function(x, z, kappa) {
  nx <- length(x)
  # if (kernel == "vonmises") {
  y <- sapply(z, dvm, mean = x, kappa = kappa)
  # }
  # else if (kernel == "wrappednormal") {
  #   rho <- exp(-bw^2/2)
  #   y <- sapply(z, DwrappednormalRad, mu = x, rho = rho,
  #               K = K, min.k = min.k)
  # }
  # else {
  #   stop("other kernels not implemented yet")
  # }
  apply(y, 2, sum) / nx
}


circular_density <- function(x, z = NULL, kappa, na.rm = TRUE, from = 0, to = 360, n = 512, axial = TRUE) {
  f <- as.numeric(axial) + 1
  x <- x * f

  if (is.null(z)) {
    z <- seq(from = from, to = to, length = n)
  } else {
    if (!is.numeric(z)) {
      stop("argument 'z' must be numeric")
    }
    namez <- deparse(substitute(z))
    z.na <- is.na(z)
    if (any(z.na)) {
      if (na.rm) {
        z <- z[!z.na]
      } else {
        stop("z contains missing values")
      }
    }
    z.finite <- is.finite(z)
    if (any(!z.finite)) {
      z <- z[z.finite]
    }
  }

  calc_circular_density(x, z, kappa)
}

circular_lines <- function(x, y, join = FALSE, nosort = FALSE, offset = 1.1, shrink = 1, axial = TRUE, ...) {
  x <- deg2rad(90 - x)

  if (axial) {
    x <- c(x, x + pi)
    y <- rep(y, 2)
  }


  n <- length(x)
  if (!nosort) {
    xorder <- order(x)
    x <- x[xorder]
    y <- y[xorder]
    spacings <- c(diff(x), x[1] - x[n] + 2 * pi)
    pos <- which.max(spacings)[1]
    if (pos == n) {
      xorder <- 1:n
    } else {
      xorder <- c((pos + 1):n, 1:pos)
    }
  } else {
    xorder <- 1:n
  }
  z <- (y / shrink + offset) * cos(x)
  w <- (y / shrink + offset) * sin(x)
  z <- z[xorder]
  w <- w[xorder]
  if (join) {
    z <- c(z, z[1])
    w <- c(w, w[1])
  }
  graphics::lines(x = z, y = w, ...)
  invisible(list(x = z, y = w))
}


#' Circular density plot
#'
#' Plot the multiples of a von Mises density distribution
#'
#' @param x Data to be plotted. A numeric vector containing angles (in degrees).
#' @param kappa Concentration parameter for the von Mises distribution.
#' Small kappa gives smooth density lines.
#' @param axial Logical. Whether data are uniaxial (`axial=FALSE`)
#' or biaxial (`TRUE`, the default).
#' @param n the number of equally spaced points at which the density is to be estimated.
#' @param norm_density logical. Normalize the density?
#' @param scale radius of plotted circle. Default is `1.1`.
#' @param shrink parameter that controls the size of the plotted function. Default is 1.
#' @param ... Further graphical parameters may also be supplied as arguments.
#' @param add logical. Add to existing plot? (`TRUE` by default).
#' @inheritParams circular_plot
#'
#' @seealso [dvm()]
#' @return plot or calculated densities as numeric vector
#' @export
#'
#' @examples
#' rose(san_andreas$azi, dots = TRUE, stack = TRUE, dot_cex = 0.5, dot_pch = 21)
#' plot_density(san_andreas$azi, kappa = 10, col = "seagreen", shrink = 1.5)
#' plot_density(san_andreas$azi, kappa = 10, col = "seagreen", add = FALSE, scale = .6)
plot_density <- function(x, kappa, axial = TRUE, n = 512, norm_density = TRUE, ...,
                         scale = 1.1, shrink = 1,
                         add = TRUE, main = NULL, labels = TRUE,
                         at = seq(0, 360 - 45, 45), cborder = TRUE) {
  if (!add) {
    if (missing(main) || is.null(main)) {
      main <- spatstat.utils::short.deparse(substitute(x))
    }
    circular_plot(main = main, labels = labels, at = at, cborder = cborder)
  }


  f <- as.numeric(axial) + 1
  d <- circular_density(x, kappa = kappa, n = n, axial = axial)
  if (norm_density) d / max(d)
  circular_lines(seq(0, 360, length = f * n), rep(d, f), axial = FALSE, n, offset = scale, shrink = shrink, ...)
  invisible(d)
}




#' Plotting stress analysis results
#'
#' Creates a set of plots including
#' the azimuth as a function of the distance to the plate boundary,
#' the Norm Chi-squared as a function of the distance to the plate boundary,
#' the circular distance (and dispersion) a function of the distance to the
#' plate boundary, and a rose diagram of the frequency distribution of the
#' azimuths.
#'
#' @param azi numeric. Azimuth of \eqn{\sigma_{Hmax}}{SHmax}
#' @param distance numeric. Distance to plate boundary
#' @param prd numeric. the predicted direction of \eqn{\sigma_{Hmax}}{SHmax}
#' @param unc numeric. Uncertainty of observed \eqn{\sigma_{Hmax}}{SHmax},
#' either a numeric vector or a number
#' @param regime character vector. The stress
#' regime (following the classification of the World Stress Map)
#' @param width integer. window width (in number of observations) for moving
#' average of the azimuths, circular dispersion, and Norm Chi-square statistics.
#' If `NULL`, an optimal width will be estimated.
#'
#' @importFrom dplyr arrange mutate
#'
#' @seealso [PoR_shmax()], [distance_from_pb()], [circular_mean()],
#' [circular_dispersion()], [confidence_interval_fisher()], [norm_chisq()],
#' [weighted_rayleigh()]
#'
#' @details
#' Plot 1 shows the transformed azimuths as a function of the distance to the
#' plate boundary. The red line indicates the rolling circular mean, stippled
#' red lines indicate the 95% confidence interval about the mean.
#'
#' Plot 2 shows the normalized \eqn{\chi^2}{chi-squared} statistics as a
#' function of the distance to the plate boundary. The red line shows the
#' rolling \eqn{\chi^2}{chi-squared} statistic.
#'
#' Plot 3 shows the circular distance of the transformed azimuths to the
#' predicted azimuth, as a function of the distance to the plate boundary. The
#' red line shows the rolling circular dispersion about the prediction.
#'
#' Plot 4 give the rose diagram of the transformed azimuths.
#'
#' @returns four R base plots
#'
#' @export
#'
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
#' quick_plot(res$azi.PoR, d, res$prd, san_andreas$unc, san_andreas$regime)
quick_plot <- function(
    azi,
    distance,
    prd,
    unc = NULL,
    regime,
    width = 51) {
  if (missing(regime)) {
    regime <- rep(NA, length(azi))
  }
  nchisq_i <- numeric()
  regime <- ifelse(is.na(regime), "U", regime)

  t <- data.frame(azi, distance, prd, unc,
    regime = factor(regime, levels = c("U", "N", "NS", "S", "TS", "T"))
  ) |>
    dplyr::filter(!is.na(azi)) |>
    dplyr::arrange(distance) |>
    dplyr::mutate(
      nchisq_i = (deviation_norm(azi, prd) / unc)^2 / (90 / unc)^2,
      cdist = circular_distance(azi, prd),
      roll_mean = roll_circstats(
        azi,
        w = 1 / unc,
        FUN = circular_mean,
        width = width
      ),
      roll_sd = roll_circstats(
        azi,
        w = 1 / unc,
        FUN = circular_sd,
        width = width
      ) / 2,
      roll_nchisq = roll_normchisq(
        azi,
        prd,
        unc,
        width = width
      ),
      roll_disp = roll_dispersion(
        azi,
        prd,
        w = 1 / unc,
        width = width
      )
    )

  # add lower and upper period to data for plotting
  tmin <- dplyr::mutate(t, azi = azi - 180)
  tmax <- dplyr::mutate(t, azi = azi + 180)
  t2 <- rbind(tmin, t, tmax)

  nchisq <- norm_chisq(azi, prd, unc)
  rt <- weighted_rayleigh(azi, mu = prd, w = 1 / unc, quiet = TRUE)
  azi.PoR.mean <- circular_mean(azi, 1 / unc)
  azi.PoR.sd <- circular_sd(azi, 1 / unc)
  disp <- circular_dispersion(azi, prd, 1 / unc)
  CIF <- confidence_interval_fisher(azi, w = 1 / unc, quiet = TRUE)
  CI <- CIF$conf.interval
  CI_ang <- CIF$conf.angle

  subtitle <-
    bquote(95 * "% CI [" * .(round(CI[1])) * degree * "," ~ .(round(CI[2])) * degree * "] | R" == .(signif(rt$statistic, 2)) ~ ("p" == .(signif(rt$p.value, 2))))

  subtitle_rose <- bquote(atop(
    "N" == .(length(azi)),
    bar(alpha) == .(round(azi.PoR.mean, 1)) * degree * "" %+-% "" * .(round(CI_ang, 1)) * ~degree
  ))
  # subtitle_rose <- do.call(expression, subtitle_rose)

  grDevices::palette(c("grey60", "#D55E00", "#E69F00", "#009E73", "#56B4E9", "#0072B2"))

  # distance plot
  ## create empty plot
  graphics::plot(
    0,
    type = "n",
    xlab = "Distance from plate boundary",
    ylab = expression("Azimuth wrt. PoR " ~ alpha ~ "(" * degree * ")"),
    sub = subtitle,
    main = "Distance from plate boundary vs. azimuth",
    xlim = range(distance),
    ylim = c(0, 180),
    yaxp = c(0, 180, 8)
  )

  ## 95% confidence interval
  graphics::polygon(
    x = c(rev(t$distance), t$distance),
    y = c(rev(t$roll_mean + t$roll_sd), t$roll_mean - t$roll_sd),
    col = "grey90", border = FALSE
  )

  ## points
  graphics::arrows(
    y0 = t2$azi - t2$unc, x0 = t2$distance,
    y1 = t2$azi + t2$unc, x1 = t2$distance,
    code = 0, lwd = .25, col = t2$regime
  )
  graphics::points(azi ~ distance, data = t2, col = t2$regime)

  ## roll statistics
  graphics::lines(roll_mean ~ distance, data = t, type = "S", col = "#85112AFF")

  ## predicted az
  graphics::abline(h = unique(prd), col = "black", lty = 2)
  graphics::legend("bottomright",
    inset = .05, cex = .75,
    legend = names(stress_colors()), title = "Stress regime",
    fill = stress_colors()
  )

  # Norm chisq plot
  grDevices::dev.new()
  graphics::plot(nchisq_i ~ distance,
    data = t, col = t$regime,
    xlab = "Distance from plate boundary", ylab = expression("Norm" ~ chi[i]^2),
    main = expression(bold("Deviation from prediction" ~ beta)),
    xlim = range(distance),
    ylim = c(0, 1), yaxp = c(0, 1, 4),
    sub = bquote("Norm" ~ chi^2 == .(round(nchisq, 2)))
  )
  graphics::lines(
    roll_nchisq ~ distance,
    data = t,
    type = "S",
    col = "#85112AFF"
  )
  graphics::abline(h = .15, col = "black", lty = 2)

  # Dispersion plot
  grDevices::dev.new()
  graphics::plot(0,
    type = "n",
    xlab = "Distance from plate boundary", ylab = expression("Circular distance " ~ "d(" * alpha[i] * "," ~ beta * ")"),
    main = expression(bold("Circular dispersion around prediction" ~ beta)),
    xlim = range(distance),
    ylim = c(0, 1), yaxp = c(0, 1, 4),
    sub = bquote("D(" * alpha * "," ~ beta * ")" == .(round(disp, 3)))
  )
  ## 95% confidence interval
  # graphics::polygon(
  #   x = c(rev(t$distance), t$distance),
  #   y = c(rev(t$roll_disp_CI[, 1]), t$roll_disp_CI[, 2]),
  #   col = "grey90", border = FALSE, lty = 3
  # )
  graphics::points(cdist ~ distance,
    data = t, col = t$regime
  )
  graphics::lines(roll_disp ~ distance, data = t, type = "S", col = "#85112AFF")
  # graphics::abline(h = disp, col = "black", lty = 2) # dispersion


  # rose plot
  grDevices::dev.new()
  rose(
    azi,
    weights = 1 / unc,
    sub = subtitle_rose,
    main = "Rose diagram",
    mtext = "PoR"
  )
  # rose_stats(azi, weights = 1 / unc)
  rose_line(prd, radius = 1.1, col = "#009E73") # show the predicted direction
  grDevices::palette("default")
}

#' Plot data in PoR map
#'
#' @param x,pb `sf` objects of the data points and the plate
#' boundary geometries in the geographical coordinate system
#' @param PoR Pole of Rotation. \code{"data.frame"} or object of class
#' \code{"euler.pole"}
#' containing the geographical coordinates of the Pole of Rotation
#' @param type Character. Type of plate boundary (optional). Can be
#' \code{"out"}, \code{"in"}, \code{"right"}, or
#' \code{"left"} for outward, inward, right-lateral, or left-lateral
#' moving plate boundaries, respectively. If \code{"none"} (the default), only
#' the PoR-equivalent azimuth is returned.
#'
#' @param show.deviation logical.
#' Whether the data should be color-coded according to the deviation from the
#' prediction, or according to the stress regime? Is ignored if `type=='none'`.
#' @param ... optional arguments passed to [tectonicr.colors()]
#'
#' @returns plot
#'
#' @importFrom sf st_coordinates st_geometry
#' @importFrom dplyr arrange
#' @export
#'
#' @seealso [PoR_shmax()], [axes()], [tectonicr.colors()]
#'
#' @examples
#' data("nuvel1")
#' na_pa <- subset(nuvel1, nuvel1$plate.rot == "na")
#'
#' data("plates")
#' plate_boundary <- subset(plates, plates$pair == "na-pa")
#'
#' data("san_andreas")
#' PoR_map(san_andreas, PoR = na_pa, pb = plate_boundary, type = "right", show.deviation = TRUE)
PoR_map <- function(x, PoR, pb = NULL, type = c("none", "in", "out", "right", "left"),
                    show.deviation = FALSE, ...) {
  val <- val2 <- character()
  type <- match.arg(type)
  x_por_df <- PoR_shmax(x, PoR, type = type)
  if (type == "none") {
    x_por_df <- data.frame(azi.PoR = x_por_df)
  }

  x_por_coords <- geographical_to_PoR_sf(x, PoR) |>
    sf::st_coordinates()
  por_crs <- PoR_crs(PoR)

  pb_por <- geographical_to_PoR_sf(pb, PoR)

  if (show.deviation & type != "none") {
    cols <- tectonicr.colors(abs(x_por_df$cdist), categorical = FALSE, ...)
    legend.title <- "Circular distance"
  } else {
    cols <- tectonicr.colors(
      x$regime,
      pal = stress_colors(),
      categorical = TRUE,
      ...
    )
    legend.title <- "Stress regime"
  }

  col.legend <- data.frame(col = cols, val = names(cols)) |>
    mutate(val2 = gsub("\\(", "", val), val2 = gsub("\\[", "", val2)) |>
    unique() |>
    dplyr::arrange(val2)

  plot(
    x_por_coords[, 1],
    x_por_coords[, 2],
    cex = 0,
    xlab = expression("PoR longitude (" * degree * ")"),
    ylab = expression("PoR latitude (" * degree * ")"),
    asp = 1
  )
  graphics::abline(
    h = seq(-90, 90, 5),
    v = seq(-180, 180, 5),
    col = "grey",
    lty = 2
  )
  axes(
    x_por_coords[, 1],
    x_por_coords[, 2],
    x_por_df$azi.PoR,
    col = cols,
    add = TRUE
  )
  plot(sf::st_geometry(pb_por), add = TRUE)
  graphics::legend("bottomleft",
    inset = .05, cex = .75,
    legend = col.legend$val, title = legend.title, fill = col.legend$col,
    bty = "o", bg = "white"
  )
}
