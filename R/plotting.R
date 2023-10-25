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
    dat <- data.frame(val = x) # |> mutate(val = ifelse(is.na(val), "NA", val))
    if (is.null(pal)) {
      val <- unique(dat$val)
      n <- length(val)
      if (n < 9) {
        pal <- RColorBrewer::brewer.pal(n, name = "Set2")
      } else {
        pal <- RColorBrewer::brewer.pal(n, name = "Set3")
      }
      names(pal) <- val
    }
    # preliminary solution for NA values in data.set
    cols <- data.frame(code = pal, val = names(pal), row.names = NULL)
    xpal <- dplyr::left_join(dat, cols) |>
      dplyr::mutate(code = ifelse(is.na(code), na.value, code))
    ret <- xpal$code
    names(ret) <- dat$val
  } else {
    breaks <- pretty(x, n = n + 1)
    n2 <- length(breaks) - 1
    # order = findInterval(x, sort(x))
    if (is.null(pal)) {
      cols <- viridis::viridis(n = n2, ...) # [order]
    } else {
      cols <- pal(n = n2, categorical = categorical, ...) # [order]
    }

    ret <- cut(x, breaks = breaks, labels = cols, include.lowest = TRUE)
    ret <- as.character(ret)
    names(ret) <- cut(x, breaks = breaks, include.lowest = TRUE)
  }
  ret
}


#' Color palette for stress regime
#'
#' @return function
#' @export
#'
#' @examples
#' stress_colors()
stress_colors <- function() {
  sc <- c("#D55E00", "#E69F00", "#009E73", "#56B4E9", "#0072B2", "grey60")
  names(sc) <- c("N", "NS", "S", "TS", "T", "U")
  sc
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

#' Rose Diagram
#'
#' @param clockwise Logical. Whether angles increase in the
#' clockwise direction (`clockwise=TRUE`, the default) or anti-clockwise,
#' counter-clockwise direction (`FALSE`).
#' @param main Character string specifying the title of the plot.
#' @param at Optional vector of angles at which tick marks should be plotted.
#' Set `at=numeric(0)` to suppress tick marks.
#' @param start The starting direction for measurement of angles, that is, the
#' spatial direction which corresponds to a measured angle of zero. Either a
#' character string giving a compass direction (`"N"` for north, `"S"` for south,
#' `"E"` for east, or `"W"` for west) or a number giving the angle from the the
#' horizontal (East) axis to the starting direction. For example, if
#' `clockwise=FALSE`, then `start=90` and `start="N"` are
#' equivalent. The default is to measure angles anti-clockwise from the
#' horizontal axis (East direction).
#' @param labels Either a logical value indicating whether to plot labels
#' next to the tick marks, or a vector of labels for the tick marks.
#' @param cborder logical. Border of rose plot.
#' @param ... optional arguments passed to `plot.default()`
#' @importFrom spatstat.geom disc
#' @importFrom spatstat.utils dont.complain.about do.call.matched resolve.defaults graphicsPars
#' @importFrom spatstat.explore circticks
#'
#' @return none
#'
#' @keywords internal
rose_baseplot <- function(start = -90,
                          clockwise = TRUE, main = NULL, labels = TRUE,
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
    at = at, unit = unit, start = start, clockwise = clockwise, labels = labels
  )
  return(invisible(result))
}


rose_histogram <- function(x, ..., start = -90,
                           clockwise = TRUE, main, labels = TRUE, at = NULL,
                           cborder = TRUE, axial = FALSE) {
  if (missing(main) || is.null(main)) {
    main <- spatstat.utils::short.deparse(substitute(x))
  }

  bks <- x$mids
  bw <- x$binwidth
  y <- x$density

  rose_baseplot(start, clockwise, main, labels, at, cborder)

  for (i in seq_along(y)) {
    rose_fan(bks[i], d = bw, radius = y[i], axial = axial, add = TRUE, ...)
  }

  # ang <- spatstat.explore::ang2rad(
  #     bks, unit = unit, start = start, clockwise = clockwise
  #     )
  #   ang <- deg2rad(bks)
  #   eps <- min(diff(ang), pi / 128) / 2
  #   eps <- min(bw, pi / 128) / 2
  #   for (i in seq_along(y)) {
  #     aa <- seq(ang[i], ang[i + 1], by = eps)
  #     aa[length(aa)] <- ang[i + 1]
  #     yi <- y[i]
  #     xx <- c(0, yi * cos(aa), 0)
  #     yy <- c(0, yi * sin(aa), 0)
  #     spatstat.utils::do.call.matched(polygon, list(x = xx, y = yy, ...))
  #   }
  # }
  # return(invisible(result))
}

graphicsAargh <- c(
  "density", "angle", "col", "border",
  "xlim", "ylim", "xlab", "ylab", "axes"
)

rose_freq <- function(x, bins = NULL, ..., weights = NULL, binwidth = NULL,
                      round_binwidth = 0, equal_area = TRUE,
                      main = NULL, axial = TRUE) {
  if (missing(main) || is.null(main)) {
    main <- short.deparse(substitute(x))
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
  otherargs <- h$otherargs

  freqs <- spatstat.geom::whist(
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


symmetric_bw <- function(x) {
  allowed <- c(2, 4, 6, 8, 10, 12, 18, 20, 24, 30, 36, 40, 60, 72, 90, 120, 180)
  target.index <- which(abs(allowed - x) == min(abs(allowed - x)))
  allowed[target.index] |> min()
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
#' @param clockwise Logical. Whether angles increase in the
#' clockwise direction (`clockwise=TRUE`, the default) or anti-clockwise,
#' counter-clockwise direction (`FALSE`).
#' @param muci logical. Whether the mean and its 95% CI are added to the plot
#' or not.
#' @param round_binwidth integer. Number of decimal places of bin width (0 by
#' default).
#' @param mtext character. String to be drawn at the top margin of the plot
#' (`"N"` by default)
#' @param main,sub Character string specifying the title and subtitle of the
#' plot. If `sub = NULL`, it will show the bin width.
#' @param at Optional vector of angles at which tick marks should be plotted.
#' Set `at=numeric(0)` to suppress tick marks.
#' @param col fill color of bins
#' @param dots logical. Whether a circular dot plot should be added
#' (`FALSE` is the default).
#' @param dot_cex,dot_pch,dot_col Plotting arguments for circular dot plot
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
#' @return A window (class `"owin"`) containing the plotted region.
#'
#' @importFrom spatstat.explore rose
#' @importFrom graphics hist title points
#' @importFrom stats na.omit
#'
#' @export
#'
#' @examples
#' x <- rvm(100, mean = 90, k = 5)
#' rose(x, axial = FALSE, border = TRUE)
#'
#' data("san_andreas")
#' rose(san_andreas$azi, dots = TRUE, main = "dot plot")
#' rose(san_andreas$azi, weights = 1 / san_andreas$unc, main = "weighted")
rose <- function(x, weights = NULL, binwidth = NULL, bins = NULL, axial = TRUE,
                 equal_area = TRUE, clockwise = TRUE, muci = TRUE,
                 round_binwidth = 0, mtext = "N", main = NULL, sub = NULL,
                 at = seq(0, 360 - 45, 45),
                 col = "grey", dots = FALSE, dot_pch = 1, dot_cex = 1,
                 dot_col = "grey", ...) {
  if (missing(main) || is.null(main)) {
    main <- short.deparse(substitute(x))
  }

  freqs <- rose_freq(
    x,
    bins = bins, ..., weights = weights, binwidth = binwidth,
    round_binwidth = round_binwidth, equal_area = equal_area,
    main = main, axial = axial
  )

  rose_histogram(freqs, ...,
    clockwise = clockwise,
    col = col, axial = axial,
    main = main, labels = TRUE, at = at, cborder = TRUE
  )

  if (dots) {
    rose_dots(x, axial, cex = dot_cex, pch = dot_pch, col = dot_col)
  }

  if (is.null(sub)) sub <- paste("Bin width:", freqs$binwidth)
  graphics::title(main = NULL, sub = sub, ylab = NULL)
  graphics::mtext(mtext)

  if (muci) rose_stats(x, weights = weights, axial = axial)
}

rose_dots <- function(x, axial, ...) {
  if (axial) {
    x_shift <- (x + 180) %% 360
    x <- c(x, x_shift)
  }

  scale <- 1.1 #* max(freqs$density)
  u <- deg2rad(90 - x)
  n <- length(x)
  z <- cos(u) * scale
  y <- sin(u) * scale
  # if (stack == FALSE) {
  graphics::points(z, y, ...)
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

#' Lines and fans in rose diagram
#'
#' @param x angles in degrees
#' @param d width of a fan (in degrees)
#' @param radius of the rose diagram
#' @param axial Logical. Whether x are uniaxial (`axial=FALSE`)
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
#' a <- c(0, 10, 45)
#' radius <- c(.7, 1, .2)
#' lwd <- c(2, 1, .75)
#' col <- c(1, 2, 3)
#' rose_line(c(0, 10, 45), radius = radius, axial = FALSE, add = FALSE, lwd = lwd, col = col)
NULL

#' @rdname rose_geom
#' @export
rose_line <- function(x, radius = 1, axial = TRUE, add = TRUE, ...) {
  xrad <- deg2rad(90 - x)
  tx <- radius * cos(xrad)
  ty <- radius * sin(xrad)

  if (!add) {
    rose_baseplot()
  }
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

  if (!add) {
    rose_baseplot()
  }

  graphics::polygon(x = xx, y = yy, ...)
  if (axial) {
    graphics::polygon(x = -xx, y = -yy, ...)
  }
  invisible()
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
#' (`"mean"`, the default) or the circular Quasi Median (`"median"`)
#' @param spread character. The measure of spread to be plotted as a fan.
#' Either the 95% confidence interval (`"CI"`, the default), the circular
#' standard deviation (`"sd"`), or the Quasi interquartile range on the circle
#' (`"IQR"`). `NULL` if no fan should be drawn.
#' @param avg.col color for the average line
#' @param avg.lty line type of the average line
#' @param avg.lwd  line width of the average line
#' @param spread.col color of the spread fan
#' @param spread.border logical. Whether to draw a border of the fan or not.
#' @param spread.lty line type of the spread fan's border
#' @param spread.lwd line width of the spread fan's border
#' @param add logical.
#' @param ... optional arguments to `rose_baseplot()` if add is `FALSE`.
#' @importFrom ggplot2 alpha
#'
#' @seealso [rose()] for plotting the rose diagram, and
#' [circular_mean()], [circular_median()], [confidence_interval()],
#' [circular_sd()], [circular_IQR()] for statistical parameters.
#'
#' @returns No return value, called for side effects
#' @export
#'
#' @examples
#' data("san_andreas")
#' rose(san_andreas$azi, weights = 1 / san_andreas$unc, muci = FALSE)
#' rose_stats(san_andreas$azi, weights = 1 / san_andreas$unc, avg = "median", spread = "IQR")
rose_stats <- function(x, weights = NULL, axial = TRUE, avg = c("mean", "median"), spread = c("CI", "sd", "IQR"),
                       avg.col = "#85112AFF", avg.lty = 2, avg.lwd = 1.5,
                       spread.col = ggplot2::alpha("#85112AFF", .2), spread.border = FALSE, spread.lty = NULL, spread.lwd = NULL, add = TRUE, ...) {
  avg <- match.arg(avg)
  mu <- switch(avg,
    mean = circular_mean(x, weights, axial),
    median = circular_median(x, weights, axial)
  )
  mu_text <- switch(avg,
    mean = "Mean: ",
    median = "Median: "
  )

  if (!is.null(spread)) {
    spread <- match.arg(spread)
    ci <- switch(spread,
      CI = confidence_angle(x, w = weights, axial = axial),
      sd = circular_sd(x, weights, axial),
      IQR = circular_IQR(x, weights, axial)
    )
    rose_fan(mu, 2 * ci,
      radius = 1.1, axial = axial, col = spread.col,
      border = spread.border, lty = spread.lty, lwd = spread.lwd,
      add = add, ...
    )
  }

  rose_line(mu,
    radius = 1.1, axial = axial, col = avg.col, lty = avg.lty,
    lwd = avg.lwd, add = add, ...
  )
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
#' [circular_dispersion()], [confidence_angle()], [norm_chisq()],
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
#' quick_plot(res$azi.PoR, (d), res$prd, san_andreas$unc, san_andreas$regime)
quick_plot <- function(azi, distance, prd, unc = NULL, regime, width = 51) {
  if (missing(regime)) {
    regime <- rep(NA, length(azi))
  }
  nchisq_i <- numeric()
  regime <- ifelse(is.na(regime), "U", regime)

  t <- data.frame(azi, distance, prd, unc,
    regime = factor(regime, levels = c("U", "N", "NS", "S", "TS", "T"))
  ) |>
    dplyr::arrange(distance) |>
    dplyr::mutate(
      nchisq_i = (deviation_norm(azi - prd) / unc)^2 / (90 / unc)^2,
      cdist = circular_distance(azi, prd),
      roll_mean = roll_circstats(azi, w = 1 / unc, FUN = circular_mean, width = width),
      # roll_conf95 = roll_confidence(azi, .95, 1 / unc, width = width),
      roll_sd = roll_circstats(azi, w = 1 / unc, FUN = circular_sd, width = width) / 2,
      roll_nchisq = roll_normchisq(azi, prd, unc, width = width),
      roll_disp = roll_dispersion(azi, prd, w = 1 / unc, width = width) # ,
      # roll_disp_CI = roll_dispersion_CI(azi, prd, w = 1 / unc, R = 100, width = width)
    )

  # add lower and upper period to data for plotting
  tmin <- dplyr::mutate(t, azi = azi - 180)
  tmax <- dplyr::mutate(t, azi = azi + 180)
  t2 <- rbind(tmin, t, tmax)

  nchisq <- norm_chisq(azi, prd, unc)
  suppressMessages(
    rt <- weighted_rayleigh(azi, prd = prd, unc = unc)
  )
  azi.PoR.mean <- circular_mean(azi, 1 / unc)
  azi.PoR.sd <- circular_sd(azi, 1 / unc)
  disp <- circular_dispersion(azi, prd, 1 / unc)
  CI <- confidence_interval(azi, w = 1 / unc)
  CI_ang <- confidence_angle(azi, w = 1 / unc)

  subtitle <-
    paste0(
      "Disp: ", round(disp, 3), " | 95% CI: ", round(CI$conf.interval[1]), "\u00B0 - ",
      round(CI$conf.interval[2]), "\u00B0 | R: ",
      signif(rt$statistic, 2), " (", signif(rt$p.value, 2), ")"
    )
  subtitle_rose <- paste0(
    "N: ", length(azi),
    "\nMean azimuth: ", round(azi.PoR.mean, 1), "\u00B0 \u00B1 ", round(CI_ang, 1),
    "\u00B0"
  )
  grDevices::palette(c("grey60", "#D55E00", "#E69F00", "#009E73", "#56B4E9", "#0072B2"))

  # distance plot
  ## create empty plot
  graphics::plot(0,
    type = "n",
    xlab = "Distance from plate boundary", ylab = "Azimuth wrt. PoR (\u00B0)",
    sub = subtitle,
    main = "Distance from plate boundary vs. azimuth",
    xlim = range(distance),
    ylim = c(0, 180), yaxp = c(0, 180, 8)
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
    xlab = "Distance from plate boundary", ylab = expression(Norm ~ chi[i]^2),
    main = "Deviation from prediction",
    xlim = range(distance),
    ylim = c(0, 1), yaxp = c(0, 1, 4),
    sub = paste0("Norm \u03C7\u00B2: ", round(nchisq, 2))
  )
  graphics::lines(roll_nchisq ~ distance, data = t, type = "S", col = "#85112AFF")
  graphics::abline(h = .15, col = "black", lty = 2)

  # Dispersion plot
  grDevices::dev.new()
  graphics::plot(0,
    type = "n",
    xlab = "Distance from plate boundary", ylab = "Circular distance",
    main = "Circular dispersion around prediction",
    xlim = range(distance),
    ylim = c(0, 1), yaxp = c(0, 1, 4),
    sub = paste0("Disp: ", round(disp, 3))
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
  rose(azi, weights = 1 / unc, sub = subtitle_rose, main = "Rose diagram")
  # rose_stats(azi, weights = 1 / unc)
  # rose_line(prd, radius = 1.1, col = "#009E73") # show the predicted direction
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
#' @param deviation logical.
#' Whether the data should be color-coded according to the deviation from the
#' prediction, or according to the stress regime?
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
#' PoR_map(san_andreas, PoR = na_pa, pb = plate_boundary, type = "right", deviation = TRUE)
PoR_map <- function(x, PoR, pb = NULL, type = c("none", "in", "out", "right", "left"),
                    deviation = FALSE, ...) {
  val <- val2 <- NULL
  x_por_df <- PoR_shmax(x, PoR, type = type)

  x_por_sf <- geographical_to_PoR_sf(x, PoR)
  x_por_coords <- sf::st_coordinates(x_por_sf)
  por_crs <- PoR_crs(PoR)

  pb_por <- geographical_to_PoR_sf(pb, PoR)

  if (deviation) {
    cols <- tectonicr.colors(abs(x_por_df$cdist), categorical = FALSE, ...)
    legend.title <- "|Circular distance|"
  } else {
    cols <- tectonicr.colors(x$regime, pal = stress_colors(), categorical = TRUE, ...)
    legend.title <- "Stress regime"
  }

  col.legend <- data.frame(col = cols, val = names(cols)) |>
    mutate(val2 = gsub("\\(", "", val), val2 = gsub("\\[", "", val2)) |>
    unique() |>
    dplyr::arrange(val2)

  plot(x_por_coords[, 1], x_por_coords[, 2],
    cex = 0,
    xlab = "PoR longitude (\u00B0)", ylab = "PoR latitude (\u00B0)", asp = 1
  )
  graphics::abline(h = seq(-90, 90, 5), v = seq(-180, 180, 5), col = "grey", lty = 2)
  axes(x_por_coords[, 1], x_por_coords[, 2], x_por_df$azi.PoR, col = cols, add = TRUE)
  plot(sf::st_geometry(pb_por), add = TRUE)
  graphics::legend("bottomleft",
    inset = .05, cex = .75,
    legend = col.legend$val, title = legend.title, fill = col.legend$col
  )
}
