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
#' @importFrom graphics par
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
  graphics::par(xpd = TRUE)
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

  if (isFALSE(add)) circular_plot(smain = main, labels, at = at, cborder = cborder)

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

  if (isTRUE(axial)) {
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
  if (isTRUE(axial)) {
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

rose_grid <- function(angles, radii, add = TRUE) {
  rose_line(angles, col = "grey80", lty = 2, add = add)
  for (i in radii) {
    plot(spatstat.geom::disc(i), col = NA, border = "grey80", lty = 2, add = add)
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
#' @param origin.text character. String to be drawn at the top margin of the plot
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
#' @param grid logical. Whether to add a grid. Default is `FALSE`.
#' @param grid.lines,grid.circles numeric. Adds a sequence of straight grid
#' lines and circles based on angles and radii, respectively. Ignored when
#' `grid=FALSE`
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
#' @importFrom graphics hist title points text
#' @importFrom stats na.omit
#'
#' @family rose-plot
#'
#' @export
#'
#' @examples
#' x <- rvm(100, mean = 90, k = 5)
#' rose(x, axial = FALSE, border = TRUE, grid = TRUE)
#'
#' data("san_andreas")
#' rose(san_andreas$azi, main = "equal area")
#' rose(san_andreas$azi, equal_area = FALSE, main = "equal angle")
#'
#' # weighted frequencies:
#' rose(san_andreas$azi, weights = 1 / san_andreas$unc, main = "weighted")
#'
#' # add dots:
#' rose(san_andreas$azi, dots = TRUE, main = "dot plot", jitter = .2)
#'
#' # stack dots:
#' rose(san_andreas$azi,
#'   dots = TRUE, stack = TRUE, dot_cex = 0.5, dot_pch = 21,
#'   main = "stacked dot plot"
#' )
rose <- function(x, weights = NULL, binwidth = NULL, bins = NULL, axial = TRUE,
                 equal_area = TRUE, muci = TRUE,
                 round_binwidth = 0, origin.text = "N", main = NULL, sub = NULL,
                 at = seq(0, 360 - 45, 45), cborder = TRUE, labels = TRUE,
                 col = "grey", dots = FALSE, dot_pch = 1, dot_cex = 1,
                 dot_col = "slategrey", stack = FALSE, jitter_factor = 0,
                 grid = FALSE, grid.lines = seq(0, 135, 45), grid.circles = seq(.2, 1, .2),
                 add = FALSE, ...) {
  if (isFALSE(add)) {
    if (missing(main) || is.null(main)) {
      main <- spatstat.utils::short.deparse(substitute(x))
    }
    circular_plot(main = main, labels = labels, at = at, cborder = cborder)
  }

  if (grid) {
    rose_grid(angles = grid.lines, radii = grid.circles)
  }

  f <- if (isTRUE(axial)) 2 else 1
  x <- x %% (360 / f)

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
  #graphics::mtext(mtext, font = 2)

  graphics::text(0, 1 + ifelse(labels, 0.5, 0.3), origin.text, font = 2)

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
#' @family rose-plot
#' @examples
#' angles <- c(0, 10, 45)
#' radius <- c(.7, 1, .2)
#' lwd <- c(2, 1, .75)
#' col <- c(1, 2, 3)
#' rose_line(angles, radius = radius, axial = FALSE, add = FALSE, lwd = lwd, col = col)
#'
#'
NULL

#' @rdname rose_geom
#' @export
rose_line <- function(x, radius = 1, axial = TRUE, add = TRUE, ...) {
  xrad <- deg2rad(90 - x)
  tx <- radius * cos(xrad)
  ty <- radius * sin(xrad)

  if (isFALSE(add)) circular_plot()
  graphics::segments(0, 0, tx, ty, ...)
  if (isTRUE(axial)) {
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

  if (isFALSE(add)) circular_plot()

  graphics::polygon(x = xx, y = yy, ...)
  if (isTRUE(axial)) {
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
#' Either Batchelet's 95% confidence interval by (`"CI"`, the default),
#' Fisher's 95% confidence interval (`"fisher"`), the circular
#' standard deviation (`"sd"`), the Quasi interquartile range on the circle
#' (`"IQR"`), or the sample median deviation (`"mdev"`). `NULL` if no fan should be drawn.
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
#' @family rose-plot
#'
#' @seealso [circular_mean()], [circular_median()], [circular_sample_median()],
#' [confidence_interval()], [confidence_interval_fisher()],
#' [circular_sd()], [circular_IQR()], [circular_sample_median_deviation()]
#' for statistical parameters.
#'
#' @returns plot or a two-element vector containing the calculated average and
#' spread when assigned.
#' @export
#'
#' @examples
#' data("san_andreas")
#' rose(san_andreas$azi, weights = 1 / san_andreas$unc, muci = FALSE)
#' rose_stats(san_andreas$azi, weights = 1 / san_andreas$unc, avg = "sample_median", spread = "mdev")
rose_stats <- function(x, weights = NULL, axial = TRUE, avg = c("mean", "median", "sample_median"), spread = c("CI", "fisher", "sd", "IQR", "mdev"),
                       avg.col = "#B63679FF", avg.lty = 2, avg.lwd = 1.5,
                       spread.col = ggplot2::alpha("#B63679FF", .2), spread.border = FALSE, spread.lty = NULL, spread.lwd = NULL, add = TRUE, ...) {
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
    rose_fan(mu, 2 * ci,
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

# Dot Plot ---------------------------------------------------------------------


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
#' @family rose-plot
#'
#' @return A list with information on the plot
#' @export
#'
#' @examples
#' x <- rvm(100, mean = 90, k = 5)
#'
#' # plot poinit without jitter
#' plot_points(x, add = FALSE)
#'
#' # with some jitter
#' plot_points(x, jitter_factor = .2, add = FALSE)
#'
#' # stacked dots:
#' plot_points(x, stack = TRUE, binwidth = 3, add = FALSE, xpd = TRUE)
plot_points <- function(x, axial = TRUE, stack = FALSE, binwidth = 1, cex = 1, sep = 0.025, jitter_factor = 0, ..., scale = 1.1, add = TRUE,
                        main = NULL, labels = TRUE,
                        at = seq(0, 360 - 45, 45), cborder = TRUE) {
  if (isFALSE(add)) {
    if (missing(main) || is.null(main)) {
      main <- spatstat.utils::short.deparse(substitute(x))
    }
    circular_plot(main = main, labels = labels, at = at, cborder = cborder)
  }

  f <- if (isTRUE(axial)) 2 else 1

  if (isFALSE(stack)) {
    if (isTRUE(axial)) {
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
    if (isTRUE(axial)) {
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

calc_circular_density <- function(x, z, kappa, axial) {
  nx <- length(x)
  # if (kernel == "vonmises") {
  y <- sapply(z, FUN = dvm, mean = x, kappa = kappa, axial = axial, log = FALSE)

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
  f <- if (isTRUE(axial)) 2 else 1

  if (is.null(kappa)) kappa <- est.kappa(f * x)

  if (is.null(z)) {
    z <- seq(from = from, to = to, length = n)
  } else {
    if (!is.numeric(z)) {
      stop("argument 'z' must be numeric")
    }
    namez <- deparse(substitute(z))
    z.na <- is.na(z)
    if (any(z.na)) {
      if (isTRUE(na.rm)) {
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

  calc_circular_density(x, z, kappa = kappa, axial = axial)
}

circular_lines <- function(x, y, join = FALSE, nosort = FALSE, offset = 1.1, shrink = 1, axial = TRUE, ...) {
  x <- deg2rad(90 - x)
  if (isTRUE(axial)) {
    x <- c(x, x + pi)
    y <- rep(y, 2)
  }
  n <- length(x)

  if (isFALSE(nosort)) {
    ord <- order(x)
    x <- x[ord]
    y <- y[ord]

    spacings <- diff(c(x, x[1] + 2 * pi))
    pos <- which.max(spacings)

    rot <- if (pos == n) NULL else c((pos + 1):n, seq_len(pos))
  } else {
    rot <- NULL
  }

  r <- y / shrink + offset
  z <- r * cos(x)
  w <- r * sin(x)

  if (!is.null(rot)) {
    z <- z[rot]
    w <- w[rot]
  }

  if (join) {
    z <- c(z, z[1])
    w <- c(w, w[1])
  }

  graphics::lines(x = z, y = w, ...)
  invisible(list(x = z, y = w))
}

circular_polygon <- function(x, y, nosort = FALSE, offset = 1.1, shrink = 1, axial = TRUE, ...) {
  x <- deg2rad(90 - x)
  if (isTRUE(axial)) {
    x <- c(x, x + pi)
    y <- rep(y, 2)
  }
  n <- length(x)
  if (isFALSE(nosort)) {
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
  graphics::polygon(x = z, y = w, ...)
  invisible(list(x = z, y = w))
}

#' Circular Kernel Density Plot
#'
#' Plots multiples of a von Mises density distribution in a circular plot
#'
#' @param x numeric. Data to be plotted, i.e. vector containing angles (in degrees).
#' @param kappa numeric. Concentration parameter for the von Mises distribution.
#' Small kappa gives smooth density lines. Will be estimated using [est.kappa()] if not specified.
#' @param axial Logical. Whether data are uniaxial (`axial=FALSE`)
#' or biaxial (`TRUE`, the default).
#' @param n integer. the number of equally spaced points at which the density is to be estimated.
#' @param norm.density logical. Normalize the density?
#' @param scale numeric. radius of plotted circle. Default is `1.1`.
#' @param shrink numeric. parameter that controls the size of the plotted function. Default is 1.
#' @param grid logical. Whether a grid should be added.
#' @param fill logical. Whether to fill the density curve or draw just a line (the default)
#' @param ... Further graphical parameters may also be supplied as arguments.
#' @param add logical. Add to existing plot? (`TRUE` by default).
#' @inheritParams circular_plot
#'
#' @seealso [dvm()]
#' @family rose-plot
#' @return plot or calculated densities as numeric vector
#' @export
#'
#' @examples
#' # Filled density curve inside the plot
#' plot_density(san_andreas$azi,
#'   kappa = 100,
#'   fill = TRUE, col = "#51127C80", border = "#51127CFF",
#'   grid = TRUE,
#'   add = FALSE
#' )
#'
#' # Superimpose a density curve on a rose diagram:
#' rose(san_andreas$azi, grid = TRUE)
#' plot_density(san_andreas$azi,
#'   kappa = 100, col = "#51127CFF",
#'   add = TRUE, lwd = 3
#' )
#'
#' # Corona plot: Density curve outside of a rose diagram plot:
#' rose(san_andreas$azi, dots = TRUE, stack = TRUE, dot_cex = 0.5, dot_pch = 21)
#' plot_density(san_andreas$azi,
#'   kappa = 100,
#'   scale = 1.1, shrink = 3, xpd = NA,
#'   col = "#51127CFF"
#' )
plot_density <- function(x, kappa = NULL, axial = TRUE, n = 512L, norm.density = TRUE,
                         ...,
                         fill = FALSE,
                         scale = 0, shrink = 1,
                         add = TRUE, main = NULL, labels = TRUE,
                         at = seq(0, 360 - 45, 45), cborder = TRUE, grid = FALSE) {
  if (isFALSE(add)) {
    if (missing(main) || is.null(main)) {
      main <- spatstat.utils::short.deparse(substitute(x))
    }
    circular_plot(main = main, labels = labels, at = at, cborder = cborder)
  }

  if (grid) {
    rose_grid(seq(0, 135, 45), seq(.2, 1, .2))
  }

  f <- 1
  d <- circular_density(x, kappa = kappa, n = n, axial = axial)
  if (norm.density) {
    d <- d / max(d)
    # shrink = 1/scale
  }

  if(isFALSE(fill)){
    circular_lines(seq(0, 360, length = f * n), rep(d, f), axial = FALSE, n, offset = scale, shrink = shrink, ...)
  } else {
    circular_polygon(seq(0, 360, length = f * n), rep(d, f), axial = FALSE, offset = scale, shrink = shrink, ...)
  }
  invisible(d)
}
