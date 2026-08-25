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
#' Show direction axes in a map
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
  if (isFALSE(add)) plot(x, y, cex = 0, ...)
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
#' @source https://stackoverflow.com/questions/55474143/how-to-center-geom-spoke-around-their-origin/
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


#' Azimuth visualization
#'
#' @description
#' `geom_azimuth()` visualizes axial-directional vector fields using a geom to
#' produce a new graphical layer, which allows aesthetic options.
#' This layer can be overlaid on a map to improve visualisation of mapped data.
#' The geom draws line segments (spokes) centered at (x, y) with a given
#' orientation (`angle` in degrees) and length (`radius`). By default the spoke
#' is centered using [`PositionCenterSpoke`], so that the given coordinates mark
#' the middle of the line. The azimuths are given as angles in degrees increasing clockwise from North.
#'
#' @param mapping Set of aesthetic mappings created by [ggplot2::aes()].
#' @param data A data frame. If `NULL`, the default, the data is inherited from
#'   the plot data as specified in the call to [ggplot2::ggplot()].
#' @param stat The statistical transformation to use on the data. Defaults to
#'   `"identity"`.
#' @param center Logical; if `TRUE` (the default) spokes are centered on (x, y) using
#'   [`PositionCenterSpoke`]  - useful for axial data. If `FALSE`, behaves like
#'   [ggplot2::geom_spoke()] (line starts at (x, y)) - useful for directional data
#'   (especially when in combination with `arrow()`).
#' @param na.rm If `FALSE`, the default, missing values are removed with a
#'   warning. If `TRUE`, missing values are silently removed.
#' @param show.legend Logical. Should this layer be included in the legends?
#' @param inherit.aes If `FALSE`, overrides the default aesthetics, rather than
#'   combining with them.
#' @param ... Other arguments passed on to [ggplot2::layer()]. These are often
#'   aesthetics (e.g. `colour`, `linetype`, `linewidth`, `alpha`).
#' @param radius Length of spoke
#'
#' @section Aesthetics:
#' `geom_azimuth()` understands the following aesthetics (required aesthetics in **bold**):
#' \itemize{
#'   \item **x**
#'   \item **y**
#'   \item angle (in degrees, transformed internally)
#'   \item radius
#'   \item colour
#'   \item alpha
#'   \item linewidth
#'   \item linetype
#' }
#'
#' @return A ggplot2 layer that adds axis-like spokes.
#' @seealso [ggplot2::geom_spoke()], [geom_azimuthpoint()]
#' @examples
#' set.seed(20250411)
#' df <- data.frame(
#'   x = runif(5), y = runif(5),
#'   angle_deg = rvm(5, mean = 90, kappa = 10),
#'   radius = runif(5, 0.1, 2)
#' )
#'
#' if (require("ggplot2")) {
#'   ggplot(df, aes(x, y)) +
#'     geom_azimuth(aes(angle = angle_deg), radius = .1, linewidth = 1.2, colour = "blue")
#'   if (require("grid")) {
#'     ggplot(df, aes(x, y, radius = radius)) +
#'       geom_azimuth(aes(angle = angle_deg), center = FALSE, colour = "red", arrow = grid::arrow())
#'   }
#' }
#' @export
geom_azimuth <- function(mapping = NULL, data = NULL,
                         stat = "azimuth", center = TRUE,
                         radius = NULL,
                         na.rm = FALSE, show.legend = NA,
                         inherit.aes = TRUE, ...) {
  position <- if (isTRUE(center)) "center_spoke" else "identity"
  radius_fact <- if (isTRUE(center)) .5 else 1

  ggplot2::layer(
    geom = ggplot2::GeomSpoke,
    mapping = mapping,
    data = data,
    stat = StatAzimuth,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      radius_fact = radius_fact,
      default_radius = radius,
      na.rm = na.rm,
      ...
    )
  )
}

# A simple internal stat that ensures 'radius' is always numeric
StatAzimuth <- ggplot2::ggproto("StatAzimuth", ggplot2::Stat,
  required_aes = c("x", "y", "angle"),
  compute_panel = function(data, scales, radius_fact = 1, default_radius = NULL) {
    data$angle <- (90 - as.numeric(data$angle)) * pi / 180

    # If radius is mapped, use it
    if ("radius" %in% names(data)) {
      data$radius <- as.numeric(data$radius) * radius_fact
    } else {
      # Otherwise use fixed radius
      fixed <- if (is.null(default_radius)) 1 else default_radius
      data$radius <- rep(fixed * radius_fact, nrow(data))
    }
    data
  }
)


#' Azimuth + point visualization
#'
#' @description
#' `geom_azimuthpoint()` draws line segments (spokes) like [geom_azimuth()], but also
#' places a point (marker) at the spoke's center `(x, y)`.
#'
#' Aesthetic rules:
#' - `linewidth`, `linetype` affect the spoke only
#' - `shape` affects the point only
#' - `colour`, `alpha` affect both spoke and point
#' - `size` sets the size of the point only
#'
#' @param mapping Set of aesthetic mappings created by [ggplot2::aes()].
#' @param data A data frame. If `NULL`, the default, the data is inherited from
#'   the plot data as specified in the call to [ggplot2::ggplot()].
#' @param stat The statistical transformation to use on the data. Defaults to
#'   `"identity"`.
#' @param center Logical; if `TRUE` spokes are centered on (x, y) using
#'   [`PositionCenterSpoke`]. If `FALSE`, behaves like
#'   [ggplot2::geom_spoke()] (line starts at (x, y)).
#' @param na.rm If `FALSE`, the default, missing values are removed with a
#'   warning. If `TRUE`, missing values are silently removed.
#' @param show.legend Logical. Should this layer be included in the legends?
#' @param inherit.aes If `FALSE`, overrides the default aesthetics, rather than
#'   combining with them.
#' @param size Size of the point marker (default = 2).
#' @param ... Other arguments passed on to [geom_azimuth()] and
#'   [ggplot2::geom_point()]. These may include `arrow`, `fill`, etc.
#'
#' @section Aesthetics:
#' `geom_azimuthpoint()` understands the following aesthetics (required aesthetics in **bold**):
#' \itemize{
#'   \item **x**
#'   \item **y**
#'   \item angle (in degrees, transformed internally; spoke only)
#'   \item radius (spoke only)
#'   \item colour (shared)
#'   \item alpha (shared)
#'   \item linewidth (spoke only)
#'   \item linetype (spoke only)
#'   \item shape (point only)
#'   \item size (point only, or via argument)
#'   \item fill (point only, for shapes that accept fill)
#' }
#'
#' @return A list of ggplot2 layers (spokes + points).
#' @seealso [geom_azimuth()], [ggplot2::geom_spoke()], [ggplot2::geom_point()]
#' @examples
#' set.seed(20250411)
#' df <- data.frame(
#'   x = runif(5), y = runif(5),
#'   angle_deg = rvm(5, mean = 90, kappa = 10),
#'   radius = runif(5, 0.1, 1),
#'   group = rep(1:2, length.out = 5)
#' )
#'
#' if (require("ggplot2")) {
#'   ggplot(data = df, aes(x, y, angle = angle_deg, radius = radius)) +
#'     geom_azimuthpoint(aes(colour = factor(group), shape = factor(group)),
#'       linewidth = 1.1, linetype = "dashed",
#'       size = 3, alpha = 0.8
#'     )
#' }
#'
#' @export
geom_azimuthpoint <- function(mapping = NULL, data = NULL,
                              stat = "identity", center = TRUE,
                              na.rm = FALSE, show.legend = NA,
                              inherit.aes = TRUE, size = 2,
                              ...) {
  # ----- Spoke layer mapping -----
  line_mapping <- mapping
  if (!is.null(line_mapping)) {
    # drop point-specific aesthetics from line mapping
    line_mapping <- line_mapping[!names(line_mapping) %in% c("shape")]
  }

  line_layer <- geom_azimuth(
    mapping = line_mapping, data = data,
    stat = stat, center = center,
    na.rm = na.rm, show.legend = show.legend,
    inherit.aes = inherit.aes, ...
  )

  # ----- Point layer mapping -----
  point_mapping <- mapping
  if (!is.null(point_mapping)) {
    # drop spoke-specific aesthetics from point mapping
    point_mapping <- point_mapping[!names(point_mapping) %in% c("angle", "radius", "linewidth", "linetype")]
  }

  suppressWarnings(
    point_layer <- ggplot2::geom_point(
      mapping = point_mapping, data = data,
      size = size,
      inherit.aes = inherit.aes, ...
    )
  )

  list(line_layer, point_layer)
}




#' Quantile-Quantile Linearised Plot for Circular Distributions
#'
#'
#' Uniformly distributed orientations should yield a straight line through the
#' origin. Systematic departures from linearity will indicate preferred
#' orientation.
#'
#' @param x numeric. Angles in degrees
#' @param axial Logical. Whether data are uniaxial (`axial=FALSE`)
#' @param xlab,ylab,main plot labels.
#' @param col color for the dots.
#' @param add_line logical. Whether to connect the points by straight lines?
#' @param ... graphical parameters
#'
#' @return plot
#' @importFrom graphics plot abline lines points
#' @export
#'
#' @family circ-qqplot
#' @seealso [cunif()]
#'
#' @references Borradaile, G. J. (2003). Statistics of earth
#' science data: their distribution in time, space, and orientation (Vol. 351,
#' p. 329). Berlin: Springer.
#'
#' @examples
#' set.seed(20250411)
#'
#' # von Mises distribution
#' x_vm <- rvm(100, mean = 0, kappa = 2)
#' circular_qqplot(x_vm, pch = 20)
#'
#' x_wcauchy <- rwcauchy(100, mean = 0, rho = 0.5)
#' circular_qqplot(x_wcauchy, pch = 20)
#'
#' # circular uniform data
#' x_cunif <- rcunif(100)
#' circular_qqplot(x_cunif, pch = 20)
circular_qqplot <- function(x, axial = TRUE,
                            xlab = paste("i/(n+1)"),
                            ylab = NULL, main = "Circular Quantile-Quantile Plot",
                            add_line = TRUE,
                            col = "#B63679FF", ...) {
  if (isTRUE(axial)) {
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

  graphics::plot(1,
    type = "n", xlim = c(0, 1), ylim = c(0, 1), asp = 1,
    xlab = xlab, ylab = ylab
  )
  graphics::abline(a = 0, b = 1, col = "slategrey")
  if (isTRUE(add_line)) {
    graphics::lines(xin, x)
  }
  graphics::points(xin, x, col = col, ...)
  graphics::title(main = main, sub = bquote("N" == .(n)))

  invisible(xin)
  # graphics::points(xin, x, col = "slategrey")
}

#' von Mises Quantile-Quantile Plot
#'
#' Produces a Q-Q plot of the data against a specified von Mises distribution
#' to graphically assess the goodness of fit of the model.
#'
#' @inheritParams circular_qqplot
#' @param w numeric. optional weightings for `x` to estimate `mean` and `kappa`.
#' @param mean numeric. Circular mean of the von Mises distribution. If `NULL`,
#' it will be estimated from `x`.
#' @param kappa numeric. Concentration parameter of the von Mises distribution.
#' If `NULL`, it will be estimated from `x`.
#'
#' @return plot
#' @importFrom stats ecdf
#' @importFrom graphics plot lines
#'
#' @family circ-qqplot
#' @seealso [vonmises()]
#'
#' @export
#'
#' @examples
#' set.seed(20250411)
#'
#' # von Mises distribution
#' x_vm <- rvm(100, mean = 0, kappa = 2)
#' vm_qqplot(x_vm, pch = 20)
#'
#' x_wcauchy <- rwcauchy(100, mean = 0, rho = 0.5)
#' vm_qqplot(x_wcauchy, pch = 20)
#'
#' # circular uniform data
#' x_cunif <- rcunif(100)
#' vm_qqplot(x_cunif, pch = 20)
vm_qqplot <- function(x, w = NULL, axial = TRUE, mean = NULL, kappa = NULL,
                      xlab = "von Mises quantile function",
                      ylab = "Empirical quantile function",
                      main = "von Mises Q-Q Plot",
                      col = "#B63679FF", add_line = TRUE, ...) {
  if (isTRUE(axial)) {
    f <- 2
  } else {
    f <- 1
  }

  n <- length(x)

  if (is.null(mean)) mean <- circular_mean(x, w = w, axial = axial)
  if (is.null(kappa)) kappa <- est.kappa(f * x, w = w)

  caption <- bquote(
    bar(alpha) == .(round(mean, 1)) * degree ~ "|" ~ widehat(kappa) == .(round(kappa, 1))
  )

  xf <- (x * f) %% 360
  xf <- sort(xf)

  edf <- stats::ecdf(xf)

  tqf <- qvm(edf(xf), mean * f, kappa, from = 0)

  graphics::plot(1, type = "n", xlab = xlab, ylab = ylab, xlim = c(0, 360 / f), ylim = c(0, 360 / f))
  graphics::abline(a = 0, b = 1, col = "slategrey")
  if (isTRUE(add_line)) graphics::lines(tqf / f, xf / f)
  graphics::points(tqf / f, xf / f, col = col, ...)
  graphics::mtext(caption, cex = .75)
  graphics::title(main = main, sub = bquote("N" == .(n)))
  invisible(tqf)
}


#' Plotting Stress Analysis Results
#'
#' Creates a set of plots including
#' the azimuth as a function of the distance to the plate boundary,
#' the Norm Chi-squared as a function of the distance to the plate boundary,
#' the circular distance (and dispersion) a function of the distance to the
#' plate boundary, a von-Mises QQ plot, and a rose diagram of the
#' quality-weighted frequency distribution of the azimuths.
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
#' [weighted_rayleigh()], [vm_qqplot()]
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
#' quick_plot(res$azi.PoR,
#'   distance = d, prd = res$prd, unc = san_andreas$unc,
#'   regime = san_andreas$regime
#' )
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
  graphics::lines(roll_mean ~ distance, data = t, type = "S", col = "#B63679FF")

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
    col = "#B63679FF"
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
  graphics::lines(roll_disp ~ distance, data = t, type = "S", col = "#B63679FF")
  # graphics::abline(h = disp, col = "black", lty = 2) # dispersion

  # von Mises QQ plot
  grDevices::dev.new()
  vm_qqplot(azi, axial = FALSE, pch = 20)

  # rose plot
  grDevices::dev.new()
  rose(
    azi,
    weights = 1 / unc,
    sub = subtitle_rose,
    main = "Rose diagram",
    origin.text = "PoR"
  )
  # rose_stats(azi, weights = 1 / unc)
  rose_line(prd, radius = 1.1, col = "#FB8861FF") # show the predicted direction
  grDevices::palette("default")
}

#' Map of data in Pole of Rotation reference frame
#'
#' Transforms the spatial data and azimuths into the PoR reference frame and
#' shows them in a map
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
  if (type == "none" | is.null(type)) {
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
