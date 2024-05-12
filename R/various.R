#' Numerical values to World Stress Map Quality Ranking
#'
#' Assigns numeric values of the precision of each measurement to the
#' categorical quality ranking of the World Stress Map (A, B, C, D).
#'
#' @param x Either a string or a character vector of WSM quality ranking
#'
#' @returns \code{"integer"} or vector of type \code{"integer"}
#'
#' @references Heidbach, O., Barth, A., M<U+00FC>ller, B., Reinecker, J.,
#' Stephansson, O., Tingay, M., Zang, A. (2016). WSM quality
#' ranking scheme, database description and analysis guidelines for stress
#' indicator. *World Stress Map Technical Report* **16-01**, GFZ German Research
#' Centre for Geosciences. \doi{10.2312/wsm.2016.001}
#'
#' @name parse_wsm
#' @examples
#' parse_wsm_quality(c("A", "B", "C", "D", NA))
#' data("san_andreas")
#' parse_wsm_quality(san_andreas$quality)
parse_wsm_quality <- function(x) {
  c(
    "A" = 15,
    "B" = 20,
    "C" = 25,
    "D" = 40
  )[x]
}

#' @rdname parse_wsm
#' @export
quantise_wsm_quality <- function(x) {
  .Deprecated(parse_wsm_quality)
  as.numeric(sapply(X = x, FUN = parse_wsm_quality))
}


#' Quick analysis of a stress data set
#'
#' Returns the converted azimuths, distances to the plate boundary,
#' statistics of the model, and some plots.
#'
#' @param x \code{data.frame} or `sf` object containing the coordinates of the point(s)
#' (\code{lat}, \code{lon}), the direction of
#' \eqn{\sigma_{Hmax}}{SHmax} \code{azi} and its standard deviation
#' \code{unc} (optional)
#' @param PoR Pole of Rotation. \code{data.frame} or object of class \code{"euler.pole"}
#' containing the geographical coordinates of the Euler pole
#' @param type Character. Type of plate boundary (optional). Can be
#' \code{"out"}, \code{"in"}, \code{"right"}, or
#' \code{"left"} for outward, inward, right-lateral, or left-lateral
#' moving plate boundaries, respectively. If \code{"none"} (the default), only
#' the PoR-equivalent azimuth is returned.
#' @param pb (optional) `sf` object of the plate boundary geometries in the geographical
#' coordinate system
#' @param plot (logical). Whether to produce a plot additional to output.
#' @param ... optional arguments to [distance_from_pb()]
#'
#' @returns list containing the following values:
#' \describe{
#' \item{`results`}{data.frame showing the the coordinate and azimuth conversions
#' (`lat.PoR`, `lon.PoR`, and `azi.PoR`), the predicted azimuths (`prd`),
#' deviation angle from predicted (`dev`), circular distance (`cdist`),
#' misfit to predicted stress direction (`nchisq`) and, if given, distance to tested
#' plate boundary (`distance`)}
#' \item{`stats`}{array with circular (weighted) mean, circular standard deviation, circular variance, circular dispersion, the 95% confidence angle, and the normalized Chi-squared test statistic}
#' \item{`test`}{list containting the test results of the (weighted) Rayleigh test against the uniform distribution about the predicted  orientation.}
#' }
#'
#' @export
#'
#' @seealso [PoR_shmax()], [distance_from_pb()], [norm_chisq()], [quick_plot()]
#'
#' @examples
#' \donttest{
#' data("nuvel1")
#' na_pa <- subset(nuvel1, nuvel1$plate.rot == "na")
#'
#' data("plates")
#' plate_boundary <- subset(plates, plates$pair == "na-pa")
#'
#' data("san_andreas")
#' stress_analysis(san_andreas, na_pa, type = "right", plate_boundary, plot = TRUE)
#' }
stress_analysis <- function(x, PoR, type = c("none", "in", "out", "right", "left"), pb, plot = TRUE, ...) {
  type <- match.arg(type)
  stopifnot(is.logical(plot))
  tangential <- ifelse(type %in% c("right", "left"), TRUE, FALSE)
  res <- PoR_shmax(x, PoR, type)
  res <- cbind(res, PoR_coordinates(x, PoR))
  if (!missing(pb)) {
    res$distance <- distance_from_pb(x, PoR, pb, tangential, ...)
  }
  prd <- res$prd

  mean <- circular_mean(res$azi.PoR, 1 / x$unc)
  sd <- circular_sd(res$azi.PoR, 1 / x$unc)
  var <- circular_var(res$azi.PoR, 1 / x$unc)
  disp <- circular_dispersion(res$azi.PoR, prd, 1 / x$unc)
  conf <- confidence_angle(res$azi.PoR, w = 1 / x$unc)
  nchisq <- norm_chisq(res$azi.PoR, prd, unc = x$unc)
  rayleigh <- weighted_rayleigh(res$azi.PoR, prd, unc = x$unc)

  if (plot) {
    PoR_map(x, PoR, pb, type = type, deviation = TRUE)
    grDevices::dev.new()
    quick_plot(azi = res$azi.PoR, distance = res$distance, unc = x$unc, regime = x$regime, prd = prd)
  }

  list(
    result = res,
    stats =
      rbind(mean = mean, sd = sd, var = var, dispersion = disp, conf95 = conf, norm_chisq = nchisq),
    test = rayleigh
  )
}

#' Extract azimuths of line segments
#'
#' @param x sf object of type `"LINESTRING"` or `"MULTILINESTRING"`
#'
#' @return sf object of type `"POINT"` with the columns and entries of the first row of `x`
#'
#' @details
#' It is recommended to perform `line_azimuth()` on single line objects, i.e.
#' type `"LINESTRING"`, instead of `"MULTILINESTRING"`. This is because the azimuth
#' of the last point of a line will be calculated to the first point of the
#' next line otherwise. This will cause a warning message. For `MULTILINESTRING`
#' objects, use `lines_azimuths()`.
#'
#' @importFrom sf st_cast st_coordinates st_as_sf st_crs st_drop_geometry
#' @export
#'
#' @name line_azimuth
#' @examples
#' data("plates")
#' subset(plates, pair == "af-eu") |>
#'   smoothr::densify() |>
#'   line_azimuth()
#'
#' \dontrun{
#' lines_azimuths(plates)
#' }
NULL

#' @rdname line_azimuth
#' @export
line_azimuth <- function(x) {
  if (nrow(x) > 1) warning("It is recommended to only use single line objects")
  if (any(sf::st_geometry_type(x) == "MULTILINESTRING")) warning("It is recommended to only use single line objects")
  mat <- x |>
    sf::st_cast("POINT") |>
    sf::st_coordinates()

  n <- nrow(mat)

  a <- numeric()
  for (i in 1:(n - 1)) {
    a[i] <- get_azimuth(mat[i, 2], mat[i, 1], mat[i + 1, 2], mat[i + 1, 1])
  }
  data.frame(
    x = mat[1:(n - 1), 1],
    y = mat[1:(n - 1), 2],
    azi = a
  ) |>
    cbind(sf::st_drop_geometry(x[1, ])) |>
    sf::st_as_sf(coords = c("x", "y"), crs = sf::st_crs(x))
}

#' @rdname line_azimuth
#' @export
lines_azimuths <- function(x) {
  for (i in 1:nrow(x)) {
    ai <- line_azimuth(x[i, ])
    if (i == 1) {
      a <- ai
    } else {
      a <- rbind(a, ai)
    }
  }
  return(a)
}
