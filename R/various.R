#' Numerical values to World Stress Map Quality Ranking
#'
#' Assigns numeric values of the precision of each measurement to the
#' categorical quality ranking of the World Stress Map (A, B, C, D).
#'
#' @param regime Either a string or a character vector of WSM quality ranking
#'
#' @returns \code{"integer"} or vector of type \code{"integer"}
#'
#' @references Heidbach, O., Barth, A., Müller, B., Reinecker, J.,
#' Stephansson, O., Tingay, M., Zang, A. (2016). WSM quality
#' ranking scheme, database description and analysis guidelines for stress
#' indicator. *World Stress Map Technical Report* **16-01**, GFZ German Research
#' Centre for Geosciences. \doi{10.2312/wsm.2016.001}
#'
#' @export
#'
#' @examples
#' quantise_wsm_quality(c("A", "B", "C", "D", NA))
#' data("san_andreas")
#' quantise_wsm_quality(san_andreas$quality)
quantise_wsm_quality <- function(regime) {
  as.numeric(sapply(X = regime, FUN = regime2unc))
}

regime2unc <- function(x) {
  c(
    "A" = 15,
    "B" = 20,
    "C" = 25,
    "D" = 40
  )[x]
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
