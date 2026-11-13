## Kriging of SHmax axes after van den Boogaart & Schaeben (2002),
## "Kriging of regionalized directions, axes, and orientations I",
## Mathematical Geology 34(5), 479-503. Implemented in the stress2grid() style.
##
## For axes (q = 2) the paper embeds the axis u in the tensor u u^T. Up to a
## constant, this is the doubled-angle vector w = (cos 2a, sin 2a), and the
## matrix-valued kriging estimate is M = 1/2 (I + [[w1, w2], [w2, -w1]]).
## Under the isotropy-of-measurements assumption, the kriging weights are
## those of ordinary kriging with the trace variogram
##     gamma(h) = E[1 - cos(2 * (a_i - a_j))]   (sill = 1)
## Predicted axis: dominant eigenvector of M, i.e. 1/2 * atan2(w2, w1).
## Confidence: |w_hat| (eigenvalues of M are (1 +/- |w_hat|) / 2).

#' Variogram model on the doubled-angle scale (sill = 1)
#' @keywords internal
.krig_gamma <- function(h, model, nugget, range_km) {
  r <- h / range_km
  g <- switch(model,
              spherical   = ifelse(r < 1, 1.5 * r - 0.5 * r^3, 1),
              exponential = 1 - exp(-3 * r),
              gaussian    = 1 - exp(-3 * r^2)
  )
  nugget + (1 - nugget) * g
}

#' Empirical trace variogram of axial data, for choosing the model
#'
#' Bins the empirical variogram gamma(h) = mean(1 - cos(2 * (a_i - a_j))) over
#' pairs of observations. Use it to choose `nugget` and `range_km` for
#' [stress2grid_kriging()], e.g. by overlaying the fitted model.
#'
#' @param x `sf` object with `azi` (degrees).
#' @param breaks numeric. Distance bin edges in km.
#' @returns data.frame with `lag` (km), `gamma` and `n` (pairs per bin).
#' @export
variogram_axial <- function(x, breaks = seq(0, 500, 50)) {
  xy <- sf::st_coordinates(x)
  lat2r <- deg2rad(xy[, 2])
  lon2r <- deg2rad(xy[, 1])
  th2 <- deg2rad(2 * x$azi)
  n <- nrow(x)
  idx <- which(upper.tri(diag(n)), arr.ind = TRUE)
  i <- idx[, 1]
  j <- idx[, 2]
  h <- .hav_km(lat2r[i], lon2r[i], lat2r[j], lon2r[j])
  g <- 1 - cos(th2[i] - th2[j])
  bins <- cut(h, breaks, include.lowest = TRUE)
  data.frame(
    lag = as.numeric(tapply(h, bins, mean)),
    gamma = as.numeric(tapply(g, bins, mean)),
    n = as.integer(table(bins))
  )
}

#' Spatial kriging of SHmax axes after van den Boogaart & Schaeben (2002)
#'
#' Ordinary kriging of the doubled-angle embedding of SHmax axes, with a
#' distance-only trace variogram (isotropy of measurements). The
#' per-observation uncertainty, given as a standard deviation in degrees, is
#' used as measurement error on the diagonal of the kriging matrix.
#'
#' @inheritParams stress2grid
#' @param R_range numeric value or vector of search radii in km. Only
#'   observations within the radius are used.
#' @param max_data integer. The number of nearest observations used for each
#'   prediction. Default is `50`.
#' @param vgm_model Variogram model: `"spherical"` (default), `"exponential"`
#'   or `"gaussian"`. The sill is fixed at 1 (see details).
#' @param nugget numeric in \eqn{[0, 1)}. Nugget of the trace variogram.
#' @param range_km numeric. Range of the variogram model in km.
#' @param use_uncertainty logical. If `TRUE` (default), `unc` (standard
#'   deviation in degrees) is used as measurement error.
#' @param axial Must be `TRUE`. Only axial data are supported.
#'
#' @returns `sf` object containing
#' \describe{
#' \item{lon,lat}{longitude and latitude in degrees}
#' \item{azi}{kriged axis direction in degrees (0 to 180), from the dominant
#'   eigenvector of the orientation tensor}
#' \item{sd}{circular standard deviation of the axis in degrees, computed
#'   from `conf` as `rad2deg(sqrt(-2 log(conf))) / 2`}
#' \item{R}{search radius in km}
#' \item{mdr}{mean distance between grid point and data / R}
#' \item{N}{number of observations used}
#' \item{conf}{confidence in \eqn{[0, 1]}: length of the kriged doubled-angle
#'   vector, capped at 1. Equals \eqn{2 \times \text{eig}_1 - 1}}
#' \item{kvar}{kriging variance on the trace-variogram scale}
#' \item{eig1}{larger eigenvalue of the kriged orientation tensor, in \eqn{[0.5, 1]}.
#'   The smaller eigenvalue is `1 - eig1`}
#' }
#' Attribute `"vgm"` records the variogram parameters used.
#'
#' @details The method is the simplified ordinary-kriging form derived by van
#' den Boogaart & Schaeben (2002) for axes in the embedding u u^T. Their
#' isotropy-of-measurements assumption gives the trace variogram used here.
#' The matrix-valued estimate is projected back onto the axes by taking the
#' dominant eigenvector, which equals the direction of the predicted
#' doubled-angle vector. The length of that vector is the confidence measure
#' described in the paper.
#'
#' Fit `nugget` and `range_km` to the empirical variogram from
#' [variogram_axial()]. The variogram has sill 1, because the trace of the
#' covariance of unit vectors is 1 (eq. 36 of the paper).
#'
#' @references van den Boogaart, K.G. and Schaeben, H. (2002).
#' Kriging of regionalized directions, axes, and orientations I.
#' Mathematical Geology, 34(5), 479-503.
#'
#' @examples
#' \dontrun{
#' data("san_andreas")
#' v <- variogram_axial(san_andreas)
#' res <- stress2grid_kriging(san_andreas, R_range = 300,
#'                            vgm_model = "spherical",
#'                            nugget = 0.1, range_km = 200)
#' head(res)
#' }
#' @export
stress2grid_kriging <- function(x,
                               grid = NULL,
                               lon_range = NULL,
                               lat_range = NULL,
                               gridsize = 2,
                               R_range = 300,
                               max_data = 50,
                               min_data = 3L,
                               min_dist_threshold = 200,
                               vgm_model = c("spherical", "exponential", "gaussian"),
                               nugget = 0,
                               range_km = 300,
                               use_uncertainty = TRUE,
                               axial = TRUE) {
  stopifnot(
    inherits(x, "sf"),
    is.numeric(gridsize) && length(gridsize) == 1,
    is.numeric(R_range) && length(R_range) >= 1,
    is.numeric(max_data) && max_data >= min_data,
    is.numeric(nugget) && length(nugget) == 1 && nugget >= 0 && nugget < 1,
    is.numeric(range_km) && length(range_km) == 1 && range_km > 0,
    is.numeric(min_dist_threshold) && min_dist_threshold > 0,
    is.logical(use_uncertainty)
  )
  if (!isTRUE(axial)) {
    stop("stress2grid_kriging() requires axial data (`axial = TRUE`).", call. = FALSE)
  }
  vgm_model <- match.arg(vgm_model)
  min_data <- as.integer(ceiling(min_data))

  x_coords <- sf::st_coordinates(x)
  lon <- x_coords[, 1]
  lat <- x_coords[, 2]
  azi <- x$azi
  n <- length(azi)
  if (n < min_data) stop("Not enough observations (fewer than `min_data`).", call. = FALSE)

  lon2r <- deg2rad(lon)
  lat2r <- deg2rad(lat)
  th2 <- deg2rad(2 * azi)     # doubled angles (w = cos, sin of 2a)

  # measurement variance on the doubled-angle scale: an angle sd s gives a
  # doubled-angle sd of 2s, and the trace variogram of a single point is 4 s^2
  sig2 <- if (use_uncertainty && "unc" %in% colnames(x)) {
    (2 * deg2rad(x$unc))^2
  } else {
    rep(0, n)
  }
  sig2[is.na(sig2)] <- 0

  # grid, identical construction to stress2grid()
  if (is.null(grid)) {
    if (is.null(lon_range) | is.null(lat_range)) {
      lon_range <- range(lon, na.rm = TRUE)
      lat_range <- range(lat, na.rm = TRUE)
    }
    grid <- sf::st_bbox(
      c(xmin = lon_range[1], xmax = lon_range[2], ymin = lat_range[1], ymax = lat_range[2]),
      crs = sf::st_crs("WGS84")
    ) |>
      sf::st_make_grid(cellsize = gridsize, what = "centers", offset = c(lon_range[1], lat_range[1])) |>
      sf::st_as_sf()
  }
  stopifnot(inherits(grid, "sf"), any(sf::st_is(grid, "POINT")))
  G <- unname(sf::st_coordinates(grid))
  nG <- nrow(G)
  glon1r <- deg2rad(G[, 1])
  glat1r <- deg2rad(G[, 2])

  results <- vector("list", length(R_range))
  cols <- c("azi", "sd", "R", "mdr", "N", "conf", "kvar", "eig1")

  for (kr in seq_along(R_range)) {
    R_k <- R_range[kr]
    mat <- matrix(NA_real_, nrow = nG, ncol = length(cols), dimnames = list(NULL, cols))

    for (i in seq_len(nG)) {
      d <- .hav_km(glat1r[i], glon1r[i], lat2r, lon2r)
      if (min(d) > min_dist_threshold) next

      # nearest observations within the search radius, at most max_data
      ord <- order(d)
      ord <- ord[d[ord] <= R_k]
      ord <- ord[seq_len(min(length(ord), max_data))]
      k <- length(ord)
      if (k < min_data) next

      # observation-observation variogram matrix with measurement error
      Dm <- matrix(.hav_km(rep(lat2r[ord], k), rep(lon2r[ord], k),
                           rep(lat2r[ord], each = k), rep(lon2r[ord], each = k)), k, k)
      Gm <- .krig_gamma(Dm, vgm_model, nugget, range_km)
      diag(Gm) <- sig2[ord]

      # ordinary kriging system with unbiasedness constraint
      A <- rbind(cbind(Gm, 1), c(rep(1, k), 0))
      g0 <- c(.krig_gamma(d[ord], vgm_model, nugget, range_km), 1)
      sol <- tryCatch(solve(A, g0), error = function(e) NULL)
      if (is.null(sol)) next

      lam <- sol[seq_len(k)]
      W <- cbind(cos(th2[ord]), sin(th2[ord]))   # k x 2 doubled-angle vectors
      wh <- colSums(lam * W)                      # kriged doubled-angle vector (uncapped)

      # orientation tensor of the kriged vector, M = 1/2 (I + [[w1, w2], [w2, -w1]])
      M <- 0.5 * matrix(c(1 + wh[1], wh[2], wh[2], 1 - wh[1]), 2, 2)
      eg <- eigen(M, symmetric = TRUE)
      v1 <- eg$vectors[, 1]
      eig1 <- eg$values[1]                        # larger eigenvalue; the smaller is 1 - eig1
      azi_i <- (atan2(v1[2], v1[1]) * 180 / pi) %% 180

      # confidence and circular standard deviation
      conf <- min(sqrt(sum(wh^2)), 1)
      sd_i <- if (conf >= 1) {
        0
      } else if (conf <= 0) {
        NA_real_
      } else {
        rad2deg(sqrt(-2 * log(conf))) / 2
      }
      kvar <- sum(lam * g0[seq_len(k)]) + sol[k + 1]

      mat[i, ] <- c(azi_i, sd_i, R_k, mean(d[ord]) / R_k, k, conf, kvar, eig1)
    }

    res_k <- data.frame(lon = G[, 1], lat = G[, 2], mat)
    results[[kr]] <- res_k[!is.na(res_k$azi), , drop = FALSE]
  }

  out <- do.call(rbind, results)
  out$N <- as.integer(out$N)
  rownames(out) <- NULL
  out <- sf::st_as_sf(out, coords = c("lon", "lat"), crs = sf::st_crs(x), remove = FALSE)
  attr(out, "vgm") <- list(model = vgm_model, nugget = nugget, range_km = range_km,
                           sill = 1, max_data = max_data)
  out
}

