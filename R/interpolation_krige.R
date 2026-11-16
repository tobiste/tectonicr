## Kriging of SHmax axes after van den Boogaart & Schaeben (2002),
## "Kriging of regionalized directions, axes, and orientations I",
## Mathematical Geology 34(5), 479-503. Implemented in the stress2grid() style.
##
## The paper embeds an axis u in the tensor u u^T and krieges that tensor.
## With ordinary kriging weights lambda (summing to 1), the kriged tensor is
##     T_hat = sum_i lambda_i * u_i u_i^T
## which is exactly the orientation tensor computed by ortensor2d() with
## weights lambda. Its dominant eigenvector is the predicted axis and its
## eigenvalues sum to 1, so the azimuth and eigenvalues come from
## ot_eigen2d(). The angles are not doubled here.
##
## The kriging weights depend on the trace variogram of the tensors,
##     gamma(h) = E || T(x) - T(x + h) ||_F^2 = E[1 - cos(2 * (a_i - a_j))],
## which is a function of distance only (isotropy of measurements, sill = 1).
## The doubling appears only in this variogram, where it is the Frobenius
## distance between tensors.

#' Variogram model of the trace variogram
#'
#' Fit the variogram gamma \eqn{\gamma} values for given distances and variogram models.
#'
#' @param h numeric. Distance values (in km) for which the gamma values should be fitted.
#' @param model Variogram model: `"spherical"` (default), `"exponential"`
#'   or `"gaussian"`.
#' @param nugget numeric in \eqn{[0, 1)}. Nugget of the trace variogram.
#' @param range_km numeric. Range of the variogram model in km.
#' @param sill numeric in \eqn{(\text{nugget}, 1]}. Sill of the trace variogram.
#'
#' @return numeric vector of the fitted gamma \eqn{\gamma} values in the variogram.
#' Same length as `h`
#'
#' @export
#'
#' @seealso [circular_krige()] and [variogram_axial()]
krige_gamma <- function(h, model, nugget, range_km, sill = 1) {
  r <- h / range_km
  g <- switch(model,
              spherical   = ifelse(r < 1, 1.5 * r - 0.5 * r^3, 1),
              exponential = 1 - exp(-3 * r),
              gaussian    = 1 - exp(-3 * r^2)
  )
  nugget + (sill - nugget) * g
}

#' Empirical trace variogram of axial data, for choosing the model
#'
#' Bins the empirical variogram \eqn{\gamma(h) = \text{mean}(||T_i - T_j||_F^2)}, where
#' `T = ortensor2d(a, norm = TRUE)` for each axis. For axes this equals
#' \eqn{\text{mean}(1 - \cos(2  (a_i - a_j)))}. Use it to choose `nugget` and `range_km`
#' for [circular_krige()].
#'
#' @param x `sf` object with `azi` (degrees).
#' @param breaks numeric. Distance bin edges in km.
#'
#' @returns data.frame with `lag` (km), `gamma` and `n` (pairs per bin).
#'
#' @seealso [circular_krige()] and [krige_gamma()]
#'
#' @export
#'
#' @examples
#' v <- variogram_axial(san_andreas, breaks = seq(0, 1000, 25))
#' v <- v[!is.na(v$gamma), ] # drop empty bins
#'
#' plot(v$lag, v$gamma, pch = 16, las = 1,
#'   xlab = "Distance h (km)", ylab = expression(gamma(h)),
#'   ylim = c(0, 1.05), main = "Trace variogram of SHmax axes")
#' abline(h = 1, lty = 2) # sill: trace of a unit tensors
#'
#' h <- seq(0, 1000, length.out = 300)
#' lines(h, krige_gamma(h, model = "spherical",
#'   nugget = 0.25, range_km = 500, sill = 0.6),
#'   col = "red", lwd = 2)#'
#' lines(h, krige_gamma(h, model = "exponential",
#'   nugget = 0.25, range_km = 500, sill = 0.6),
#'   col = "blue", lwd = 2)
#' lines(h, krige_gamma(h, model = "gaussian",
#'   nugget = 0.25, range_km = 500, sill = 0.6),
#'   col = "green", lwd = 2)
#' legend("bottomright", legend = c("empirical", "spherical model", "exponential", "gaussian"),
#'  pch = c(16, NA, NA, NA), lty = c(NA, 1, 1, 1), col = c("black", "red", "blue", "green"), bty = "n")
variogram_axial <- function(x, breaks = seq(0, 500, 50)) {
  xy <- sf::st_coordinates(x)
  keep <- !is.na(x$azi) & stats::complete.cases(xy)
  x <- x[keep, ]
  xy <- xy[keep, , drop = FALSE]
  if (nrow(x) < 2) stop("At least two observations with valid azimuth and coordinates are needed.", call. = FALSE)

  lat2r <- deg2rad(xy[, 2])
  lon2r <- deg2rad(xy[, 1])
  azi <- x$azi
  n <- nrow(x)
  idx <- which(upper.tri(diag(n)), arr.ind = TRUE)
  i <- idx[, 1]
  j <- idx[, 2]
  h <- .hav_km(lat2r[i], lon2r[i], lat2r[j], lon2r[j])

  g <- 2 * circular_distance(azi[i], azi[j], axial = TRUE, na.rm = FALSE)

  bins <- cut(h, breaks, include.lowest = TRUE)
  data.frame(
    lag = as.numeric(tapply(h, bins, mean)),
    gamma = as.numeric(tapply(g, bins, mean)),
    n = as.integer(table(bins))
  )
}

#' Spatial kriging of SHmax axes after van den Boogaart & Schaeben (2002)
#'
#' Ordinary kriging of the axial orientation tensors, with a distance-only
#' trace variogram (isotropy of measurements). The kriged tensor is decomposed
#' with the eigenvalue method of [ot_eigen2d()]. The per-observation
#' uncertainty, given as a standard deviation in degrees, is used as
#' measurement error on the diagonal of the kriging matrix.
#'
#' @inheritParams PoR_circular_smooth
#' @param R_range numeric value or vector of search radii in km. Only
#'   observations within the radius are used.
#' @param max_data integer. The number of nearest observations used for each
#'   prediction. Default is `50`.
#' @param vgm_model Variogram model: `"spherical"` (default), `"exponential"`
#'   or `"gaussian"`.
#' @param nugget numeric in \eqn{[0, 1)}. Nugget of the trace variogram.
#' @param range_km numeric. Range of the variogram model in km.
#' @param sill numeric in \eqn{(\text{nugget}, 1]}. Sill of the trace variogram.
#'   If `NULL` (default), it is set to \eqn{1 - \bar{R}^2}, where \eqn{\bar{R}}
#'   is the mean resultant length of the doubled angles of `x$azi`.
#' @param use_uncertainty logical. If `TRUE` (default), `unc` (standard
#'   deviation in degrees) is used as measurement error.
#' @param ... Arguments passed to [circular_krige()]
#'
#' @note  Only axial (180&deg; periodic) data are
#'   supported.
#'
#' @returns `sf` object containing
#' \describe{
#' \item{lon,lat}{longitude and latitude in degrees}
#' \item{azi}{kriged axis direction in degrees (0 to 180), the principal
#'   direction from [ot_eigen2d()]}
#' \item{sd}{circular standard deviation of the axis in degrees, computed
#'   from `conf` as `rad2deg(sqrt(-2 log(conf))) / 2`}
#' \item{R}{search radius in km}
#' \item{mdr}{mean distance between grid point and data / R}
#' \item{N}{number of observations used}
#' \item{conf}{confidence in \eqn{[0, 1]}, \eqn{2  \lambda_1 - 1}, i.e. the axial strength
#'   relative to the trace of the kriged tensor}
#' \item{kvar}{kriging variance on the trace-variogram scale}
#' \item{eig1}{largest eigenvalue of the kriged tensor (axial strength), in
#'   \eqn{[0.5, 1]}. The smaller eigenvalue, \eqn{1 = \lambda_1}, is the axial dispersion}
#' }
#'
#' Attribute `"vgm"` records the variogram parameters used.
#'
#' @details The method is the simplified ordinary-kriging form derived by van
#' den Boogaart & Schaeben (2002) for axes in the embedding \eqn{u u^T}. Their
#' isotropy-of-measurements assumption gives the trace variogram used here.
#' The kriged tensor is \eqn{\sum_i \lambda_i  u_i u_i^T}, computed by
#' [ot_eigen2d()] with the kriging weights as weights. The principal
#' direction is the eigenvector of the largest eigenvalue, and the eigenvalues
#' sum to one. Their values give the strength and dispersion of the prediction.
#'
#' Fit `nugget`, `sill`, and `range_km` to the empirical variogram from
#' [variogram_axial()].
# #' The variogram has sill 1, because the trace of the
# #' tensor of a unit vector is 1 (eq. 36 of the paper).
#'
#' @references van den Boogaart, K.G. and Schaeben, H. (2002).
#' Kriging of regionalized directions, axes, and orientations I.
#' Mathematical Geology, 34(5), 479-503.
#'
#' @family smoothing
#' @seealso [variogram_axial()] for computing the empirical variogram of axial data;
#' [ortensor2d()] and [ot_eigen2d()] for the orientation tensor and its eigenvalue decomposition.
#'
#' @name krige
#'
#' @examples
#' data("san_andreas")
#' v <- variogram_axial(san_andreas)
#' res <- circular_krige(san_andreas, R_range = 300,
#'                            vgm_model = "spherical",
#'                            nugget = 0.25, range_km = 1000)
#' head(res)
#' attributes(res)$vgm
NULL

#' @rdname krige
#' @export
circular_krige <- function(x,
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
                           sill = NULL,
                           use_uncertainty = TRUE
) {
  stopifnot(
    inherits(x, "sf"),
    is.numeric(gridsize) && length(gridsize) == 1,
    is.numeric(R_range) && length(R_range) >= 1,
    is.numeric(max_data) && max_data >= min_data,
    is.numeric(nugget) && length(nugget) == 1 && nugget >= 0,
    is.numeric(range_km) && length(range_km) == 1 && range_km > 0,
    is.null(sill) || (is.numeric(sill) && length(sill) == 1 && sill <= 1),
    is.numeric(min_dist_threshold) && min_dist_threshold > 0,
    is.logical(use_uncertainty)
  )
  vgm_model <- match.arg(vgm_model)
  min_data <- as.integer(ceiling(min_data))

  x_coords <- sf::st_coordinates(x)
  lon <- x_coords[, 1]
  lat <- x_coords[, 2]
  azi <- x$azi
  n <- length(azi)
  if (n < min_data) stop("Not enough observations (fewer than `min_data`).", call. = FALSE)

  # sill: 1 - R_bar^2 from the data, where R_bar = lambda_1 - lambda_2 of the
  # orientation tensor of all observations (mean resultant length of 2a)
  if (is.null(sill)) {
    ot_all <- ot_eigen2d(azi, scale = TRUE)
    R_bar <- ot_all$values[1] - ot_all$values[2]
    sill <- 1 - R_bar^2
  }
  if (sill <= nugget) {
    stop("`sill` must be larger than `nugget` (data sill is ", round(sill, 3), ").",
         call. = FALSE)
  }

  lon2r <- deg2rad(lon)
  lat2r <- deg2rad(lat)

  # measurement variance on the trace-variogram scale: an angle sd s gives
  # a variogram contribution of (2 s)^2 between two noisy observations
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
  cols <- c("azi", "sd", "R", "mdr", "N", "conf", "kvar", "eig1", "eig2")

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
      Gm <- krige_gamma(Dm, vgm_model, nugget, range_km, sill)
      diag(Gm) <- sig2[ord]

      # ordinary kriging system with unbiasedness constraint
      A <- rbind(cbind(Gm, 1), c(rep(1, k), 0))
      g0 <- c(krige_gamma(d[ord], vgm_model, nugget, range_km, sill), 1)
      sol <- tryCatch(solve(A, g0), error = function(e) NULL)
      if (is.null(sol)) next
      lam <- sol[seq_len(k)]

      # kriged orientation tensor sum_i lam_i * u_i u_i^T and its eigen
      # decomposition, via the tectonicr eigenvalue method (no doubling)
      ot <- ot_eigen2d(azi[ord], w = lam, scale = TRUE)
      eig1 <- ot$values[1]                    # largest eigenvalue; the other is 1 - eig1
      eig2 <- ot$values[2]                    # smaller eigenvalue
      azi_i <- ot$vectors[1] %% 180           # principal direction, degrees

      conf <- min(max(2 * eig1 - 1, 0), 1)
      sd_i <- if (conf >= 1) {
        0
      } else if (conf <= 0) {
        NA_real_
      } else {
        rad2deg(sqrt(-2 * log(conf))) / 2
      }
      kvar <- sum(lam * g0[seq_len(k)]) + sol[k + 1]

      mat[i, ] <- c(azi_i, sd_i, R_k, mean(d[ord]) / R_k, k, conf, kvar, eig1, eig2)
    }

    res_k <- data.frame(lon = G[, 1], lat = G[, 2], mat)
    results[[kr]] <- res_k[!is.na(res_k$azi), , drop = FALSE]
  }

  out <- do.call(rbind, results)
  out$N <- as.integer(out$N)
  rownames(out) <- NULL
  out <- sf::st_as_sf(out, coords = c("lon", "lat"), crs = sf::st_crs(x), remove = FALSE)
  attr(out, "vgm") <- list(model = vgm_model, nugget = nugget, range_km = range_km,
                           sill = sill, max_data = max_data)
  out
}

#' @rdname krige
#' @export
PoR_circular_krige <- function(x, PoR, grid = NULL, PoR_grid = TRUE,
                                   lon_range = NULL, lat_range = NULL, gridsize = 2.5,
                                   remove_PoR = FALSE, ...) {
  if (!is.null(grid)) {
    lon_range <- lat_range <- gridsize <- NULL
    PoR_grid <- FALSE
  } else {
    if (isFALSE(PoR_grid)) {
      if (is.null(lon_range) || is.null(lat_range)) {
        coords <- sf::st_coordinates(x)
        lon_range <- range(coords[, 1], na.rm = TRUE)
        lat_range <- range(coords[, 2], na.rm = TRUE)
      }

      grid <- sf::st_bbox(
        c(
          xmin = lon_range[1],
          xmax = lon_range[2],
          ymin = lat_range[1],
          ymax = lat_range[2]
        )
      ) |>
        sf::st_make_grid(
          cellsize = gridsize,
          what = "centers",
          offset = c(lon_range[1], lat_range[1])
        )
    }
  }

  grid_PoR <- if (isFALSE(PoR_grid)) {
    sf::st_as_sf(grid, crs = sf::st_crs(x)) |>
      geographical_to_PoR_sf(PoR)
  } else {
    NULL
  }

  x_PoR <- geographical_to_PoR_sf(x, PoR)
  x_PoR_coords <- sf::st_coordinates(x_PoR) |>
    as.data.frame() |>
    dplyr::rename(lat = Y, lon = X)

  azi <- lat <- lon <- lat.PoR <- lon.PoR <- X <- Y <- R <- numeric() # pre allocating:

  x_PoR$lat <- x_PoR_coords$lat
  x_PoR$lon <- x_PoR_coords$lon
  x_PoR$azi <- PoR_shmax(x, PoR)

  int <- circular_krige(x_PoR, grid = grid_PoR, lon_range = lon_range, lat_range = lat_range, gridsize = gridsize, ...) |>
    dplyr::rename(azi.PoR = azi, lat.PoR = lat, lon.PoR = lon) |>
    PoR_to_geographical_sf(PoR)
  int_coords <- sf::st_coordinates(int) |>
    as.data.frame() |>
    dplyr::rename(lat = Y, lon = X)
  int$lat <- int_coords$lat
  int$lon <- int_coords$lon
  int$azi <- PoR2Geo_azimuth(int, PoR)

  if (remove_PoR) int <- dplyr::select(int, !dplyr::ends_with(".PoR"))

  return(int)
}
