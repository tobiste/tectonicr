#' Earth's radius in km
#'
#' IERS mean radius of Earth in km (based on WGS 84)
#'
#' @returns numeric value
#'
#' @export
earth_radius <- function() 6371.0087714

#' @keywords internal
wcmean <- function(x, w, f = 2) {
  Z <- sum(w, na.rm = TRUE)
  if (Z != 0) {
    m <- mean_SC(2 * x, w = w, na.rm = TRUE)
    meanR <- sqrt(m["C"]^2 + m["S"]^2)
    # sd_s <- if (meanR > 1) {
    #   0
    # } else {
    #   sqrt(-2 * log(meanR)) / 2
    # }
    # mean_s <- (atan2(m["S"], m["C"]) / 2) %% pi
    # unname(rad2deg(c(mean_s, sd_s)))
    sd_s <- if (is.na(meanR)) NA_real_ else if (meanR > 1) 0 else rad2deg(sqrt(-2 * log(meanR))) / f
    mean_s <- rad2deg(atan2(m["S"], m["C"]) / f) %% (360 / f)
    unname(c(mean_s, sd_s))
  } else {
    c(NA_real_, NA_real_)
  }
}

#' @keywords internal
wcmedian <- function(x, w, f = 2) {
  is.axial <- f == 2
  Z <- sum(w, na.rm = TRUE)
  if (Z > 3) {
    quantiles <- circular_quantiles(x, w, axial = is.axial)
    median_s <- (quantiles[3])
    iqr_s <- deviation_norm(quantiles[4], quantiles[2])
  } else if (Z > 0 && Z <= 3) {
    median_s <- circular_median(x, w, axial = is.axial)
    iqr_s <- ceiling(deviation_norm(max(x), min(x)))
  } else {
    median_s <- iqr_s <- NA_real_
  }
  unname(c(median_s, iqr_s))
}

wprincipal <- function(x, w, f) {
  ot <- ot_eigen2d(x, w, scale = TRUE)
  c(ot$vectors[1], ot$values[2])
}


#' @keywords internal
dist_weighting_linear <- function(R_search, dist_threshold, distij, idp = 0) {
  dist_threshold_scal <- R_search * dist_threshold
  R_search + 1 - pmax(dist_threshold_scal, distij)
}

#' @keywords internal
dist_weighting_inverse <- function(R_search, dist_threshold, distij, idp = 0) {
  dist_threshold_scal <- R_search * dist_threshold
  1 / (pmax(dist_threshold_scal, distij))^idp
}

#' Indices of n smallest values in array
#' @keywords internal
which.nsmallest <- function(x, n) {
  # nsmallest <- utils::head(sort(x), n)
  # which(x %in% nsmallest)
  order(x)[seq_len(min(n, length(x)))]
}

#' Spatial Interpolation of SHmax
#'
#' Stress field interpolation and wavelength analysis using a kernel (weighted)
#' mean/median and standard deviation/IQR of stress data.
#' Parameters can be adjusted to have inverse-distance-weighting (IDW) or
#' nearest-neighbor interpolations (NN).
#'
#' @param x `sf` object containing
#' \describe{
#' \item{azi}{\eqn{\sigma_\text{Hmax}}{SHmax} in degree}
#' \item{unc}{(optional) Uncertainties of SHmax in degree}
#' \item{type}{(optional) Methods used for the determination of the direction
#' of \eqn{\sigma_\text{Hmax}}{SHmax}}
#' }
#' @param grid (optional) Point object of class `sf`.
#' @param lon_range,lat_range (optional) numeric vector specifying the minimum
#' and maximum longitudes and latitudes (ignored if `grid` is specified).
#' @param gridsize numeric. Target spacing of the regular grid in decimal
#' degree. Default is `2.5`. (is ignored if `grid` is specified)
#' @param stat whether the direction of interpolated \eqn{\sigma_\text{Hmax}}{SHmax} is based on the
#' circular mean and standard deviation (`"mean"`, the default), the
#' quasi-circular median and quasi-interquartile range (`"median"`), or the
#' orientation tensor based principal direction and dispersion ("tensor").
#' @param min_data integer. If the number of observations within distance
#' `R_range` is less than `min_data`, a missing value `NA` will be generated.
#' Default is `3` for [stress2grid()] and `4` for [stress2grid_stats()].
#' @param max_data integer. The number of nearest observations that should be
#' used for prediction, where "nearest" is defined in terms of the space of the
#' spatial locations. Default is `Inf`.
#' @param max_sd numeric. Threshold for deviation of direction in degrees;
#' if exceeds, missing values will be generated.
#' @param threshold `r lifecycle::badge("deprecated")` is no
#'   longer supported; use `max_sd` instead.
#' @param min_dist_threshold numeric. Distance threshold for smallest distance
#' of the prediction location to the next observation location.
#' Default is `200` km.
#' @param arte_thres `r lifecycle::badge("deprecated")` is no
#'   longer supported; use `min_dist_threshold` instead.
#' @param dist_weighting Distance weighting method which should be used. One of
#' `"none"`, `"linear"`, or `"inverse"` (the default).
#' @param idp,qp,mp numeric. The weighting power of inverse distance, quality
#' and method (the higher the value, the more weight).
#' Default is `1`. When set to `0`, no weighting is applied. Only effective when
#' `dist_weighting=="inverse"`.
#' @param dist_threshold numeric. Distance weight to prevent overweight of data
#' nearby (0 to 1). Default is `0.1`
#' @param method_weighting logical. If a method weighting should be applied:
#' Default is `FALSE`. If `FALSE`, overwrites `mp`.
#' @param quality_weighting logical. If a quality weighting should be applied:
#' Default is `TRUE`. If `FALSE`, overwrites `qp`.
#' @param R_range numeric value or vector specifying the kernel half-width(s)
#' search radii,
#' i.e. the maximum distance from the prediction location to be used for
#' prediction (in km). Default is `seq(50, 1000, 50)`. If combined with
#' `max_data`, both criteria apply.
#' @param mode logical. Should the circular mode be included in the statistical summary (slow)?
#' @param kappa  numeric. von Mises distribution concentration parameter used
#' for the circular mode. Will be estimated using [est.kappa()] if not provided.
#' @param axial Logical. Whether angles in `x$azi` are uniaxial (`FALSE`)
#' or biaxial (`TRUE`, the default).
#' @param ... (optional) arguments to [circular_summary()]
#'
#' @importFrom sf st_coordinates st_bbox st_make_grid st_crs st_as_sf
#'
#' @returns `sf` object containing
#' \describe{
#' \item{lon,lat}{longitude and latitude in degrees}
#' \item{azi}{Circular mean od median SHmax in degree}
#' \item{sd}{Circular standard deviation or Quasi-IQR on the Circle of
#'  \eqn{\sigma_\text{Hmax}}{SHmax} in degrees}
#' \item{R}{Search radius in km}
#' \item{mdr}{Mean distance between grid point and datapoints per search radius}
#' \item{N}{Number of data points in search radius}
#' }
#' When [stress2grid_stats()], `azi` and `sd` are replaced by the output of
#' [circular_summary()].
#'
#' @details [stress2grid()] is originally based on the MATLAB script
#' "stress2grid" by Ziegler and Heidbach (2019):
#' \url{https://github.com/MorZieg/Stress2Grid}.
#' The tectonicr version has been significantly modified to provide better
#' performance and more flexibility.
#'
#' [stress2grid_stats()] is based on [stress2grid()] but calculates circular
#' summary statistics (see [circular_summary()]).
#'
#' @note Although specialized for stress fields, this spatial interpolation
#' algorithm can be applied to type of data. Please adjust the `axial` argument
#' when applied to different, directional data sets.
#'
#' @seealso [dist_greatcircle()], [compact_grid()],
#' [circular_mean()], [circular_median()], [circular_sd()], [circular_summary()]
#'
#' @family smoothing
#'
#' @references Ziegler, M. and Heidbach, O. (2019).
#' Matlab Script Stress2Grid v1.1. GFZ Data Services. \doi{10.5880/wsm.2019.002}
#'
#' @name stress2grid
#'
#' @examples
#' data("san_andreas")
#'
#' # Inverse Distance Weighting interpolation:
#' stress2grid(san_andreas, stat = "median") |> head()
#'
#' stress2grid(san_andreas, stat = "tensor") |> head()
#'
#' # Nearest Neighbor interpolation:
#' stress2grid(san_andreas, stat = "median", max_data = 5) |> head()
#'
#' \dontrun{
#' stress2grid_stats(san_andreas, mode = TRUE) |> head()
#' }
NULL

#' @rdname stress2grid
#' @export
# stress2grid <- function(x,
#                         stat = c("mean", "median", "tensor"),
#                         grid = NULL,
#                         lon_range = NULL,
#                         lat_range = NULL,
#                         gridsize = 2,
#                         min_data = 3L,
#                         max_data = Inf,
#                         max_sd = Inf,
#                         threshold = deprecated(),
#                         min_dist_threshold = 200,
#                         arte_thres = deprecated(),
#                         method_weighting = FALSE,
#                         quality_weighting = TRUE,
#                         dist_weighting = c("inverse", "linear", "none"),
#                         idp = 1,
#                         qp = 1,
#                         mp = 1,
#                         dist_threshold = 0.1,
#                         R_range = seq(50, 1000, 50),
#                         ...) {
#   stopifnot(
#     inherits(x, "sf"),
#     is.numeric(gridsize) && length(gridsize) == 1,
#     is.numeric(max_sd) | is.infinite(max_sd),
#     is.numeric(max_data) | is.infinite(max_data),
#     is.numeric(min_data) | is.infinite(min_data),
#     max_data >= min_data,
#     is.numeric(min_dist_threshold),
#     is.numeric(dist_threshold),
#     min_dist_threshold > 0 && length(min_dist_threshold) == 1,
#     is.numeric(R_range),
#     is.logical(method_weighting),
#     is.logical(quality_weighting),
#     is.numeric(idp) && length(idp) == 1,
#     is.numeric(qp) && length(qp) == 1,
#     is.numeric(mp) && length(mp) == 1
#   )
#
#   if (lifecycle::is_present(arte_thres)) {
#     lifecycle::deprecate_warn(
#       when = "0.4.6.9002",
#       what = "stress2grid(arte_thres)",
#       details = "Ability to specify arte_thres will be dropped in next release."
#     )
#   }
#   arte_thres <- min_dist_threshold
#
#   if (lifecycle::is_present(threshold)) {
#     lifecycle::deprecate_warn(
#       when = "0.4.6.9002",
#       what = "stress2grid(threshold)",
#       details = "Ability to specify threshold will be dropped in next release."
#     )
#   }
#   threshold <- max_sd
#
#   min_data <- as.integer(ceiling(min_data))
#
#   dist_weighting <- match.arg(dist_weighting)
#   w_distance_fun <- if (dist_weighting == "linear") dist_weighting_linear else dist_weighting_inverse
#
#   stat <- match.arg(stat)
#   stats_fun <- if (stat == "median") wcmedian else if (stat == "tensor") wprincipal else wcmean
#
#   colnames_x <- colnames(x)
#
#   if (quality_weighting & "unc" %in% colnames_x) {
#     x <- subset(x, !is.na(unc))
#   }
#
#   # pre-allocating
#   azi <- x$azi
#   length_azi <- length(azi)
#   unc <- lat <- lon <- numeric(length_azi)
#   type <- character(9)
#   N <- md <- R <- numeric()
#
#
#   if (isFALSE(quality_weighting)) qp <- 0
#   if (isFALSE(method_weighting)) mp <- 0
#   if (dist_weighting == "none") idp <- 0
#
#   # WSM method weighting (from 0 to 5)
#   if ("type" %in% colnames_x) {
#     parse_method <- setNames(
#       c(4, 5, 5, 5, 4, 5, 4, 2, 1) / 5,
#       c("FMS", "FMF", "BO", "DIF", "HF", "GF", "GV", "OC", NA)
#     )
#     w_method <- parse_method[x$type]
#   } else {
#     w_method <- rep(1, length_azi)
#   }
#
#   w_quality <- if ("unc" %in% colnames_x) {
#     weighting(x$unc)
#   } else {
#     rep(1, length_azi)
#   }
#
#   x_coords <- sf::st_coordinates(x)
#
#   datas <- cbind(
#     lon = x_coords[, 1],
#     lat = x_coords[, 2],
#     azi = azi,
#     w_method = ifelse(is.na(w_method), 1 / 5, w_method)^mp,
#     w_quality = w_quality^qp
#   )
#
#   if (is.null(grid)) {
#     # Regular grid
#     if (is.null(lon_range) | is.null(lat_range)) {
#       lon_range <- range(datas[, 1], na.rm = TRUE)
#       lat_range <- range(datas[, 2], na.rm = TRUE)
#     }
#
#     grid <- sf::st_bbox(
#       c(
#         xmin = lon_range[1],
#         xmax = lon_range[2],
#         ymin = lat_range[1],
#         ymax = lat_range[2]
#       ),
#       crs = sf::st_crs("WGS84")
#     ) |>
#       sf::st_make_grid(
#         cellsize = gridsize,
#         what = "centers",
#         offset = c(lon_range[1], lat_range[1])
#       ) |>
#       sf::st_as_sf()
#   }
#   stopifnot(inherits(grid, "sf"), any(sf::st_is(grid, "POINT")))
#   G <- unname(sf::st_coordinates(grid))
#   R_seq <- seq_along(R_range)
#
#   lapply(seq_along(G[, 1]), function(i) {
#     distij <- dist_greatcircle(G[i, 2], G[i, 1], datas[, 2], datas[, 1])
#     if (max_data < Inf) distij <- distij[which.nsmallest(distij, max_data)] # select the `max_data` nearest locations
#
#     # min_dist_thresholdij <- min(c(max(distij), min_dist_threshold))
#
#     if (min(distij) <= min_dist_threshold) {
#       t(vapply(R_seq, function(k) {
#         R_search <- R_range[k]
#         ids_R <- (distij <= R_search) # select those that are in search radius
#         N_in_R <- sum(ids_R)
#         # if(is.null(N_in_R)) N_in_R <- 0L
#
#         if (N_in_R < min_data) {
#           # not enough data within search radius
#           sdSH <- 0
#           meanSH <- md <- NA_real_
#         } else if (N_in_R == 1) {
#           sdSH <- 0
#           meanSH <- datas[ids_R, 3]
#           md <- distij[ids_R]
#         } else {
#           md <- mean(distij[ids_R], na.rm = TRUE)
#
#           # distance weighting
#           w_distance <- w_distance_fun(R_search, dist_threshold, distij[ids_R], idp)
#           w <- w_distance * datas[ids_R, 5] * datas[ids_R, 4]
#
#           # mean value
#           stats <- stats_fun(x = datas[ids_R, 3], w = w)
#           meanSH <- stats[1]
#           sdSH <- stats[2]
#         }
#
#         c(
#           lon = G[i, 1], # lon
#           lat = G[i, 2], # lat
#           azi = meanSH, # azi
#           sd = sdSH, # sd
#           R = R_search, # R_search
#           md = md, # mdr
#           N = N_in_R # N_in_R
#         )
#       }, FUN.VALUE = numeric(7)))
#     }
#   }) |>
#     lapply(as.data.frame) |>
#     dplyr::bind_rows() |>
#     dplyr::mutate(
#       N = as.integer(N),
#       mdr = md / R
#     ) |>
#     dplyr::select(-md) |>
#     # dplyr::filter(!is.na(azi), sd <= max_sd, !is.na(sd)) |>
#     sf::st_as_sf(coords = c("lon", "lat"), crs = sf::st_crs(x), remove = FALSE)
# }
stress2grid <- function(x,
                        stat = c("mean", "median", "tensor"),
                        grid = NULL,
                        lon_range = NULL,
                        lat_range = NULL,
                        gridsize = 2,
                        min_data = 3L,
                        max_data = Inf,
                        max_sd = Inf,
                        threshold = deprecated(),
                        min_dist_threshold = 200,
                        arte_thres = deprecated(),
                        method_weighting = FALSE,
                        quality_weighting = TRUE,
                        dist_weighting = c("inverse", "linear", "none"),
                        idp = 1,
                        qp = 1,
                        mp = 1,
                        dist_threshold = 0.1,
                        R_range = seq(50, 1000, 50),
                        axial = TRUE) {
  stopifnot(
    inherits(x, "sf"),
    is.numeric(gridsize) && length(gridsize) == 1,
    is.numeric(max_sd) | is.infinite(max_sd),
    is.numeric(max_data) | is.infinite(max_data),
    is.numeric(min_data) | is.infinite(min_data),
    max_data >= min_data,
    is.numeric(min_dist_threshold),
    is.numeric(dist_threshold),
    min_dist_threshold > 0 && length(min_dist_threshold) == 1,
    is.numeric(R_range),
    is.logical(method_weighting),
    is.logical(quality_weighting),
    is.numeric(idp) && length(idp) == 1,
    is.numeric(qp) && length(qp) == 1,
    is.numeric(mp) && length(mp) == 1
  )


  if (lifecycle::is_present(arte_thres)) {
    lifecycle::deprecate_warn(
      when = "0.4.6.9002", what = "stress2grid(arte_thres)",
      details = "Ability to specify arte_thres will be dropped in next release."
    )
  }
  arte_thres <- min_dist_threshold
  if (lifecycle::is_present(threshold)) {
    lifecycle::deprecate_warn(
      when = "0.4.6.9002", what = "stress2grid(threshold)",
      details = "Ability to specify threshold will be dropped in next release."
    )
  }
  threshold <- max_sd

  # Periodicity of angles
  f <- if (axial) 2 else 1

  min_data <- as.integer(ceiling(min_data))
  dist_weighting <- match.arg(dist_weighting)
  stat <- match.arg(stat)
  stats_fun <- if (stat == "median") wcmedian else if (stat == "tensor") wprincipal else wcmean
  fast_mean <- identical(stats_fun, wcmean) # only "mean" gets the fully inlined fast path
  w_distance_fun <- if (dist_weighting == "linear") dist_weighting_linear else dist_weighting_inverse

  colnames_x <- colnames(x)
  if (quality_weighting & "unc" %in% colnames_x) {
    x <- x[!is.na(x$unc), ]
  }
  azi <- x$azi
  length_azi <- length(azi)
  if (isFALSE(quality_weighting)) qp <- 0
  if (isFALSE(method_weighting)) mp <- 0
  if (dist_weighting == "none") idp <- 0

  # WSM method weighting (from 0 to 5)
  if ("type" %in% colnames_x) {
    parse_method <- setNames(
      c(4, 5, 5, 5, 4, 5, 4, 2, 1) / 5,
      c("FMS", "FMF", "BO", "DIF", "HF", "GF", "GV", "OC", NA)
    )
    w_method <- parse_method[x$type]
  } else {
    w_method <- rep(1, length_azi)
  }
  w_quality <- if ("unc" %in% colnames_x) weighting(x$unc) else rep(1, length_azi)

  x_coords <- sf::st_coordinates(x)
  datas <- cbind(
    lon = x_coords[, 1],
    lat = x_coords[, 2],
    azi = azi,
    w_method = ifelse(is.na(w_method), 1 / 5, w_method)^mp,
    w_quality = w_quality^qp
  )

  if (is.null(grid)) {
    if (is.null(lon_range) | is.null(lat_range)) {
      lon_range <- range(datas[, "lon"], na.rm = TRUE)
      lat_range <- range(datas[, "lat"], na.rm = TRUE)
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

  ## ------------------------------------------------------------------
  ## optimized computational engine
  ##   1. per-data-point trig (lat/lon, and for stat="mean" also 2*azi)
  ##      is computed ONCE, not once per grid cell
  ##   2. within each grid cell, distances are sorted ONCE and every
  ##      R_range cutoff is located in one vectorized findInterval()
  ##      call, instead of one O(n_data) boolean scan per R value
  ##   3. `datas` is reordered together with `distij` (fixes the
  ##      max_data<Inf row-misalignment bug in the original)
  ##   4. results accumulate as plain numeric matrices and are only
  ##      converted to a data.frame/sf object once, at the very end
  ## ------------------------------------------------------------------
  n_data <- nrow(datas)
  lat2r <- deg2rad(datas[, "lat"])
  lon2r <- deg2rad(datas[, "lon"])
  r_earth <- earth_radius()
  if (fast_mean) {
    cos2azi <- cos(deg2rad(f * datas[, "azi"]))
    sin2azi <- sin(deg2rad(f * datas[, "azi"]))
  }
  R_range_sorted <- sort(R_range)
  R_order <- order(R_range)
  n_R <- length(R_range)
  R_thr <- R_range * dist_threshold
  n_G <- nrow(G)
  results <- vector("list", n_G)

  for (i in seq_len(n_G)) {
    lat1r <- deg2rad(G[i, 2])
    lon1r <- deg2rad(G[i, 1])
    dlon <- lon2r - lon1r
    hd <- sin((lat2r - lat1r) / 2)^2
    hl <- sin(dlon / 2)^2
    hs <- sin((lat2r + lat1r) / 2)^2
    distij <- ahav(hd + (1 - hd - hs) * hl) * r_earth

    idx <- seq_len(n_data)
    if (max_data < Inf && max_data < n_data) idx <- which.nsmallest(distij, max_data)
    ord <- idx[order(distij[idx])]
    distij_s <- distij[ord]

    if (length(distij_s) == 0 || distij_s[1] > min_dist_threshold) next

    cutoffs_sorted <- findInterval(R_range_sorted, distij_s)
    cutoffs <- integer(n_R)
    cutoffs[R_order] <- cutoffs_sorted

    mat <- matrix(NA_real_,
      nrow = n_R, ncol = 7,
      dimnames = list(NULL, c("lon", "lat", "azi", "sd", "R", "md", "N"))
    )
    mat[, 1] <- G[i, 1]
    mat[, 2] <- G[i, 2]
    mat[, 5] <- R_range
    mat[, 7] <- cutoffs

    if (fast_mean) {
      c2_s <- cos2azi[ord]
      s2_s <- sin2azi[ord]
      wq_wm_s <- datas[ord, "w_quality"] * datas[ord, "w_method"]
      for (k in seq_len(n_R)) {
        N_in_R <- cutoffs[k]
        if (N_in_R < min_data) {
          mat[k, 4] <- 0
        } else if (N_in_R == 1) {
          mat[k, 3] <- rad2deg(atan2(s2_s[1], c2_s[1]) / f) %% (360 / f)
          mat[k, 4] <- 0
          mat[k, 6] <- distij_s[1]
        } else {
          sl <- seq_len(N_in_R)
          d_sub <- distij_s[sl]
          mat[k, 6] <- sum(d_sub) / N_in_R
          w_distance <- if (idp == 0) rep(1, N_in_R) else 1 / (pmax(R_thr[k], d_sub))^idp
          w <- w_distance * wq_wm_s[sl]
          Z <- sum(w)
          Cc <- sum(w * c2_s[sl]) / Z
          Ss <- sum(w * s2_s[sl]) / Z
          meanR <- sqrt(Cc * Cc + Ss * Ss)
          mat[k, 4] <- if (is.na(meanR)) {
            NA_real_
          } else if (meanR > 1) {
            0
          } else {
            rad2deg(sqrt(-2 * log(meanR))) / f
          }
          mat[k, 3] <- rad2deg(atan2(Ss, Cc) / f) %% (360 / f)
        }
      }
    } else {
      # median / tensor: general path -- still gets the sort-once +
      # findInterval speedup and the datas/distij alignment fix, just
      # not the azimuth-trig precompute (wcmedian/wprincipal need the
      # raw per-point values, not summable trig terms)
      datas_s <- datas[ord, , drop = FALSE]
      for (k in seq_len(n_R)) {
        N_in_R <- cutoffs[k]
        if (N_in_R < min_data) {
          mat[k, 4] <- 0
        } else if (N_in_R == 1) {
          mat[k, 3] <- datas_s[1, "azi"]
          mat[k, 4] <- 0
          mat[k, 6] <- distij_s[1]
        } else {
          sl <- seq_len(N_in_R)
          d_sub <- distij_s[sl]
          mat[k, 6] <- sum(d_sub) / N_in_R
          w_distance <- w_distance_fun(R_range[k], dist_threshold, d_sub, idp)
          w <- w_distance * datas_s[sl, "w_quality"] * datas_s[sl, "w_method"]
          st <- stats_fun(x = datas_s[sl, "azi"], w = w, f = f)
          mat[k, 3] <- st[1]
          mat[k, 4] <- st[2]
        }
      }
    }
    results[[i]] <- mat
  }

  out <- as.data.frame(do.call(rbind, results))
  out$N <- as.integer(out$N)
  out$mdr <- out$md / out$R
  out$md <- NULL
  sf::st_as_sf(out, coords = c("lon", "lat"), crs = sf::st_crs(x), remove = FALSE)
}


#' @rdname stress2grid
#' @export
# stress2grid_stats <- function(x,
#                               grid = NULL,
#                               lon_range = NULL,
#                               lat_range = NULL,
#                               gridsize = 2,
#                               min_data = 4L,
#                               max_data = Inf,
#                               threshold = deprecated(),
#                               min_dist_threshold = 200,
#                               arte_thres = deprecated(),
#                               method_weighting = FALSE,
#                               quality_weighting = TRUE,
#                               dist_weighting = c("inverse", "linear", "none"),
#                               idp = 1,
#                               qp = 1,
#                               mp = 1,
#                               dist_threshold = 0.1,
#                               R_range = seq(50, 1000, 50),
#                               mode = FALSE,
#                               kappa = 10,
#                               ...) {
#   stopifnot(
#     inherits(x, "sf"),
#     is.numeric(gridsize),
#     is.numeric(min_dist_threshold),
#     min_dist_threshold > 0,
#     is.numeric(max_data) | is.infinite(max_data),
#     is.numeric(min_data) | is.infinite(min_data),
#     max_data >= min_data,
#     is.numeric(dist_threshold),
#     is.numeric(R_range),
#     is.logical(method_weighting),
#     is.logical(quality_weighting),
#     is.numeric(idp),
#     is.numeric(qp),
#     is.numeric(mp),
#     is.logical(mode)
#   )
#
#   if (lifecycle::is_present(arte_thres)) {
#     lifecycle::deprecate_warn(
#       when = "0.4.6.9002",
#       what = "stress2grid_stats(arte_thres)",
#       details = "Ability to specify arte_thres will be dropped in next release."
#     )
#   }
#   arte_thres <- min_dist_threshold
#
#   if (lifecycle::is_present(threshold)) {
#     lifecycle::deprecate_warn(
#       when = "0.4.6.9002",
#       what = "stress2grid_stats(threshold)",
#       details = "Ability to specify threshold will be dropped in next release."
#     )
#   }
#   # threshold <- max_sd
#
#
#
#   min_data <- as.integer(ceiling(min_data))
#
#   dist_weighting <- match.arg(dist_weighting)
#   if (dist_weighting == "linear") {
#     w_distance_fun <- dist_weighting_linear
#   } else {
#     w_distance_fun <- dist_weighting_inverse
#   }
#
#   colnames_x <- colnames(x)
#
#   if (quality_weighting & "unc" %in% colnames_x) {
#     x <- subset(x, !is.na(unc))
#   }
#
#   # pre-allocating
#   azi <- x$azi
#   length_azi <- length(azi)
#   unc <- lat <- lon <- numeric(length_azi)
#   type <- character(9)
#   # num_r <- length(R_range)
#   n <- N <- md <- R <- numeric()
#
#
#   if (isFALSE(quality_weighting)) qp <- 0
#   if (isFALSE(method_weighting)) mp <- 0
#   if (dist_weighting == "none") idp <- 0
#
#   # WSM method weighting (from 0 to 5)
#   if ("type" %in% colnames_x) {
#     parse_method <- setNames(
#       c(4, 5, 5, 5, 4, 5, 4, 2, 1) / 5,
#       c("FMS", "FMF", "BO", "DIF", "HF", "GF", "GV", "OC", NA)
#     )
#     w_method <- parse_method[x$type]
#   } else {
#     w_method <- rep(1, length_azi)
#   }
#
#   w_quality <- if ("unc" %in% colnames_x) {
#     weighting(x$unc)
#   } else {
#     rep(1, length_azi)
#   }
#
#   x_coords <- sf::st_coordinates(x)
#
#   datas <- cbind(
#     lon = x_coords[, 1],
#     lat = x_coords[, 2],
#     azi = azi,
#     w_method = ifelse(is.na(w_method), 1 / 5, w_method)^mp,
#     w_quality = w_quality^qp
#   )
#
#   if (is.null(grid)) {
#     # Regular grid
#     if (is.null(lon_range) || is.null(lat_range)) {
#       lon_range <- range(datas[, 1], na.rm = TRUE)
#       lat_range <- range(datas[, 2], na.rm = TRUE)
#     }
#
#     grid <- sf::st_bbox(
#       c(
#         xmin = lon_range[1],
#         xmax = lon_range[2],
#         ymin = lat_range[1],
#         ymax = lat_range[2]
#       ),
#       crs = sf::st_crs("WGS84")
#     ) |>
#       sf::st_make_grid(
#         cellsize = gridsize,
#         what = "centers",
#         offset = c(lon_range[1], lat_range[1])
#       ) |>
#       sf::st_as_sf()
#   }
#   stopifnot(inherits(grid, "sf"), any(sf::st_is(grid, "POINT")))
#   G <- sf::st_coordinates(grid)
#   # num_G <- nrow(G)
#
#   # r <- R <- N <- n <- numeric(num_G)
#   R_seq <- seq_along(R_range)
#   # nR <- length(R_seq)
#
#   cols <- c(
#     "lon", "lat", "n", "mean", "sd", "var",
#     "25%", "quasi-median", "75%", "median", "CI",
#     "skewness", "kurtosis", "meanR", "R", "md", "N"
#   )
#   if (mode) {
#     cols <- append(cols, "mode", after = 10)
#   }
#
#   lapply(seq_along(G[, 1]), function(i) {
#     distij <- dist_greatcircle(G[i, 2], G[i, 1], datas[, 2], datas[, 1], ...)
#     distij <- distij[which.nsmallest(distij, max_data)] # select the `max_data` nearest locations
#
#     if (min(distij) <= min_dist_threshold) {
#       t(vapply(R_seq, function(k) {
#         R_search <- R_range[k]
#         ids_R <- (distij <= R_search) # select those that are in search radius
#         N_in_R <- sum(ids_R)
#
#         if (N_in_R < min_data) {
#           # not enough data within search radius
#           stats <- rep(NA_real_, length(cols) - 5)
#           md <- NA
#         } else if (N_in_R == 1) {
#           stats <- rep(NA_real_, length(cols) - 5)
#           stats[2] <- datas[ids_R, 3]
#           md <- distij[ids_R]
#         } else {
#           md <- mean(distij[ids_R], na.rm = TRUE)
#
#           # distance weighting
#           w_distance <- w_distance_fun(R_search, dist_threshold, distij[ids_R], idp)
#
#           w <- w_distance * datas[ids_R, 5] * datas[ids_R, 4]
#
#           # mean value
#           stats <- circular_summary(x = datas[ids_R, 3], w = w, axial = TRUE, mode = mode, kappa = kappa, na.rm = TRUE) |> unname()
#         }
#         c(
#           lon = G[i, 1], # lon
#           lat = G[i, 2], # lat
#           stats,
#           R = R_search, # R_search
#           md = md, # mdr
#           N = N_in_R # N_in_R
#         )
#       }, FUN.VALUE = numeric(length(cols))))
#     }
#   }) |>
#     lapply(as.data.frame) |>
#     dplyr::bind_rows() |>
#     setNames(cols) |>
#     dplyr::mutate(N = as.integer(N), mdr = md / R) |>
#     dplyr::select(-c(md, n)) |>
#     sf::st_as_sf(coords = c("lon", "lat"), crs = sf::st_crs(x), remove = FALSE)
# }
stress2grid_stats <- function(x,
                              grid = NULL,
                              lon_range = NULL,
                              lat_range = NULL,
                              gridsize = 2,
                              min_data = 4L,
                              max_data = Inf,
                              threshold = deprecated(),
                              min_dist_threshold = 200,
                              arte_thres = deprecated(),
                              method_weighting = FALSE,
                              quality_weighting = TRUE,
                              dist_weighting = c("inverse", "linear", "none"),
                              idp = 1,
                              qp = 1,
                              mp = 1,
                              dist_threshold = 0.1,
                              R_range = seq(50, 1000, 50),
                              mode = FALSE,
                              kappa = 10,
                              axial = TRUE,
                              ...) {
  stopifnot(
    inherits(x, "sf"),
    is.numeric(gridsize),
    is.numeric(min_dist_threshold),
    min_dist_threshold > 0,
    is.numeric(max_data) | is.infinite(max_data),
    is.numeric(min_data) | is.infinite(min_data),
    max_data >= min_data,
    is.numeric(dist_threshold),
    is.numeric(R_range),
    is.logical(method_weighting),
    is.logical(quality_weighting),
    is.numeric(idp),
    is.numeric(qp),
    is.numeric(mp),
    is.logical(mode)
  )
  if (lifecycle::is_present(arte_thres)) {
    lifecycle::deprecate_warn(
      when = "0.4.6.9002", what = "stress2grid_stats(arte_thres)",
      details = "Ability to specify arte_thres will be dropped in next release."
    )
  }
  arte_thres <- min_dist_threshold
  if (lifecycle::is_present(threshold)) {
    lifecycle::deprecate_warn(
      when = "0.4.6.9002", what = "stress2grid_stats(threshold)",
      details = "Ability to specify threshold will be dropped in next release."
    )
  }
  # threshold <- max_sd   # (no max_sd parameter on this function -- unchanged from original)

  min_data <- as.integer(ceiling(min_data))
  dist_weighting <- match.arg(dist_weighting)
  w_distance_fun <- if (dist_weighting == "linear") dist_weighting_linear else dist_weighting_inverse

  colnames_x <- colnames(x)
  if (quality_weighting & "unc" %in% colnames_x) {
    x <- x[!is.na(x$unc), ]
  }
  azi <- x$azi
  length_azi <- length(azi)
  if (isFALSE(quality_weighting)) qp <- 0
  if (isFALSE(method_weighting)) mp <- 0
  if (dist_weighting == "none") idp <- 0

  if ("type" %in% colnames_x) {
    parse_method <- setNames(
      c(4, 5, 5, 5, 4, 5, 4, 2, 1) / 5,
      c("FMS", "FMF", "BO", "DIF", "HF", "GF", "GV", "OC", NA)
    )
    w_method <- parse_method[x$type]
  } else {
    w_method <- rep(1, length_azi)
  }
  w_quality <- if ("unc" %in% colnames_x) weighting(x$unc) else rep(1, length_azi)

  x_coords <- sf::st_coordinates(x)
  datas <- cbind(
    lon = x_coords[, 1],
    lat = x_coords[, 2],
    azi = azi,
    w_method = ifelse(is.na(w_method), 1 / 5, w_method)^mp,
    w_quality = w_quality^qp
  )

  if (is.null(grid)) {
    if (is.null(lon_range) || is.null(lat_range)) {
      lon_range <- range(datas[, "lon"], na.rm = TRUE)
      lat_range <- range(datas[, "lat"], na.rm = TRUE)
    }
    grid <- sf::st_bbox(
      c(xmin = lon_range[1], xmax = lon_range[2], ymin = lat_range[1], ymax = lat_range[2]),
      crs = sf::st_crs("WGS84")
    ) |>
      sf::st_make_grid(cellsize = gridsize, what = "centers", offset = c(lon_range[1], lat_range[1])) |>
      sf::st_as_sf()
  }
  stopifnot(inherits(grid, "sf"), any(sf::st_is(grid, "POINT")))
  G <- sf::st_coordinates(grid)

  cols <- c(
    "lon", "lat", "n", "mean", "sd", "var",
    "25%", "quasi-median", "75%", "median", "CI",
    "skewness", "kurtosis", "meanR", "R", "md", "N"
  )
  if (mode) cols <- append(cols, "mode", after = 10)
  n_stat_cols <- length(cols) - 5 # everything between lat and R: n,mean,sd,var,...(,mode)

  ## ------------------------------------------------------------------
  ## optimized computational engine (see stress2grid() for the rationale):
  ##  - per-data-point distance trig computed once, not once per grid cell
  ##  - distances sorted once per grid cell; every R_range cutoff located
  ##    via one vectorized findInterval() call instead of n_R boolean scans
  ##  - datas reordered together with distij (fixes the max_data<Inf
  ##    row-misalignment bug present in the original)
  ##  - results accumulate as numeric matrices, converted to a data.frame
  ##    only once at the very end
  ## circular_summary() itself is left untouched and called once per
  ## (grid cell, R_range value) exactly as before -- it still needs the
  ## real per-point values (quantiles aren't summable), so it isn't
  ## inlined the way wcmean() was in stress2grid().
  ## ------------------------------------------------------------------
  n_data <- nrow(datas)
  lat2r <- deg2rad(datas[, "lat"])
  lon2r <- deg2rad(datas[, "lon"])
  r_earth <- earth_radius()
  R_range_sorted <- sort(R_range)
  R_order <- order(R_range)
  n_R <- length(R_range)
  n_G <- nrow(G)
  results <- vector("list", n_G)

  for (i in seq_len(n_G)) {
    lat1r <- deg2rad(G[i, 2])
    lon1r <- deg2rad(G[i, 1])
    dlon <- lon2r - lon1r
    hd <- sin((lat2r - lat1r) / 2)^2
    hl <- sin(dlon / 2)^2
    hs <- sin((lat2r + lat1r) / 2)^2
    distij <- ahav(hd + (1 - hd - hs) * hl) * r_earth

    idx <- seq_len(n_data)
    if (max_data < Inf && max_data < n_data) idx <- which.nsmallest(distij, max_data)
    ord <- idx[order(distij[idx])]
    distij_s <- distij[ord]

    if (length(distij_s) == 0 || distij_s[1] > min_dist_threshold) next

    cutoffs_sorted <- findInterval(R_range_sorted, distij_s)
    cutoffs <- integer(n_R)
    cutoffs[R_order] <- cutoffs_sorted
    datas_s <- datas[ord, , drop = FALSE]

    mat <- matrix(NA_real_,
      nrow = n_R, ncol = n_stat_cols + 5,
      dimnames = list(NULL, cols)
    )
    mat[, "lon"] <- G[i, 1]
    mat[, "lat"] <- G[i, 2]
    mat[, "R"] <- R_range
    mat[, "N"] <- cutoffs

    for (k in seq_len(n_R)) {
      N_in_R <- cutoffs[k]
      if (N_in_R < min_data) {
        # stats stay NA, md stays NA
      } else if (N_in_R == 1) {
        mat[k, "mean"] <- datas_s[1, "azi"]
        mat[k, "md"] <- distij_s[1]
      } else {
        sl <- seq_len(N_in_R)
        d_sub <- distij_s[sl]
        mat[k, "md"] <- sum(d_sub) / N_in_R
        w_distance <- w_distance_fun(R_range[k], dist_threshold, d_sub, idp)
        w <- w_distance * datas_s[sl, "w_quality"] * datas_s[sl, "w_method"]
        st <- circular_summary(
          x = datas_s[sl, "azi"], w = w, axial = axial,
          mode = mode, kappa = kappa, na.rm = TRUE
        ) |> unname()
        mat[k, 3:(3 + n_stat_cols - 1)] <- st
      }
    }
    results[[i]] <- mat
  }

  out <- as.data.frame(do.call(rbind, results))
  out$N <- as.integer(out$N)
  out$mdr <- out$md / out$R
  out$md <- NULL
  out$n <- NULL
  sf::st_as_sf(out, coords = c("lon", "lat"), crs = sf::st_crs(x), remove = FALSE)
}


#' Spatial Interpolation of SHmax in PoR Coordinate Reference System
#'
#' Stress field and wavelength analysis in PoR system and back-transformed
#'
#' @param x \code{sf} object containing
#' \describe{
#' \item{azi}{\eqn{\sigma_\text{Hmax}}{SHmax} in degree}
#' \item{unc}{Uncertainties of \eqn{\sigma_\text{Hmax}}{SHmax} in degree}
#' \item{type}{Methods used for the determination of the orientation of \eqn{\sigma_\text{Hmax}}{SHmax}}
#' }
#' @param PoR Pole of Rotation. `data.frame` or object of class
#' \code{"euler.pole"} containing the geographical coordinates of the Euler pole
#' @param grid (optional) Point object of class `sf`.
#' @param PoR_grid logical. Whether the grid should be generated based on the
#' coordinate range in the PoR (`TRUE`, the default) CRS or the geographical CRS
#' (`FALSE`). Is ignored if `grid` is specified.
#' @param lon_range,lat_range (optional) numeric vector specifying the minimum
#' and maximum longitudes and latitudes (are ignored if `grid` is specified).
#' @param gridsize Numeric. Target spacing of the regular grid in decimal
#' degree. Default is `2.5` (is ignored if `grid` is specified)
#' @param remove_PoR logical. Whether PoR azimuths and coordinates will be
#' removed from final output or not (the default.)
#' @param ... Arguments passed to [stress2grid()]
#'
#' @description The data is transformed into the PoR system before the
#' interpolation. The interpolation grid is returned in geographical coordinates
#'  and azimuths.
#'
#' @importFrom dplyr rename group_by
#' @importFrom sf st_coordinates st_as_sf st_bbox st_make_grid
#'
#' @returns \code{sf} object containing
#' \describe{
#' \item{lon,lat}{longitude and latitude in geographical CRS (in degrees)}
#' \item{lon.PoR,lat.PoR}{longitude and latitude in PoR CRS (in degrees).
#' Only if `remove_PoR=TRUE`}
#' \item{azi}{geographical mean \eqn{\sigma_\text{Hmax}}{SHmax} in degree}
#' \item{azi.PoR}{PoR mean \eqn{\sigma_\text{Hmax}}{SHmax} in degree. Only if `remove_PoR=TRUE`}
#' \item{sd}{Standard deviation of \eqn{\sigma_\text{Hmax}}{SHmax} in degrees}
#' \item{R}{Search radius in km}
#' \item{mdr}{Mean distance of datapoints per search radius}
#' \item{N}{Number of data points in search radius}
#' }
#'
#' @seealso [compact_grid()]
#' @family smoothing
#'
#' @name PoR_stress2grid
#'
#' @examples
#' data("san_andreas")
#' data("nuvel1")
#' PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
#' PoR_stress2grid(san_andreas, PoR) |> head()
#'
#' \dontrun{
#' PoR_stress2grid_stats(san_andreas, PoR, mode = TRUE) |> head()
#' }
NULL

#' @rdname PoR_stress2grid
#' @export
PoR_stress2grid <- function(x, PoR, grid = NULL, PoR_grid = TRUE,
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

  int <- stress2grid(x_PoR, grid = grid_PoR, lon_range = lon_range, lat_range = lat_range, gridsize = gridsize, ...) |>
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

#' @rdname PoR_stress2grid
#' @export
PoR_stress2grid_stats <- function(x, PoR, grid = NULL, PoR_grid = TRUE,
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
    sf::st_as_sf(grid) |>
      geographical_to_PoR_sf(PoR)
  } else {
    NULL
  }

  x_PoR <- geographical_to_PoR_sf(x, PoR)
  x_PoR_coords <- sf::st_coordinates(x_PoR) |>
    as.data.frame() |>
    dplyr::rename(lat = Y, lon = X)

  # binding global variables
  azi <- lat <- lon <- lat.PoR <- lon.PoR <- X <- Y <- R <- numeric() # pre allocating:
  `25%` <- `75%` <- `25%.PoR` <- `75%.PoR` <- median.PoR <- mean.PoR <- mode.PoR <- `quasi-median` <- `quasi-median.PoR` <- `median` <- NULL

  x_PoR$lat <- x_PoR_coords$lat
  x_PoR$lon <- x_PoR_coords$lon
  x_PoR$azi <- PoR_shmax(x, PoR)

  int <- stress2grid_stats(x_PoR, grid = grid_PoR, lon_range = lon_range, lat_range = lat_range, gridsize = gridsize, ...) |>
    dplyr::rename(
      mean.PoR = mean, `25%.PoR` = `25%`, `quasi-median.PoR` = `quasi-median`,
      `75%.PoR` = `75%`, median.PoR = median,
      "mode.PoR" = dplyr::matches("mode"),
      lat.PoR = lat, lon.PoR = lon
    ) |>
    PoR_to_geographical_sf(PoR)
  int_coords <- sf::st_coordinates(int) |>
    as.data.frame() |>
    dplyr::rename(lat = Y, lon = X)
  int$lat <- int_coords$lat
  int$lon <- int_coords$lon
  int$mean <- PoR2Geo_azimuth(int |> rename(azi.PoR = mean.PoR), PoR)
  int$`25%` <- PoR2Geo_azimuth(int |> rename(azi.PoR = `25%.PoR`), PoR)
  int$`quasi-median` <- PoR2Geo_azimuth(int |> rename(azi.PoR = `quasi-median.PoR`), PoR)
  int$`75%` <- PoR2Geo_azimuth(int |> rename(azi.PoR = `75%.PoR`), PoR)
  int$median <- PoR2Geo_azimuth(int |> rename(azi.PoR = median.PoR), PoR)
  if ("mode.PoR" %in% colnames(int)) int$mode <- PoR2Geo_azimuth(int |> rename(azi.PoR = mode.PoR), PoR)

  if (remove_PoR) int <- dplyr::select(int, !dplyr::ends_with(".PoR"))

  return(int)
}


#' Compact Smoothed Stress Field
#'
#' Filter smoothed stress field containing a range of search radii or kernel
#' half widths to find shortest wavelength (R) with the least circular sd. or
#' dispersion (or any statistic) for each coordinate, respectively.
#'
#' @param x output of [stress2grid()], [PoR_stress2grid()],
#' [stress2grid_stats()], or [kernel_dispersion()]
#' @param type character. Type of the grid `x`. Either `"stress"` (when input
#' is [stress2grid()] or [PoR_stress2grid()]) or `"dispersion"` (when input
#' is [kernel_dispersion()]).
#' @param ... `<tidy-select>` One unquoted expression separated by
#' commas. Variable names can be used as if they were positions in the data
#' frame. Variable must be a column in `x`.
#' @param FUN function is used to aggregate the data using the search radius
#' `R`. Default is [min()].
#' @returns \code{sf} object
#'
#' @importFrom dplyr ungroup mutate select left_join
#' @importFrom tidyr drop_na
#' @importFrom sf st_as_sf
#' @importFrom stats aggregate
#' @seealso [stress2grid()], [stress2grid_hansen()], [PoR_stress2grid()], [kernel_dispersion()],
#' [stress2grid_stats()], [dplyr::dplyr_tidy_select()]
#'
#' @name compact-grid
#'
#' @examples
#' data("san_andreas")
#' res <- stress2grid(san_andreas)
#' compact_grid(res) |> head()
#'
#' \dontrun{
#' res2 <- stress2grid_stats(san_andreas)
#' compact_grid2(res2, var, FUN = min)
#' }
NULL

#' @rdname compact-grid
#' @export
compact_grid <- function(x, type = c("stress", "dispersion")) {
  var <- character()
  type <- match.arg(type)

  if (type == "stress") {
    var <- "azi"
  } else {
    var <- "stat"
  }
  compact_grid2(x, var, FUN = min)
}

#' @rdname compact-grid
#' @export
compact_grid2 <- function(x, ..., FUN = min) {
  lon <- lat <- R <- numeric()
  group <- character()

  data <- x |>
    tidyr::drop_na(...) |>
    dplyr::mutate(group = paste(lon, lat))

  aggregate(R ~ group, data, FUN, na.rm = TRUE) |>
    dplyr::left_join(data, by = c("group", "R")) |>
    dplyr::select(-group) |>
    sf::st_as_sf()
}


#' Adaptive Kernel Dispersion
#'
#' Stress field and wavelength analysis using circular dispersion
#' (or other statistical estimators for dispersion)
#'
#' @param x `sf` object containing
#' \describe{
#' \item{azi}{the observed \eqn{\sigma_\text{Hmax}}{SHmax} in degree}
#' \item{unc}{(optional) Uncertainties of ibserved SHmax in degree}
#' \item{type}{(optional) Methods used for the determination of the direction
#' of \eqn{\sigma_\text{Hmax}}{SHmax}}
#' \item{prd}{the predicted \eqn{\sigma_\text{Hmax}}{SHmax} in degree}
#' }
#' @inheritParams stress2grid
#' @param stat The measurement of dispersion to be calculated. Either
#' `"dispersion"` (default), `"nchisq"`, or `"rayleigh"` for circular
#' dispersion, normalized Chi-squared test statistic, or Rayleigh test
#' statistic.
#' @param stat_threshold numeric. Generates missing values when the kernel
#' `stat` value exceeds this threshold. Default is `Inf`.
#' @param ... arguments passed to `stat` functions [weighted_rayleigh()] or
#' [circular_dispersion()]
#'
#' @importFrom sf st_coordinates st_bbox st_make_grid st_crs st_as_sf
#'
#' @returns
#' \code{sf} object containing
#' \describe{
#' \item{lon,lat}{longitude and latitude in degree}
#' \item{stat}{output of function defined in `stat`}
#' \item{R}{The rearch radius in km.}
#' \item{mdr}{Mean distance of datapoints per search radius}
#' \item{N}{Number of data points in search radius}
#' }
#'
#' @seealso [circular_dispersion()], [norm_chisq()], [weighted_rayleigh()]
#'
#' @note `dispersion_grid()` was renamed to `kernel_dispersion()` to create
#'  a more consistent API.
#'
#' @name kernel_dispersion
#'
#' @examples
#' data("nuvel1")
#' PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
#' san_andreas_por <- data2PoR(san_andreas, PoR)
#' san_andreas_por$prd <- 135
#' kernel_dispersion(san_andreas_por) |> head()
NULL

#' @rdname kernel_dispersion
#' @export
# kernel_dispersion <- function(x,
#                               stat = c("dispersion", "nchisq", "rayleigh"),
#                               grid = NULL,
#                               lon_range = NULL,
#                               lat_range = NULL,
#                               gridsize = 2.5,
#                               min_data = 3L,
#                               max_data = Inf,
#                               min_dist_threshold = 200,
#                               dist_threshold = 0.1,
#                               stat_threshold = Inf,
#                               R_range = seq(100, 2000, 100),
#                               ...) {
#   stopifnot(
#     inherits(x, "sf"), is.numeric(gridsize), is.numeric(min_data), is.numeric(min_dist_threshold),
#     min_dist_threshold > 0, is.numeric(dist_threshold), is.numeric(R_range),
#     is.numeric(stat_threshold) | is.infinite(stat_threshold),
#     max_data >= min_data
#   )
#   stat <- match.arg(stat)
#   min_data <- as.integer(ceiling(min_data))
#   stat <- match.arg(stat)
#
#   N <- md <- R <- numeric()
#
#   # pre-allocating
#   azi <- x$azi
#   length_azi <- length(azi)
#   colnames_x <- colnames(x)
#   unc <- lat <- lon <- prd <- numeric(length_azi)
#   type <- character(9)
#
#   # num_r <- length(R_range)
#
#   x_coords <-
#     sf::st_coordinates(x) |>
#     as.data.frame()
#
#   datas <- cbind(
#     lon = x_coords$X,
#     lat = x_coords$Y,
#     azi = x$azi,
#     unc = x$unc,
#     prd = x$prd
#   )
#
#   if (is.null(grid)) {
#     # Regular grid
#     if (is.null(lon_range) || is.null(lat_range)) {
#       lon_range <- range(datas[, "lon"], na.rm = TRUE)
#       lat_range <- range(datas[, "lat"], na.rm = TRUE)
#     }
#
#     grid <- sf::st_bbox(
#       c(
#         xmin = lon_range[1],
#         xmax = lon_range[2],
#         ymin = lat_range[1],
#         ymax = lat_range[2]
#       ),
#       crs = sf::st_crs("WGS84")
#     ) |>
#       sf::st_make_grid(
#         cellsize = gridsize,
#         what = "centers",
#         offset = c(lon_range[1], lat_range[1])
#       ) |>
#       sf::st_as_sf()
#   }
#   stopifnot(inherits(grid, "sf"), any(sf::st_is(grid, "POINT")))
#   G <- unname(sf::st_coordinates(grid))
#
#   # R <- N <- numeric(nrow(G))
#   R_seq <- seq_along(R_range)
#   # nR <- length(R_seq)
#
#   # SH <- matrix(nrow = 0, ncol = 6, dimnames = list(NULL, c("lon", "lat", "stat", "R", "md", "N")))
#
#   lapply(seq_along(G[, 1]), function(i) {
#     # for (i in seq_along(G[, 1])) {
#     distij <- dist_greatcircle(G[i, 2], G[i, 1], datas[, "lat"], datas[, "lon"], ...)
#     if (max_data < Inf) distij <- distij[which.nsmallest(distij, max_data)] # select the `max_data` nearest locations
#
#     if (min(distij) <= min_dist_threshold) {
#       t(vapply(R_seq, function(k) {
#         R_search <- R_range[k]
#         ids_R <- which(distij <= R_search)
#         N_in_R <- length(ids_R)
#
#         if (N_in_R < min_data) {
#           # not enough data within search radius
#           y <- md <- NA
#         } else if (N_in_R == 1) {
#           y <- NA
#           md <- distij[ids_R]
#         } else {
#           md <- mean(distij[ids_R], na.rm = TRUE)
#
#           if (stat == "nchisq") {
#             y <- norm_chisq(datas[ids_R, "azi"], prd = datas[ids_R, "prd"], datas[ids_R, "unc"])
#           } else if (stat == "rayleigh") {
#             y <- weighted_rayleigh(datas[ids_R, "azi"], mu = datas[ids_R, "prd"], w = weighting(datas[ids_R, "unc"]), ...)$statistic
#           } else {
#             y <- circular_dispersion(datas[ids_R, "azi"], y = datas[ids_R, "prd"], w = weighting(datas[ids_R, "unc"]), ...)
#           }
#         }
#
#         c(
#           lon = G[i, 1],
#           lat = G[i, 2],
#           stat = y,
#           R = R_search,
#           md = md,
#           N = N_in_R
#         )
#       }, FUN.VALUE = numeric(6)))
#     }
#   }) |>
#     lapply(as.data.frame) |>
#     dplyr::bind_rows() |>
#     dplyr::mutate(
#       N = as.integer(N),
#       mdr = md / R,
#       stat = ifelse(stat >= stat_threshold, NA, stat)
#     ) |>
#     dplyr::select(-md) |>
#     sf::st_as_sf(coords = c("lon", "lat"), crs = sf::st_crs(x), remove = FALSE)
# }
kernel_dispersion <- function(x,
                              stat = c("dispersion", "nchisq", "rayleigh"),
                              grid = NULL,
                              lon_range = NULL,
                              lat_range = NULL,
                              gridsize = 2.5,
                              min_data = 3L,
                              max_data = Inf,
                              min_dist_threshold = 200,
                              dist_threshold = 0.1,
                              stat_threshold = Inf,
                              R_range = seq(100, 2000, 100),
                              ...) {
  stopifnot(
    inherits(x, "sf"), is.numeric(gridsize), is.numeric(min_data), is.numeric(min_dist_threshold),
    min_dist_threshold > 0, is.numeric(dist_threshold), is.numeric(R_range),
    is.numeric(stat_threshold) | is.infinite(stat_threshold),
    max_data >= min_data
  )
  stat <- match.arg(stat)
  min_data <- as.integer(ceiling(min_data))

  azi <- x$azi
  length_azi <- length(azi)
  colnames_x <- colnames(x)

  x_coords <- sf::st_coordinates(x) |> as.data.frame()
  datas <- cbind(
    lon = x_coords$X,
    lat = x_coords$Y,
    azi = x$azi,
    unc = x$unc,
    prd = x$prd
  )

  if (is.null(grid)) {
    if (is.null(lon_range) || is.null(lat_range)) {
      lon_range <- range(datas[, "lon"], na.rm = TRUE)
      lat_range <- range(datas[, "lat"], na.rm = TRUE)
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

  ## ------------------------------------------------------------------
  ## optimized computational engine (see stress2grid() for the rationale):
  ##  - per-data-point distance trig computed once, not once per grid cell
  ##  - distances sorted once per grid cell; every R_range cutoff located
  ##    via one vectorized findInterval() call instead of n_R boolean/which() scans
  ##  - datas reordered together with distij -- this is the actual fix:
  ##    the original's `distij <- distij[which.nsmallest(distij, max_data)]`
  ##    shrinks distij but NEVER touches `datas`, so its later
  ##    `ids_R <- which(distij <= R_search); datas[ids_R, ...]` uses indices
  ##    valid for the SHORTENED distij against the FULL, unshortened datas
  ##    matrix -- silently pulling the wrong rows whenever max_data < n_data.
  ##  - results accumulate as numeric matrices, converted to a data.frame
  ##    only once at the very end
  ## norm_chisq()/weighted_rayleigh()/circular_dispersion() are left
  ## untouched and still called once per (grid cell, R_range value),
  ## just on a correctly-aligned, cheaply-sliced subset.
  ## ------------------------------------------------------------------
  n_data <- nrow(datas)
  lat2r <- deg2rad(datas[, "lat"])
  lon2r <- deg2rad(datas[, "lon"])
  r_earth <- earth_radius()
  R_range_sorted <- sort(R_range)
  R_order <- order(R_range)
  n_R <- length(R_range)
  n_G <- nrow(G)
  results <- vector("list", n_G)

  for (i in seq_len(n_G)) {
    lat1r <- deg2rad(G[i, 2])
    lon1r <- deg2rad(G[i, 1])
    dlon <- lon2r - lon1r
    hd <- sin((lat2r - lat1r) / 2)^2
    hl <- sin(dlon / 2)^2
    hs <- sin((lat2r + lat1r) / 2)^2
    distij <- ahav(hd + (1 - hd - hs) * hl) * r_earth
    # NOTE: if you rely on non-default `...` args to dist_greatcircle() (e.g. a
    # custom `r` or `method`), pass them through to this haversine call, or
    # fall back to `distij <- dist_greatcircle(G[i,2], G[i,1], datas[,"lat"], datas[,"lon"], ...)`
    # here (losing the precomputed-trig speedup but keeping full generality).

    idx <- seq_len(n_data)
    if (max_data < Inf && max_data < n_data) idx <- which.nsmallest(distij, max_data)
    ord <- idx[order(distij[idx])]
    distij_s <- distij[ord]

    if (length(distij_s) == 0 || distij_s[1] > min_dist_threshold) next

    cutoffs_sorted <- findInterval(R_range_sorted, distij_s)
    cutoffs <- integer(n_R)
    cutoffs[R_order] <- cutoffs_sorted
    datas_s <- datas[ord, , drop = FALSE]

    mat <- matrix(NA_real_,
      nrow = n_R, ncol = 6,
      dimnames = list(NULL, c("lon", "lat", "stat", "R", "md", "N"))
    )
    mat[, "lon"] <- G[i, 1]
    mat[, "lat"] <- G[i, 2]
    mat[, "R"] <- R_range
    mat[, "N"] <- cutoffs

    for (k in seq_len(n_R)) {
      N_in_R <- cutoffs[k]
      if (N_in_R < min_data) {
        # stat, md stay NA
      } else if (N_in_R == 1) {
        mat[k, "md"] <- distij_s[1]
      } else {
        sl <- seq_len(N_in_R)
        mat[k, "md"] <- sum(distij_s[sl]) / N_in_R
        azi_sub <- datas_s[sl, "azi"]
        prd_sub <- datas_s[sl, "prd"]
        unc_sub <- datas_s[sl, "unc"]
        y <- if (stat == "nchisq") {
          norm_chisq(azi_sub, prd = prd_sub, unc_sub)
        } else if (stat == "rayleigh") {
          weighted_rayleigh(azi_sub, mu = prd_sub, w = weighting(unc_sub), ...)$statistic
        } else {
          circular_dispersion(azi_sub, y = prd_sub, w = weighting(unc_sub), ...)
        }
        mat[k, "stat"] <- y
      }
    }
    results[[i]] <- mat
  }

  out <- as.data.frame(do.call(rbind, results))
  out$N <- as.integer(out$N)
  out$mdr <- out$md / out$R
  out$stat <- ifelse(out$stat >= stat_threshold, NA, out$stat)
  out$md <- NULL
  sf::st_as_sf(out, coords = c("lon", "lat"), crs = sf::st_crs(x), remove = FALSE)
}

#' @rdname kernel_dispersion
#' @export
#' @keywords internal
dispersion_grid <- function(...) {
  lifecycle::deprecate_warn(" 0.2.97", "dispersion_grid()", "kernel_dispersion")
  kernel_dispersion(...)
}








## Smoothing of SHmax orientations after Hansen & Mount (1990) and
## Müller et al. (2003), implemented in the stress2grid() style.
##
## Two stages:
##   1. .hm_fit():   iteratively estimate a smoothed axis for every observation
##                   (dominant eigenvector of the 2x2 objective matrix, written
##                   in doubled-angle form so the 180-degree periodicity is exact)
##   2. prediction:  weighted circular mean of the fitted axes at each grid node
##                   (reuses wcmean(), so it matches stat = "mean")

#' Haversine distance in km (vectorised, radians in)
#' @keywords internal
.hav_km <- function(lat1, lon1, lat2, lon2) {
  hd <- sin((lat2 - lat1) / 2)^2
  hl <- sin((lon2 - lon1) / 2)^2
  hs <- sin((lat2 + lat1) / 2)^2
  ahav(hd + (1 - hd - hs) * hl) * earth_radius()
}

#' Kernel bandwidth. FSR: fixed radius A. FNN: distance to the (R+1)-th
#' nearest other point (Hansen & Mount 1990, eq. 7).
#' @keywords internal
.bandwidth <- function(d, method, R, A) {
  if (method == "fsr") return(A)
  if (length(d) > R) sort(d, partial = R + 1)[R + 1] else Inf
}

#' Sum of `v` per group `g` into a vector of length n (zero where absent)
#' @keywords internal
.group_sum <- function(v, g, n) {
  out <- numeric(n)
  if (length(v)) {
    s <- rowsum(v, g, reorder = TRUE)
    out[as.integer(rownames(s))] <- s[, 1]
  }
  out
}

#' Sparse neighbour triplets (i, j, w) between observations, j != i, with
#' rows rescaled to sum to one (Müller et al. 2003, eq. 9).
#' @keywords internal
.hm_neighbours <- function(lon2r, lat2r, method, R, A) {
  n <- length(lon2r)
  ii <- jj <- ww <- vector("list", n)
  for (i in seq_len(n)) {
    others <- seq_len(n)[-i]
    d <- .hav_km(lat2r[i], lon2r[i], lat2r[others], lon2r[others])
    h <- .bandwidth(d, method, R, A)
    keep <- which(d < h)
    ii[[i]] <- rep.int(i, length(keep))
    jj[[i]] <- others[keep]
    ww[[i]] <- (1 - (d[keep] / h)^3)^3
  }
  nb <- list(i = unlist(ii), j = unlist(jj), w = unlist(ww))
  G <- .group_sum(nb$w, nb$i, n)
  nb$w <- nb$w / G[nb$i]
  nb
}

#' Iterative fit of smoothed doubled angles at the observation points
#'
#' Objective per observation i (Hansen & Mount 1990, eq. 5; Müller et al. 2003,
#' eq. 3):  M_i = q_i delta_i y_i y_i^T + lambda * sum_j w~_ij delta_j f_j f_j^T.
#' The identity part of y y^T does not change the dominant eigenvector, so only
#' the doubled-angle components (cos 2a, sin 2a) are carried.
#'
#' @param robust_scale minimum robustness scale in chord units
#'   (2 * sin(deg / 2)), see stress2grid_hansen().
#' @keywords internal
.hm_fit <- function(nb, n, th2, q, lambda, n_iter, tol, robust_scale = 0) {
  c0 <- cos(th2)
  s0 <- sin(th2)
  phi2 <- th2                       # start from the observations
  delta <- rep(1, n)                # robustness weights
  for (m in seq_len(n_iter)) {
    c <- cos(phi2)
    s <- sin(phi2)
    C <- q * delta * c0 + lambda * .group_sum(nb$w * (delta * c)[nb$j], nb$i, n)
    S <- q * delta * s0 + lambda * .group_sum(nb$w * (delta * s)[nb$j], nb$i, n)
    phi2_new <- atan2(S, C)

    # observations without any weight keep their previous angle
    zero <- C^2 + S^2 == 0
    phi2_new[zero] <- phi2[zero]

    # axial misfit between observed and fitted axis, as chord length
    dif <- atan2(sin(th2 - phi2_new), cos(th2 - phi2_new))
    r <- 2 * sin(abs(dif) / 4)

    # bisquare robustness weights; the scale is floored so the cutoff cannot
    # collapse onto the fit (which would remove real local rotations)
    s_med <- max(stats::median(r), robust_scale)
    delta <- if (s_med > 1e-6) {
      ifelse(r < 6 * s_med, (1 - (r / (6 * s_med))^2)^2, 0)
    } else {
      rep(1, n)
    }

    step <- max(abs(atan2(sin(phi2_new - phi2), cos(phi2_new - phi2))))
    phi2 <- phi2_new
    if (step < tol) break
  }
  list(phi2 = phi2, delta = delta)
}

#' Smoothed SHmax on a grid after Hansen & Mount (1990) / Müller et al. (2003)
#'
#' Spatial smoothing of SHmax orientations with a distance-weighted, robust
#' smoother. `lambda` balances smoothness against fit to the data. Two kernels
#' are available: `"fsr"` (fixed search radius) and `"fnn"` (fixed number of
#' nearest neighbours).
#'
#' @inheritParams stress2grid
#' @param method Kernel type: `"fsr"` (fixed search radius) or `"fnn"`
#'   (fixed number of nearest neighbours).
#' @param R_range numeric. For `method = "fsr"`: search radii in km. For
#'   `method = "fnn"`: number of nearest neighbours. Can also be a vector of different radii.
#' @param lambda numeric. Smoothness parameter. `lambda > 1` favours
#'   smoothness; `lambda < 1` favours fidelity to the data.
#' @param n_iter integer. Maximum number of fitting iterations.
#' @param tol numeric. Convergence tolerance on the doubled angle (radians).
#' @param robust_scale_deg numeric. Minimum robustness scale in degrees. The
#'   robust cutoff is 6 times the larger of this value and the median residual,
#'   so observations further than about `6 * robust_scale_deg` from the fit are
#'   downweighted. Protects real local rotations (e.g. along faults) from being
#'   removed as outliers. Default is `10`.
#' @param min_sq numeric. Minimum sum of quality x robustness weights of the
#'   contributing observations; grid nodes below are dropped.
#'  @param qp numeric. The weighting power of the data quality (the higher the
#'  value, the more weight). Default is `1`. When set to `0`, no weighting is applied.
#'
#' @note  Only axial (180&deg; periodic) data are
#'   supported.
#'
#' @returns `sf` object with columns
#' \describe{
#' \item{lon, lat}{longitude and latitude in degrees}
#' \item{azi}{smoothed direction in degrees (0 to 180)}
#' \item{sd}{circular standard deviation in degrees}
#' \item{R}{search radius in km (for FNN, the effective radius at each node)}
#' \item{mdr}{mean distance / R}
#' \item{N}{number of contributing observations}
#' \item{sq}{sum of quality x robustness weights}
#' }
#' Attribute `"fit"` is a list (one element per `R_range`) with the
#' observation-level fit: `lon`, `lat`, `azi_obs`, `azi_fit`, `q` (quality,
#' scaled so the best observation is 1) and `delta` (robustness weight; low
#' values mark downweighted outliers).
#'
#' @details Based on Hansen & Mount (1990) and
#' Müller, Wehrle, Hettel, Sperner & Fuchs (2003). The
#' predicted direction at a grid node is the weighted circular mean of the
#' fitted observation axes, computed with the same routine as `stat = "mean"`.
#'
#' @family smoothing
#' @seealso [compact_grid()]
#'
#' @references
#' Hansen, K. M., & Mount, V. S. (1990). Smoothing and extrapolation of crustal
#' stress orientation measurements. Journal of Geophysical Research: Solid Earth,
#' 95(B2), 1155–1165. https://doi.org/10.1029/JB095iB02p01155
#'
#' Müller, B., Wehrle, V., Hettel, S., Sperner, B., & Fuchs, K. (2003). A new
#' method for smoothing orientated data and its application to stress data.
#' Geological Society, London, Special Publications, 209(1), 107–126.
#' https://doi.org/10.1144/GSL.SP.2003.209.01.11
#'
#' @examples
#' \dontrun{
#' data("san_andreas")
#' res <- stress2grid_hansen(san_andreas, method = "fsr",
#'                           R_range = 300, lambda = 0.5)
#' head(res)
#' }
#' @export
stress2grid_hansen <- function(x,
                               grid = NULL,
                               lon_range = NULL,
                               lat_range = NULL,
                               gridsize = 2,
                               method = c("fsr", "fnn"),
                               R_range = 300,
                               lambda = 0.5,
                               n_iter = 20,
                               tol = 1e-6,
                               robust_scale_deg = 10,
                               min_data = 3L,
                               min_sq = 0,
                               min_dist_threshold = Inf,
                               quality_weighting = FALSE,
                               qp = 1) {
  stopifnot(
    inherits(x, "sf"),
    is.numeric(gridsize) && length(gridsize) == 1,
    is.numeric(R_range) && length(R_range) >= 1,
    is.numeric(lambda) && length(lambda) == 1 && lambda >= 0,
    is.numeric(n_iter) && n_iter >= 1,
    is.numeric(robust_scale_deg) && robust_scale_deg >= 0,
    is.numeric(min_dist_threshold) && min_dist_threshold > 0,
    is.logical(quality_weighting),
    is.numeric(qp) && length(qp) == 1
  )
  # if (!isTRUE(axial)) {
  #   stop("stress2grid_hansen() requires axial data (`axial = TRUE`).", call. = FALSE)
  # }
  method <- match.arg(method)
  min_data <- as.integer(ceiling(min_data))
  if (method == "fnn") R_range <- as.integer(round(R_range))

  # robustness scale floor, in chord units (chord = 2 sin(angle / 2))
  robust_scale <- 2 * sin(deg2rad(robust_scale_deg) / 2)

  # quality weights: weighting() sums to 1 over the data, so rescale to a
  # maximum of 1 (best observation = 1, Müller et al. A = 1, B = 0.75, C = 0.5)
  colnames_x <- colnames(x)
  if (quality_weighting && "unc" %in% colnames_x) {
    x <- x[!is.na(x$unc), ]
  }
  if (isFALSE(quality_weighting)) qp <- 0
  w_quality <- if ("unc" %in% colnames_x) weighting(x$unc) else rep(1, nrow(x))
  q <- w_quality^qp
  q <- q / max(q, na.rm = TRUE)

  x_coords <- sf::st_coordinates(x)
  lon <- x_coords[, 1]
  lat <- x_coords[, 2]
  azi <- x$azi
  n <- length(azi)
  if (n < min_data) stop("Not enough observations (fewer than `min_data`).", call. = FALSE)

  lon2r <- deg2rad(lon)
  lat2r <- deg2rad(lat)
  th2 <- deg2rad(2 * azi)     # doubled observed angles

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
  fits <- vector("list", length(R_range))

  for (k in seq_along(R_range)) {
    R_k <- R_range[k]
    Rn <- if (method == "fnn") R_k else NA_integer_
    A  <- if (method == "fsr") R_k else NA_real_

    # stage 1: fit smoothed axes at the observation points
    nb <- .hm_neighbours(lon2r, lat2r, method, R = Rn, A = A)
    fit <- .hm_fit(nb, n, th2, q, lambda, n_iter, tol, robust_scale)
    qd <- q * fit$delta                       # weights entering the prediction
    azi_fit <- rad2deg(fit$phi2 / 2) %% 180
    fits[[k]] <- data.frame(lon = lon, lat = lat, azi_obs = azi,
                            azi_fit = azi_fit, q = q, delta = fit$delta)

    # stage 2: predict on the grid (weighted circular mean of fitted axes)
    mat <- matrix(NA_real_, nrow = nG, ncol = 6,
                  dimnames = list(NULL, c("azi", "sd", "R", "mdr", "N", "sq")))
    for (i in seq_len(nG)) {
      d <- .hav_km(glat1r[i], glon1r[i], lat2r, lon2r)
      if (min(d) > min_dist_threshold) next
      h <- .bandwidth(d, method, Rn, A)
      if (!is.finite(h)) next
      keep <- which(d < h)
      N <- length(keep)
      if (N < min_data) next
      sq <- sum(qd[keep])
      if (sq < min_sq) next
      w <- (1 - (d[keep] / h)^3)^3 * qd[keep]
      st <- wcmean(azi_fit[keep], w)
      R_eff <- if (method == "fsr") A else h
      mat[i, ] <- c(st[1], st[2], R_eff, mean(d[keep]) / R_eff, N, sq)
    }

    res_k <- data.frame(lon = G[, 1], lat = G[, 2], mat)
    results[[k]] <- res_k[!is.na(res_k$azi), , drop = FALSE]
  }

  out <- do.call(rbind, results)
  out$N <- as.integer(out$N)
  rownames(out) <- NULL
  out <- sf::st_as_sf(out, coords = c("lon", "lat"), crs = sf::st_crs(x), remove = FALSE)
  attr(out, "fit") <- setNames(fits, paste0("R=", R_range))
  out
}
