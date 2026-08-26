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
#' @param axial Logical. Whether data are uniaxial (`axial=FALSE`)
#' or biaxial (`TRUE`, the default).
# #' @param ... (optional) arguments to [dist_greatcircle()]
#'
#' @importFrom sf st_coordinates st_bbox st_make_grid st_crs st_as_sf
#' @importFrom dplyr group_by mutate filter rename mutate bind_rows select
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
#' @seealso [dist_greatcircle()], [PoR_stress2grid()], [compact_grid()],
#' [circular_mean()], [circular_median()], [circular_sd()], [circular_summary()]
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
    x <- subset(x, !is.na(unc))
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
    x <- subset(x, !is.na(unc))
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
  out <- dplyr::select(out, -c(md, n))
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
#' @seealso [stress2grid()], [compact_grid()]
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
PoR_stress2grid <- function(x, PoR, grid = NULL, PoR_grid = TRUE, lon_range = NULL, lat_range = NULL, gridsize = 2.5, remove_PoR = FALSE, ...) {
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
PoR_stress2grid_stats <- function(x, PoR, grid = NULL, PoR_grid = TRUE, lon_range = NULL, lat_range = NULL, gridsize = 2.5, remove_PoR = FALSE, ...) {
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
#' @seealso [stress2grid()], [PoR_stress2grid()], [kernel_dispersion()],
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
#' @importFrom dplyr group_by mutate
#' @importFrom tidyr drop_na
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
