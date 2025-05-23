#' Earth's radius in km
#'
#' IERS mean radius of Earth in km (based on WGS 84)
#'
#' @returns numeric value
#'
#' @export
earth_radius <- function() 6371.0087714

#' @keywords internal
wcmean <- function(x, w) {
  Z <- sum(w, na.rm = TRUE)
  if (Z != 0) {
    m <- mean_SC(2 * x, w = w, na.rm = TRUE)
    meanR <- sqrt(m["C"]^2 + m["S"]^2)
    sd_s <- if (meanR > 1) {
      0
    } else {
      sqrt(-2 * log(meanR))
    }
    mean_s <- atan2(m["S"], m["C"]) / 2
    unname(rad2deg(c(mean_s, sd_s)) %% 180)
  } else {
    c(NA, NA)
  }
}

#' @keywords internal
wcmedian <- function(x, w) {
  Z <- sum(w, na.rm = TRUE)
  if (Z > 3) {
    quantiles <- circular_quantiles(x, w)
    median_s <- (quantiles[3])
    iqr_s <- deviation_norm(quantiles[4], quantiles[2])
  } else if (Z > 0 & Z <= 3) {
    median_s <- circular_median(x, w)
    iqr_s <- ceiling(deviation_norm(max(x), min(x)))
  } else {
    median_s <- iqr_s <- NA
  }
  unname(c(median_s, iqr_s))
}

#' @keywords internal
dist_weighting_linear <- function(R_search, dist_threshold, distij, idp = 0) {
  dist_threshold_scal <- R_search * dist_threshold
  R_search + 1 - max(dist_threshold_scal, distij)
}

#' @keywords internal
dist_weighting_inverse <- function(R_search, dist_threshold, distij, idp = 0) {
  dist_threshold_scal <- R_search * dist_threshold
  1 / (max(dist_threshold_scal, distij))^idp
}

#' Indices of n smallest values in array
#' @keywords internal
which.nsmallest <- function(x, n) {
  # n <- min(c(length(x), Inf))
  # nsmallest <- sort(x)[1:n]
  nsmallest <- utils::head(sort(x), n)

  which(x %in% nsmallest)
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
#' \item{azi}{SHmax in degree}
#' \item{unc}{(optional) Uncertainties of SHmax in degree}
#' \item{type}{(optional) Methods used for the determination of the direction
#' of SHmax}
#' }
#' @param grid (optional) Point object of class `sf`.
#' @param lon_range,lat_range (optional) numeric vector specifying the minimum
#' and maximum longitudes and latitudes (ignored if `grid` is specified).
#' @param gridsize numeric. Target spacing of the regular grid in decimal
#' degree. Default is `2.5`. (is ignored if `grid` is specified)
#' @param stat whether the direction of interpolated SHmax is based on the
#' circular mean and standard deviation (`"mean"`, the default) or the
#' quasi-circular median and quasi-interquartile range (`"median"`).
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
#' @param ... (optional) arguments to [dist_greatcircle()]
#'
#' @importFrom sf st_coordinates st_bbox st_make_grid st_crs st_as_sf
#' @importFrom dplyr group_by mutate filter rename mutate bind_rows select
#'
#' @returns `sf` object containing
#' \describe{
#' \item{lon,lat}{longitude and latitude in degrees}
#' \item{azi}{Circular mean od median SHmax in degree}
#' \item{sd}{Circular standard deviation or Quasi-IQR on the Circle of SHmax in degrees}
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
#' # Nearest Neighbor interpolation:
#' stress2grid(san_andreas, stat = "median", max_data = 5) |> head()
#'
#' \dontrun{
#' stress2grid_stats(san_andreas) |> head()
#' }
NULL

#' @rdname stress2grid
#' @export
stress2grid <- function(x,
                        stat = c("mean", "median"),
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
                        ...) {
  stopifnot(
    inherits(x, "sf"),
    is.numeric(gridsize),
    is.numeric(max_sd) | is.infinite(max_sd),
    is.numeric(max_data) | is.infinite(max_data),
    is.numeric(min_data) | is.infinite(min_data),
    max_data >= min_data,
    is.numeric(min_dist_threshold),
    is.numeric(dist_threshold),
    min_dist_threshold > 0,
    is.numeric(R_range),
    is.logical(method_weighting),
    is.logical(quality_weighting),
    is.numeric(idp),
    is.numeric(qp),
    is.numeric(mp)
  )

  if (lifecycle::is_present(arte_thres)) {
    lifecycle::deprecate_warn(
      when = "0.4.6.9002",
      what = "stress2grid(arte_thres)",
      details = "Ability to specify arte_thres will be dropped in next release."
    )
  }
  arte_thres <- min_dist_threshold

  if (lifecycle::is_present(threshold)) {
    lifecycle::deprecate_warn(
      when = "0.4.6.9002",
      what = "stress2grid(threshold)",
      details = "Ability to specify threshold will be dropped in next release."
    )
  }
  threshold <- max_sd

  min_data <- as.integer(ceiling(min_data))

  dist_weighting <- match.arg(dist_weighting)
  if (dist_weighting == "linear") {
    w_distance_fun <- dist_weighting_linear
  } else {
    w_distance_fun <- dist_weighting_inverse
  }

  stat <- match.arg(stat)
  if (stat == "median") {
    stats_fun <- wcmedian
  } else {
    stats_fun <- wcmean
  }

  colnames_x <- colnames(x)

  if (quality_weighting & "unc" %in% colnames_x) {
    x <- subset(x, !is.na(unc))
  }

  # pre-allocating
  azi <- x$azi
  length_azi <- length(azi)
  unc <- lat <- lon <- numeric(length_azi)
  type <- character(9)
  N <- md <- R <- numeric()


  if (!quality_weighting) qp <- 0
  if (!method_weighting) mp <- 0
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

  w_quality <- if ("unc" %in% colnames_x) {
    1 / x$unc
  } else {
    rep(1, length_azi)
  }

  x_coords <- sf::st_coordinates(x)

  datas <- cbind(
    lon = x_coords[, 1],
    lat = x_coords[, 2],
    azi = azi,
    w_method = ifelse(is.na(w_method), 1 / 5, w_method)^mp,
    w_quality = w_quality^qp
  )

  if (is.null(grid)) {
    # Regular grid
    if (is.null(lon_range) || is.null(lat_range)) {
      lon_range <- range(datas[, 1], na.rm = TRUE)
      lat_range <- range(datas[, 2], na.rm = TRUE)
    }

    grid <- sf::st_bbox(
      c(
        xmin = lon_range[1],
        xmax = lon_range[2],
        ymin = lat_range[1],
        ymax = lat_range[2]
      ),
      crs = sf::st_crs("WGS84")
    ) |>
      sf::st_make_grid(
        cellsize = gridsize,
        what = "centers",
        offset = c(lon_range[1], lat_range[1])
      ) |>
      sf::st_as_sf()
  }
  stopifnot(inherits(grid, "sf"), any(sf::st_is(grid, "POINT")))
  G <- unname(sf::st_coordinates(grid))
  R_seq <- seq_along(R_range)

  lapply(seq_along(G[, 1]), function(i) {
    distij <- dist_greatcircle(G[i, 2], G[i, 1], datas[, 2], datas[, 1])
    if (max_data < Inf) distij <- distij[which.nsmallest(distij, max_data)] # select the `max_data` nearest locations

    # min_dist_thresholdij <- min(c(max(distij), min_dist_threshold))

    if (min(distij) <= min_dist_threshold) {
      t(vapply(R_seq, function(k) {
        R_search <- R_range[k]
        ids_R <- (distij <= R_search) # select those that are in search radius
        N_in_R <- sum(ids_R)
        # if(is.null(N_in_R)) N_in_R <- 0L

        if (N_in_R < min_data) {
          # not enough data within search radius
          sdSH <- 0
          meanSH <- md <- NA
        } else if (N_in_R == 1) {
          sdSH <- 0
          meanSH <- datas[ids_R, 3]
          md <- distij[ids_R]
        } else {
          md <- mean(distij[ids_R], na.rm = TRUE)

          # distance weighting
          w_distance <- w_distance_fun(R_search, dist_threshold, distij[ids_R], idp)

          w <- w_distance * datas[ids_R, 5] * datas[ids_R, 4]

          # mean value
          stats <- stats_fun(x = datas[ids_R, 3], w = w)
          meanSH <- stats[1]
          sdSH <- stats[2]
        }

        c(
          lon = G[i, 1], # lon
          lat = G[i, 2], # lat
          azi = meanSH, # azi
          sd = sdSH, # sd
          R = R_search, # R_search
          md = md, # mdr
          N = N_in_R # N_in_R
        )
      }, FUN.VALUE = numeric(7)))
    }
  }) |>
    lapply(as.data.frame) |>
    dplyr::bind_rows() |>
    dplyr::mutate(
      N = as.integer(N),
      mdr = md / R
    ) |>
    dplyr::select(-md) |>
    dplyr::filter(!is.na(azi), sd <= max_sd, !is.na(sd)) |>
    sf::st_as_sf(coords = c("lon", "lat"), crs = sf::st_crs(x), remove = FALSE)
}


#' @rdname stress2grid
#' @export
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
      when = "0.4.6.9002",
      what = "stress2grid_stats(arte_thres)",
      details = "Ability to specify arte_thres will be dropped in next release."
    )
  }
  arte_thres <- min_dist_threshold

  if (lifecycle::is_present(threshold)) {
    lifecycle::deprecate_warn(
      when = "0.4.6.9002",
      what = "stress2grid_stats(threshold)",
      details = "Ability to specify threshold will be dropped in next release."
    )
  }
  # threshold <- max_sd



  min_data <- as.integer(ceiling(min_data))

  dist_weighting <- match.arg(dist_weighting)
  if (dist_weighting == "linear") {
    w_distance_fun <- dist_weighting_linear
  } else {
    w_distance_fun <- dist_weighting_inverse
  }

  colnames_x <- colnames(x)

  if (quality_weighting & "unc" %in% colnames_x) {
    x <- subset(x, !is.na(unc))
  }

  # pre-allocating
  azi <- x$azi
  length_azi <- length(azi)
  unc <- lat <- lon <- numeric(length_azi)
  type <- character(9)
  # num_r <- length(R_range)
  n <- N <- md <- R <- numeric()


  if (!quality_weighting) qp <- 0
  if (!method_weighting) mp <- 0
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

  w_quality <- if ("unc" %in% colnames_x) {
    1 / x$unc
  } else {
    rep(1, length_azi)
  }

  x_coords <- sf::st_coordinates(x)

  datas <- cbind(
    lon = x_coords[, 1],
    lat = x_coords[, 2],
    azi = azi,
    w_method = ifelse(is.na(w_method), 1 / 5, w_method)^mp,
    w_quality = w_quality^qp
  )

  if (is.null(grid)) {
    # Regular grid
    if (is.null(lon_range) || is.null(lat_range)) {
      lon_range <- range(datas[, 1], na.rm = TRUE)
      lat_range <- range(datas[, 2], na.rm = TRUE)
    }

    grid <- sf::st_bbox(
      c(
        xmin = lon_range[1],
        xmax = lon_range[2],
        ymin = lat_range[1],
        ymax = lat_range[2]
      ),
      crs = sf::st_crs("WGS84")
    ) |>
      sf::st_make_grid(
        cellsize = gridsize,
        what = "centers",
        offset = c(lon_range[1], lat_range[1])
      ) |>
      sf::st_as_sf()
  }
  stopifnot(inherits(grid, "sf"), any(sf::st_is(grid, "POINT")))
  G <- sf::st_coordinates(grid)
  # num_G <- nrow(G)

  # r <- R <- N <- n <- numeric(num_G)
  R_seq <- seq_along(R_range)
  # nR <- length(R_seq)

  cols <- c(
    "lon", "lat", "n", "mean", "sd", "var",
    "25%", "quasi-median", "75%", "median", "CI",
    "skewness", "kurtosis", "meanR", "R", "md", "N"
  )
  if (mode) {
    cols <- append(cols, "mode", after = 10)
  }

  lapply(seq_along(G[, 1]), function(i) {
    distij <- dist_greatcircle(G[i, 2], G[i, 1], datas[, 2], datas[, 1], ...)
    distij <- distij[which.nsmallest(distij, max_data)] # select the `max_data` nearest locations

    if (min(distij) <= min_dist_threshold) {
      t(vapply(R_seq, function(k) {
        R_search <- R_range[k]
        ids_R <- (distij <= R_search) # select those that are in search radius
        N_in_R <- sum(ids_R)

        if (N_in_R < min_data) {
          # not enough data within search radius
          stats <- rep(NA, length(cols) - 5)
          md <- NA
        } else if (N_in_R == 1) {
          stats <- rep(NA, length(cols) - 5)
          stats[2] <- datas[ids_R, 3]
          md <- distij[ids_R]
        } else {
          md <- mean(distij[ids_R], na.rm = TRUE)

          # distance weighting
          w_distance <- w_distance_fun(R_search, dist_threshold, distij[ids_R], idp)

          w <- w_distance * datas[ids_R, 5] * datas[ids_R, 4]

          # mean value
          stats <- circular_summary(x = datas[ids_R, 3], w = w, axial = TRUE, mode = mode, kappa = kappa, na.rm = TRUE) |> unname()
        }
        c(
          lon = G[i, 1], # lon
          lat = G[i, 2], # lat
          stats,
          R = R_search, # R_search
          md = md, # mdr
          N = N_in_R # N_in_R
        )
      }, FUN.VALUE = numeric(length(cols))))
    }
  }) |>
    lapply(as.data.frame) |>
    dplyr::bind_rows() |>
    setNames(cols) |>
    dplyr::mutate(N = as.integer(N), mdr = md / R) |>
    dplyr::select(-c(md, n)) |>
    sf::st_as_sf(coords = c("lon", "lat"), crs = sf::st_crs(x), remove = FALSE)
}




#' Spatial Interpolation of SHmax in PoR Coordinate Reference System
#'
#' Stress field and wavelength analysis in PoR system and back-transformed
#'
#' @param x \code{sf} object containing
#' \describe{
#' \item{azi}{SHmax in degree}
#' \item{unc}{Uncertainties of SHmax in degree}
#' \item{type}{Methods used for the determination of the orientation of SHmax}
#' }
#' @param PoR Pole of Rotation. \code{"data.frame"} or object of class
#' \code{"euler.pole"} containing the geographical coordinates of the Euler pole
#' @param grid (optional) Point object of class \code{sf}.
#' @param PoR_grid logical. Whether the grid should be generated based on the
#' coordinate range in the PoR (`TRUE`, the default) CRS or the geographical CRS
#' (`FALSE`). Is ignored if `grid` is specified.
#' @param lon_range,lat_range (optional) numeric vector specifying the minimum
#' and maximum longitudes and latitudes (are ignored if `"grid"` is specified).
#' @param gridsize Numeric. Target spacing of the regular grid in decimal
#' degree. Default is 2.5 (is ignored if `grid` is specified)
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
#' \item{lon.PoR,lat.PoR}{longitude and latitude in PoR CRS (in degrees)}
#' \item{azi}{geographical mean SHmax in degree}
#' \item{azi.PoR}{PoR mean SHmax in degree}
#' \item{sd}{Standard deviation of SHmax in degrees}
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
PoR_stress2grid <- function(x, PoR, grid = NULL, PoR_grid = TRUE, lon_range = NULL, lat_range = NULL, gridsize = 2.5, ...) {
  if (!is.null(grid)) {
    lon_range <- lat_range <- gridsize <- NULL
    PoR_grid <- FALSE
  } else {
    if (!PoR_grid) {
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

  grid_PoR <- if (!PoR_grid) {
    sf::st_as_sf(grid) |>
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
  return(int)
}

#' @rdname PoR_stress2grid
#' @export
PoR_stress2grid_stats <- function(x, PoR, grid = NULL, PoR_grid = TRUE, lon_range = NULL, lat_range = NULL, gridsize = 2.5, ...) {
  if (!is.null(grid)) {
    lon_range <- lat_range <- gridsize <- NULL
    PoR_grid <- FALSE
  } else {
    if (!PoR_grid) {
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

  grid_PoR <- if (!PoR_grid) {
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
#' @param x \code{sf} object containing
#' \describe{
#' \item{azi}{Azimuth in degree}
#' \item{unc}{Uncertainties of azimuth in degree}
#' \item{prd}{Predicted value for azimuth}
#' }
#'
#' @param grid (optional) Point object of class \code{sf}.
#' @param lon_range,lat_range (optional) numeric vector specifying the minimum
#' and maximum longitudes and latitudes (are ignored if `"grid"` is specified).
#' @param gridsize Numeric. Target spacing of the regular grid in decimal
#' degree. Default is 2.5. (is ignored if `"grid"` is specified)
#' @param stat The measurement of dispersion to be calculated. Either
#' `"dispersion"` (default), `"nchisq"`, or `"rayleigh"` for circular
#' dispersion, normalized Chi-squared test statistic, or Rayleigh test
#' statistic.
#' @param stat_threshold numeric. Generates missing values when the kernel
#' `stat` value exceeds this threshold. Default is `Inf`.
#' @param min_data Integer. Minimum number of data per bin. Default is `3`
#' @param max_data integer. The number of nearest observations that should be
#' used for prediction, where "nearest" is defined in terms of the space of the
#' spatial locations. Default is `Inf`.
#' @param min_dist_threshold Numeric. Maximum distance (in km) of the grid point to the
#' next data point. Default is 200
#' @param dist_threshold Numeric. Distance weight to prevent overweight of data
#' nearby (`0` to `1`). Default is `0.1`
#' @param R_range Numeric value or vector specifying the (adaptive) kernel
#' half-width(s) as search radius (in km). Default is \code{seq(50, 1000, 50)}
#' @param ... optional arguments to [dist_greatcircle()]
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
#' @seealso [circular_dispersion()], [norm_chisq()], [rayleigh_test()]
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
  stat <- match.arg(stat)

  N <- md <- R <- numeric()

  # pre-allocating
  azi <- x$azi
  length_azi <- length(azi)
  colnames_x <- colnames(x)
  unc <- lat <- lon <- prd <- numeric(length_azi)
  type <- character(9)

  # num_r <- length(R_range)

  x_coords <-
    sf::st_coordinates(x) |>
    as.data.frame()

  datas <- cbind(
    lon = x_coords$X,
    lat = x_coords$Y,
    azi = x$azi,
    unc = x$unc,
    prd = x$prd
  )

  if (is.null(grid)) {
    # Regular grid
    if (is.null(lon_range) || is.null(lat_range)) {
      lon_range <- range(datas[, "lon"], na.rm = TRUE)
      lat_range <- range(datas[, "lat"], na.rm = TRUE)
    }

    grid <- sf::st_bbox(
      c(
        xmin = lon_range[1],
        xmax = lon_range[2],
        ymin = lat_range[1],
        ymax = lat_range[2]
      ),
      crs = sf::st_crs("WGS84")
    ) |>
      sf::st_make_grid(
        cellsize = gridsize,
        what = "centers",
        offset = c(lon_range[1], lat_range[1])
      ) |>
      sf::st_as_sf()
  }
  stopifnot(inherits(grid, "sf"), any(sf::st_is(grid, "POINT")))
  G <- unname(sf::st_coordinates(grid))

  # R <- N <- numeric(nrow(G))
  R_seq <- seq_along(R_range)
  # nR <- length(R_seq)

  # SH <- matrix(nrow = 0, ncol = 6, dimnames = list(NULL, c("lon", "lat", "stat", "R", "md", "N")))

  lapply(seq_along(G[, 1]), function(i) {
    # for (i in seq_along(G[, 1])) {
    distij <- dist_greatcircle(G[i, 2], G[i, 1], datas[, "lat"], datas[, "lon"], ...)
    if (max_data < Inf) distij <- distij[which.nsmallest(distij, max_data)] # select the `max_data` nearest locations

    if (min(distij) <= min_dist_threshold) {
      t(vapply(R_seq, function(k) {
        R_search <- R_range[k]
        ids_R <- which(distij <= R_search)
        N_in_R <- length(ids_R)

        if (N_in_R < min_data) {
          # not enough data within search radius
          y <- md <- NA
        } else if (N_in_R == 1) {
          y <- NA
          md <- distij[ids_R]
        } else {
          md <- mean(distij[ids_R], na.rm = TRUE)

          if (stat == "nchisq") {
            y <- norm_chisq(datas[ids_R, "azi"], prd = datas[ids_R, "prd"], datas[ids_R, "unc"])
          } else if (stat == "rayleigh") {
            y <- weighted_rayleigh(datas[ids_R, "azi"], mu = datas[ids_R, "prd"], w = 1 / datas[ids_R, "unc"], ...)$statistic
          } else {
            y <- circular_dispersion(datas[ids_R, "azi"], y = datas[ids_R, "prd"], w = 1 / datas[ids_R, "unc"], ...)
          }
        }

        c(
          lon = G[i, 1],
          lat = G[i, 2],
          stat = y,
          R = R_search,
          md = md,
          N = N_in_R
        )
      }, FUN.VALUE = numeric(6)))
    }
  }) |>
    lapply(as.data.frame) |>
    dplyr::bind_rows() |>
    dplyr::mutate(
      N = as.integer(N),
      mdr = md / R,
      stat = ifelse(stat >= stat_threshold, NA, stat)
    ) |>
    dplyr::select(-md) |>
    sf::st_as_sf(coords = c("lon", "lat"), crs = sf::st_crs(x), remove = FALSE)
}


#' @rdname kernel_dispersion
#' @export
#' @keywords internal
dispersion_grid <- function(...) {
  lifecycle::deprecate_warn(" 0.2.97", "dispersion_grid()", "kernel_dispersion")
  kernel_dispersion(...)
}
