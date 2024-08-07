#' Earth's radius in km
#'
#' IERS mean radius of Earth in km (based on WGS 84)
#'
#' @returns numeric value
#'
#' @export
earth_radius <- function() 6371.0087714

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
    rad2deg(c(mean_s, sd_s)) %% 180
  } else {
    c(NA, NA)
  }
}

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
  c(median_s, iqr_s)
}


dist_weight_linear <- function(R_search, dist_threshold, distij, idp = 0) {
  dist_threshold_scal <- R_search * dist_threshold
  R_search + 1 - max(dist_threshold_scal, distij)
}

dist_weight_inverse <- function(R_search, dist_threshold, distij, idp = 0) {
  dist_threshold_scal <- R_search * dist_threshold
  1 / (max(dist_threshold_scal, distij))^idp
}


#' Spatial interpolation of SHmax
#'
#' Stress field interpolation and wavelength analysis using a kernel (weighted)
#' mean/median and standard deviation/IQR of stress data
#'
#' @param x \code{sf} object containing
#' \describe{
#' \item{azi}{SHmax in degree}
#' \item{unc}{(optional) Uncertainties of SHmax in degree}
#' \item{type}{(optional) Methods used for the determination of the direction
#' of SHmax}
#' }
#' @param grid (optional) Point object of class \code{sf}.
#' @param lon_range,lat_range (optional) numeric vector specifying the minimum
#' and maximum longitudes and latitudes (ignored if `grid` is specified).
#' @param gridsize numeric. Target spacing of the regular grid in decimal
#' degree. Default is `2.5`. (is ignored if `grid` is specified)
#' @param stat whether the direction of interpolated SHmax is based on the
#' circular mean and standard deviation (\code{"mean"}, the default) or the
#' circular median and interquartile range (\code{"median"})
#' @param min_data integer. Minimum number of data per bin. Default is `3`
#' @param threshold numeric. Threshold for deviation of direction. Default is
#' 25
#' @param arte_thres numeric. Maximum distance (in km) of the grid point to the
#' next data point. Default is `200`
#' @param dist_weight Distance weighting method which should be used. One of
#' `"none"`, `"linear"`, or `"inverse"` (the default).
#' @param idp,qp,mp numeric. The weighting power of inverse distance, quality
#' and method. Default is `1`. The higher the value, the more weight it will
#' put. When set to `0`, no weighting is applied. `idp` is only effective if
#' inverse distance weighting (`dist_weight="inverse"`) is applied.
#' @param dist_threshold numeric. Distance weight to prevent overweight of data
#' nearby (0 to 1). Default is `0.1`
#' @param method_weighting logical. If a method weighting should be applied:
#' Default is \code{FALSE}. If `FALSE`, overwrites `mp`.
#' @param quality_weighting logical. If a quality weighting should be applied:
#' Default is \code{TRUE}. If `FALSE`, overwrites `qp`.
#' @param R_range numeric value or vector specifying the kernel half-width(s),
#' i.e. the search radius (in km). Default is \code{seq(50, 1000, 50)}
#' @param ... (optional) arguments to [dist_greatcircle()]
#'
#' @importFrom sf st_coordinates st_bbox st_make_grid st_crs st_as_sf
#' @importFrom dplyr group_by mutate filter rename mutate as_tibble
#'
#' @returns
#' \code{sf} object containing
#' \describe{
#' \item{lon,lat}{longitude and latitude in degrees}
#' \item{azi}{Mean SHmax in degree}
#' \item{sd}{Standard deviation of SHmax in degrees}
#' \item{R}{Search radius in km}
#' \item{mdr}{Mean distance of datapoints per search radius}
#' \item{N}{Number of data points in search radius}
#' }
#'
#' @details This is a modified version of the MATLAB script "stress2grid"
#'
#' @seealso [dist_greatcircle()], [PoR_stress2grid()], [compact_grid()],
#' [circular_mean()], [circular_median()], [circular_sd()]
#'
#' @source \url{https://github.com/MorZieg/Stress2Grid}
#'
#' @references Ziegler, M. and Heidbach, O. (2019).
#' Matlab Script Stress2Grid v1.1. GFZ Data Services. \doi{10.5880/wsm.2019.002}
#'
#' @export
#'
#' @examples
#' data("san_andreas")
#' stress2grid(san_andreas, stat = "median")
stress2grid <- function(x,
                        stat = c("mean", "median"),
                        grid = NULL,
                        lon_range = NULL,
                        lat_range = NULL,
                        gridsize = 2,
                        min_data = 3L,
                        threshold = 25,
                        arte_thres = 200,
                        method_weighting = FALSE,
                        quality_weighting = TRUE,
                        dist_weight = c("inverse", "linear", "none"),
                        idp = 1,
                        qp = 1,
                        mp = 1,
                        dist_threshold = 0.1,
                        R_range = seq(50, 1000, 50),
                        ...) {
  stopifnot(
    inherits(x, "sf"), is.numeric(gridsize), is.numeric(threshold), is.numeric(arte_thres),
    arte_thres > 0, is.numeric(dist_threshold), is.numeric(R_range), is.logical(method_weighting),
    is.logical(quality_weighting), is.numeric(idp), is.numeric(qp), is.numeric(mp)
  )

  min_data <- as.integer(ceiling(min_data))

  dist_weight <- match.arg(dist_weight)
  if (dist_weight == "linear") {
    w_distance_fun <- dist_weight_linear
  } else {
    w_distance_fun <- dist_weight_inverse
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
  num_r <- length(R_range)


  if (!quality_weighting) qp <- 0
  if (!method_weighting) mp <- 0
  if (dist_weight == "none") idp <- 0

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
  num_G <- nrow(G)

  R <- N <- numeric(num_G)

  # SH <- matrix(nrow = num_G * length(R_range), ncol = 7, dimnames = list(NULL, c('lon', 'lat', 'azi', 'sd', 'R', 'mdr', 'N')))
  # SH[, 1] <- rep(G[, 1], length(R_range))
  # SH[, 2] <- rep(G[, 2], length(R_range))
  SH <- matrix(nrow = 0, ncol = 7, dimnames = list(NULL, c("lon", "lat", "azi", "sd", "R", "md", "N")))

  for (i in seq_along(G[, 1])) {
    distij <- dist_greatcircle(G[i, 2], G[i, 1], datas[, 2], datas[, 1], ...)

    if (min(distij) <= arte_thres) {
      for (k in seq_along(R_range)) {
        R_search <- R_range[k]
        # ids_R <- which(distij <= R_search) # select those that are in search radius
        # N_in_R <- length(ids_R)
        ids_R <- (distij <= R_search) # select those that are in search radius
        N_in_R <- sum(ids_R)

        if (N_in_R < min_data) {
          # not enough data within search radius
          sd <- 0
          meanSH <- md <- NA
        } else if (N_in_R == 1) {
          sd <- 0
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
          sd <- stats[2]
        }
        SH.ik <- c(
          G[i, 1], # lon
          G[i, 2], # lat
          meanSH, # azi
          sd, # sd
          R_search, # R_search
          md, # mdr
          N_in_R # N_in_R
        )

        SH <- rbind(SH, SH.ik)
      }
    }
  }

  res <- dplyr::as_tibble(SH) |>
    dplyr::mutate(N = as.integer(N), sd = sd / 2, mdr = md / R) |>
    dplyr::select(-md) |>
    dplyr::filter(!is.na(azi), sd <= threshold, !is.na(sd)) |>
    sf::st_as_sf(coords = c("lon", "lat"), crs = sf::st_crs(x), remove = FALSE)

  return(res)
}



#' Spatial interpolation of SHmax in PoR coordinate reference system
#'
#' Stress field and wavelength analysis in PoR system and back-transformed
#'
#' @param x \code{sf} object containing
#' \describe{
#' \item{azi}{SHmax in degree}
#' \item{unc}{Uncertainties of SHmax in degree}
#' \item{type}{Methods used for the determination of the orientation of SHmax}
#' }
#' @param PoR Pole of Rotation. \code{"data.frame"} or object of class \code{"euler.pole"}
#' containing the geographical coordinates of the Euler  pole
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
#' @importFrom dplyr rename as_tibble group_by
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
#' @export
#'
#' @examples
#' data("san_andreas")
#' data("nuvel1")
#' PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
#' PoR_stress2grid(san_andreas, PoR)
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
    dplyr::as_tibble() |>
    dplyr::rename(lat = Y, lon = X)

  azi <- lat <- lon <- lat.PoR <- lon.PoR <- X <- Y <- R <- numeric() # pre allocating:

  x_PoR$lat <- x_PoR_coords$lat
  x_PoR$lon <- x_PoR_coords$lon
  x_PoR$azi <- PoR_shmax(x, PoR)

  int <- stress2grid(x_PoR, grid = grid_PoR, lon_range = lon_range, lat_range = lat_range, gridsize = gridsize, ...) |>
    dplyr::rename(azi.PoR = azi, lat.PoR = lat, lon.PoR = lon) |>
    PoR_to_geographical_sf(PoR)
  int_coords <- sf::st_coordinates(int) |>
    dplyr::as_tibble() |>
    dplyr::rename(lat = Y, lon = X)
  int$lat <- int_coords$lat
  int$lon <- int_coords$lon
  int$azi <- PoR2Geo_azimuth(int, PoR)
  return(int)
}

#' Compact smoothed stress field
#'
#' Filter smoothed stress field containing a range of search radii or kernel
#' half widths to find smallest wavelength (R) with the least circular sd. or
#' dispersion for each coordinate, respectively.
#'
#' @param x output of [stress2grid()], [PoR_stress2grid()], or [kernel_dispersion()]
#' @param type character. Type of the grid `x`. Either `"stress"` (when input
#' is [stress2grid()] or [PoR_stress2grid()]) or `"dispersion"` (when input
#' is [kernel_dispersion()]).
#'
#' @returns \code{sf} object
#'
#' @importFrom dplyr ungroup mutate select left_join as_tibble
#' @importFrom tidyr drop_na
#' @importFrom sf st_as_sf
#' @importFrom stats aggregate
#' @seealso [stress2grid()], [PoR_stress2grid()], [kernel_dispersion()]
#'
#' @export
#'
#' @examples
#' data("san_andreas")
#' res <- stress2grid(san_andreas)
#' compact_grid(res)
compact_grid <- function(x, type = c("stress", "dispersion")) {
  lon <- lat <- azi <- R <- stat <- numeric()
  group <- character()
  type <- match.arg(type)

  if (type == "stress") {
    data <- x |>
      # dplyr::ungroup() |>
      dplyr::as_tibble() |>
      tidyr::drop_na(azi) |>
      dplyr::mutate(group = paste(lon, lat))
  } else {
    data <- x |>
      # dplyr::ungroup() |>
      dplyr::as_tibble() |>
      tidyr::drop_na(stat) |>
      dplyr::mutate(group = paste(lon, lat))
  }

  aggregate(R ~ group, data, min, na.rm = TRUE) |>
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
#' `"dispersion"` (default), `"nchisq"`, or `"rayleigh"` for circular dispersion,
#' normalized Chi-squared test statistic, or Rayleigh test statistic.
#' @param min_data Integer. Minimum number of data per bin. Default is 3
#' @param threshold Numeric. Threshold for stat value (default is 1)
#' @param arte_thres Numeric. Maximum distance (in km) of the grid point to the
#' next data point. Default is 200
#' @param dist_threshold Numeric. Distance weight to prevent overweight of data
#' nearby
#' (0 to 1). Default is 0.1
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
#' @export
#'
#' @examples
#' data("nuvel1")
#' PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
#' san_andreas_por <- san_andreas
#' san_andreas_por$azi <- PoR_shmax(san_andreas, PoR, "right")$azi.PoR
#' san_andreas_por$prd <- 135
#' kernel_dispersion(san_andreas_por)
kernel_dispersion <- function(x,
                              stat = c("dispersion", "nchisq", "rayleigh"),
                              grid = NULL,
                              lon_range = NULL,
                              lat_range = NULL,
                              gridsize = 2.5,
                              min_data = 3,
                              threshold = 1,
                              arte_thres = 200,
                              dist_threshold = 0.1,
                              R_range = seq(100, 2000, 100),
                              ...) {
  stopifnot(
    inherits(x, "sf"), is.numeric(gridsize), is.numeric(threshold), is.numeric(arte_thres),
    arte_thres > 0, is.numeric(dist_threshold), is.numeric(R_range)
  )
  stat <- match.arg(stat)
  min_data <- as.integer(ceiling(min_data))
  stat <- match.arg(stat)

  # pre-allocating
  azi <- x$azi
  length_azi <- length(azi)
  colnames_x <- colnames(x)
  unc <- lat <- lon <- prd <- numeric(length_azi)
  type <- character(9)

  num_r <- length(R_range)

  x_coords <-
    sf::st_coordinates(x) |>
    as.data.frame()

  datas <- data.frame(
    lon = x_coords$X,
    lat = x_coords$Y,
    azi = x$azi,
    unc = x$unc,
    prd = x$prd
  )

  if (is.null(grid)) {
    # Regular grid
    if (is.null(lon_range) || is.null(lat_range)) {
      lon_range <- range(datas$lon, na.rm = TRUE)
      lat_range <- range(datas$lat, na.rm = TRUE)
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
  G <- grid |>
    sf::st_coordinates()

  R <- N <- numeric(nrow(G))


  SH <- matrix(nrow = 0, ncol = 6, dimnames = list(NULL, c("lon", "lat", "stat", "R", "md", "N")))

  for (i in seq_along(G[, 1])) {
    distij <- dist_greatcircle(G[i, 2], G[i, 1], datas$lat, datas$lon, ...)

    if (min(distij) <= arte_thres) {
      for (k in seq_along(R_range)) {
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
          # dist_threshold_scal <- R_search * dist_threshold

          if (stat == "nchisq") {
            y <- norm_chisq(datas$azi[ids_R], prd = datas$prd[ids_R], datas$unc[ids_R])
          } else if (stat == "rayleigh") {
            y <- weighted_rayleigh(datas$azi[ids_R], mu = datas$prd[ids_R], w = 1 / datas$unc[ids_R], ...)$statistic
          } else {
            y <- circular_dispersion(datas$azi[ids_R], y = datas$prd[ids_R], w = 1 / datas$unc[ids_R], ...)
          }
        }

        SH.ik <- c(
          G[i, 1],
          G[i, 2],
          y,
          R_search,
          md,
          N_in_R
        )

        SH <- rbind(SH, SH.ik)
      }
    }
  }

  res <- dplyr::as_tibble(SH) |>
    dplyr::mutate(N = as.integer(N), mdr = md / R) |>
    dplyr::select(-md) |>
    sf::st_as_sf(coords = c("lon", "lat"), crs = sf::st_crs(x), remove = FALSE)

  return(res)
}

dispersion_grid <- function(...) {
  .Deprecated("kernel_dispersion")
  kernel_dispersion(...)
}
