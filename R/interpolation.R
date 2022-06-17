wcmean <- function(x, w) {
  Z <- sum(w)
  meansin2 <- sum(w * sind(2 * x)) / Z
  meancos2 <- sum(w * cosd(2 * x)) / Z
  meanR <- sqrt(meansin2^2 + meancos2^2)

  if (meanR > 1) {
    sd_s <- 0
  } else {
    sd_s <- sqrt(-2 * log(meanR)) / 2
  }

  mean_s <- atan2d(meansin2, meancos2) / 2
  rad2deg(c(mean_s, sd_s))
}

wcmedian <- function(x, w) {
  if (length(x) > 3) {
    quantiles <- circular_weighted_quantiles(x, w)
    median_s <- quantiles[3]
    iqr_s <- deviation_norm(as.numeric(quantiles[4] - quantiles[2]))
  } else {
    median_s <- circular_weighted_median(x, w)
    iqr_s <- ceiling(
      abs(
        deviation_norm(max(x)) - deviation_norm(min(x))
        ) / 2
      )
  }
  c(median_s, iqr_s)
}

#' Stress2Grid
#'
#' Stress pattern and wavelength analysis
#'
#' @param x \code{sf} object containing
#' \describe{
#' \item{azi}{SHmax in degree}
#' \item{unc}{Uncertainties of SHmax in degree}
#' \item{type}{Methods used for the determination of the orientation of SHmax}
#' }
#' @param stat Should the orientation of interpolated SHmax be based  on the
#' circular mean and standard deviation (\code{"mean"}, the default) or on the circular median
#' and interquartile range (\code{"median"})?
#' @param lon_range,lat_range (optional) numeric vector specifying the minimum
#' and maximum longitudes and latitudes.
#' @param gridsize Numeric. Spacing of the regular grid in decimal degree.
#' @param min_data Integer. Minimum number of data per bin. Default is 3
#' @param threshold Numeric. Threshold for deviation of orientation. Default is 25
#' @param arte_thres Numeric. Maximum distance (in km) of the gridpoint to the next
#' datapoint. Default is 200
#' @param method_weighting Logical. If a method weighting should be applied:
#' Default is \code{FALSE}.
#' @param quality_weighting Logical. If a quality weighting should be applied: Default is
#' \code{TRUE}.
#' @param dist_weight Distance weighting method which should be used. One of
#' "none", "linear", or "inverse" (the default).
#' @param dist_threshold Numeric. Distance weight to prevent overweight of data nearby
#' (0 to 1). Default is 0.1
#' @param R_range Numeric value or vector specifying the search radius (im km).
#' @importFrom sf st_coordinates st_bbox st_make_grid st_crs st_distance st_as_sf
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by mutate
#' @returns
#' \code{sf} object containing
#' \describe{
#' \item{lon,lat}{longitude and latitude in degree}
#' \item{azi}{Mean SHmax in degree}
#' \item{sd}{Standard deviation of SHmax in degree}
#' \item{R}{Search radius in km}
#' \item{mdr}{Mean distance of datapoint per search radius}
#' \item{N}{Number of data points in search radius}
#' }
#' @details Calculates the weighted mean and standard deviation of the stress
#' data
#' @source \url{https://github.com/MorZieg/Stress2Grid}
#' @references Ziegler, M. O. and Heidbach, O. (2019).
#' Matlab Script Stress2Grid v1.1. GFZ Data Services. \doi{10.5880/wsm.2019.002}
#' @export
#' @examples
#' data("san_andreas")
#' stress2grid(san_andreas)
stress2grid <- function(x,
                        stat = c("mean", "median"),
                        lon_range = NULL,
                        lat_range = NULL,
                        gridsize = 2.5,
                        min_data = 3,
                        threshold = 25,
                        arte_thres = 200,
                        method_weighting = FALSE,
                        quality_weighting = TRUE,
                        dist_weight = c("inverse", "linear", "none"),
                        dist_threshold = 0.1,
                        R_range = seq(50, 1000, 50)) {
  stopifnot(inherits(x, "sf"))
  stopifnot(is.numeric(gridsize))
  stopifnot(is.numeric(threshold))
  stopifnot(is.numeric(arte_thres) & arte_thres > 0)
  stopifnot(is.numeric(dist_threshold))
  stopifnot(is.numeric(R_range))
  stopifnot(is.logical(method_weighting))
  stopifnot(is.logical(quality_weighting))

  min_data <- as.integer(ceiling(min_data))

  dist_weight <- match.arg(dist_weight)
  stat <- match.arg(stat)
  azi <- unc <- type <- lat <- lon <- R <- NULL

  num_r <- length(R_range)

  # WSM method weighting (from 0 to 5)
  if (method_weighting & "type" %in% colnames(x)) {
    method_weights <- data.frame(
      type = c("FMS", "FMF", "BO", "DIF", "HF", "GF", "GV", "OC", NA),
      w_method = c(4, 5, 5, 5, 4, 5, 4, 2, 1) / 5
    )
    x <- dplyr::left_join(x, method_weights)
  } else {
    x$w_method <- rep(1, length(x$azi))
  }

  if (quality_weighting & "unc" %in% colnames(x)) {
    x$w_quality <- 1 / x$unc
  } else {
    x$w_quality <- rep(1, length(x$azi))
  }

  x_coords <- sf::st_coordinates(x) %>% as.data.frame()

  datas <- data.frame(
    lon = x_coords$X,
    lat = x_coords$Y,
    x = x_coords$X,
    y = x_coords$Y,
    azi = x$azi,
    # sin2 = sind(2 * x$azi),
    # cos2 = cosd(2 * x$azi),
    w_method = x$w_method,
    w_quality = x$w_quality # ,
    # pb_dist = 1 # rep(1, length(lat))
  ) %>%
    mutate(
      w_method = ifelse(is.na(w_method), 1/5, w_method)
    ) %>%
    sf::st_as_sf(coords = c("x", "y"), crs = sf::st_crs(x))

  if(quality_weighting){
    datas <- filter(datas, !is.na(w_quality))
  }

  # Regular grid
  if (missing(lon_range) | missing(lat_range)) {
    lon_range <- range(datas$lon, na.rm = TRUE)
    lat_range <- range(datas$lat, na.rm = TRUE)
  }

  G <-
    sf::st_bbox(
      c(
        xmin = lon_range[1],
        xmax = lon_range[2],
        ymin = lat_range[1],
        ymax = lat_range[2]
      ),
      crs = sf::st_crs(x)
    ) %>%
    sf::st_make_grid(cellsize = gridsize, what = "centers")

  G_coords <- sf::st_coordinates(G) %>% as.data.frame()
  XG <- G_coords$X
  YG <- G_coords$Y
  n_G <- length(XG)

  SH <- c()

  for (i in 1:n_G) {
    distij <-
      sf::st_distance(datas, G[i], which = "Great Circle") %>% # in meter
      as.numeric() / 1000 # in km

    if (min(distij) <= arte_thres) {
      for (k in 1:length(R_range)) {
        R_search <- R_range[k]
        ids_R <-
          which(distij <= R_search) # select those that are in search radius

        N_in_R <- length(ids_R)
        # sq_R  <-  sum(datas$pb_dist[ids_R])

        if (N_in_R < min_data) {
          # not enough data within search radius
          sd <- 0
          meanSH <- NA
          mdr <- NA
        } else if (N_in_R == 1) {
          sd <- 0
          meanSH <- datas$azi[ids_R]
          mdr <- distij[ids_R] / R_search
        } else {
          mdr <- mean(distij[ids_R], na.rm = TRUE) / R_search

          dist_threshold_scal <- R_search * dist_threshold

          if (dist_weight == "inverse") {
            w_distance <- 1 / max(dist_threshold_scal, distij[ids_R])
          } else if (dist_weight == "linear") {
            w_distance <- R_search + 1 - max(dist_threshold_scal, distij[ids_R])
          } else {
            w_distance <- rep(1, length(ids_R))
          }

          w <- w_distance * datas$w_quality[ids_R] * datas$w_method[ids_R]

          # mean value
          if (stat == "median") {
            stats <- wcmedian(datas$azi[ids_R], w)
          } else {
            stats <- wcmean(datas$azi[ids_R], w)
          }
          meanSH <- stats[1] %% 180
          sd <- stats[2]
          #   sw_R <- sum(w_R)
          #
          #
          #   array_sin2 <- w_R * datas$sin2[ids_R]
          #   array_cos2 <- w_R * datas$cos2[ids_R]
          #
          #   sumsin2 <- sum(array_sin2)
          #   sumcos2 <- sum(array_cos2)
          #   meansin2 <- sumsin2 / sw_R
          #   meancos2 <- sumcos2 / sw_R
          #   meanR <- sqrt(meansin2^2 + meancos2^2)
          #
          #   if (meanR > 1) {
          #     sd <- -Inf
          #   } else {
          #     sd <- rad2deg(sqrt(-2 * log(meanR)) / 2)
          #   }
          #
          #   meanSH <- (atan2d(meansin2, meancos2) / 2) %% 180
          #   #meanSH <- rad2deg(meanSH_rad)
          # }
        }
        SH.ik <- c(
          lon = XG[i],
          lat = YG[i],
          azi = meanSH,
          sd = sd,
          R = R_search,
          mdr = mdr,
          N = as.integer(N_in_R)
        )

        if (SH.ik[4] <= threshold) {
          SH <- rbind(SH, SH.ik)
        }
      }
    }
  }

  res <- as.data.frame(SH) %>%
    mutate(x = lon, y = lat) %>%
    sf::st_as_sf(coords = c("x", "y"), crs = sf::st_crs(x)) %>%
    group_by(R)

  return(res)
}