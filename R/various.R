#' World Stress Map Database (WSM)
#'
#' Download WSM2025 or WSM2016 database from the GFZ sever and applies optional filters.
#' If `destdir` is specified, the data can be reloaded in a later R session
#' using `load_WSM()` using the same arguments.
#'
#' @param quality a character vectors containing the quality levels to be
#' included. Includes all quality ranks (A-X) by default.
#' @param lat_range,lon_range two-element numeric vectors giving the range of
#' latitudes and longitudes (in degrees).
#' @param depth_range two-element numeric vectors giving the depth interval
#' (in km)
#' @param type a character vectors containing the methods of stress inversion
#' to be included. Includes all methods by default. See WSM2016 manual for used
#' acronyms.
#' @param regime a character vectors containing the stress regimes to be
#' included. Acronyms: `"N"` - normal, `"T"` - thrust, `"S"` - strike-slip,
#' `"NS"` - oblique normal, `"TS"` - oblique thrust, and `NA` - unknown faulting
#' @param load `TRUE` load the dataset into R, `FALSE` return the
#' file name of the downloaded object.
#' @param version character. Version of the World stress map database. Either
#' `"2025"` (default) or `"2016"`
#' @param file 	the name of the file which the data are to be read from.
#' @param destdir where to save files, defaults to [base::tempdir()],
#' [base::getwd()] is also possible.
#' @param ... (optional) arguments passed to [load_WSM()]
#'
#' @returns `sf` object of and the parsed numeric uncertainty (`unc`) based on
#' the reported standard deviation and the quality rank. If `load=FALSE`,
#' the path to the downloaded file is returned.
#'
#' @note Because of R-compatibility and easy readability reasons, the downloaded
#' dataset is a modified version of the original, WSM server version:
#' All column names have been changed from uppercase (in the original dataset) to
#' lowercase characters.
#' Unknown azimuth values are represented by `NA` values instead of `999` in
#' the original.
#' Unknown regimes are represented by `NA` instead of "U" in the original.
#'
#' @source \url{https://datapub.gfz.de/download/10.5880.WSM.2025.001-Scbwez/WSM_Database_2025.csv}
#'
#' \url{https://datapub.gfz-potsdam.de/download/10.5880.WSM.2016.001/wsm2016.csv}
#'
#' @references Heidbach, O., M. Rajabi, X. Cui, K. Fuchs, B. M<U+00FC>ller, J.
#' Reinecker, K. Reiter, M. Tingay, F. Wenzel, F. Xie, M. O. Ziegler,
#' M.-L. Zoback, and M. D. Zoback (2018): The World Stress Map database
#' release 2016: Crustal stress pattern across scales. *Tectonophysics*,
#' **744**, 484-498, \doi{10.1016/j.tecto.2018.07.007}.
#'
#' Heidbach, Oliver; Rajabi, Mojtaba; Di Giacomo, Domenico; Harris, James;
#' Lammers, Steffi; Morawietz, Sophia; Pierdominici, Simona; Reiter, Karsten;
#' von Specht, Sebastian; Storchak, Dmitry; Ziegler, Moritz O. (2025): World
#' Stress Map Database Release 2025. GFZ Data Services.
#' \doi{10.5880/WSM.2025.001}
#'
#' @importFrom dplyr rename_with filter between mutate select case_when
#' @importFrom sf st_as_sf st_crs
#' @importFrom utils read.csv download.file
#'
#' @keywords datasets
#'
#' @name import_WSM
#' @examples
#' \dontrun{
#' download_WSM(
#'   quality = c("A", "B", "C"), lat_range = c(51, 72),
#'   lon_range = c(-180, -130), depth_range = c(0, 10), type = "FMS"
#' )
#' }
NULL

#' @rdname import_WSM
#' @export
download_WSM <- function(destdir = tempdir(), load = TRUE, version = c("2025", "2016"), ...) {
  stopifnot(is.character(destdir))
  version <- match.arg(version)

  options(timeout = 999)

  if (version == "2025") {
    address <- "https://datapub.gfz.de/download/10.5880.WSM.2025.001-Scbwez/WSM_Database_2025.csv"
  } else {
    address <- "https://datapub.gfz-potsdam.de/download/10.5880.WSM.2016.001/wsm2016.csv"
  }
  dat_file <- tempfile("wsm", destdir, fileext = ".csv")

  tryCatch(
    utils::download.file(address, destfile = dat_file)
  )

  options(timeout = 60)

  if (load) {
    load_WSM(dat_file, ...)
  } else {
    return(dat_file)
  }
}

#' @rdname import_WSM
#' @export
load_WSM <- function(file, quality = c("A", "B", "C", "D", "E", "X"),
                     lat_range = c(-90, 90), lon_range = c(-180, 180),
                     depth_range = c(-Inf, Inf),
                     type = c("BO", "BOC", "BOT", "BS", "DIF", "FMA", "FMF", "FMS", "GFI", "GFM", "GFS", "GVA", "HF", "HFG", "HFM", "HFH", "HFP", "HFS", "OC", "PC", "SWB", "SWL", "SWS"),
                     regime = c("N", "NS", "T", "TS", "S", NA)) {
  stopifnot(
    is.character(quality),
    is.character(file),
    is.character(regime) | is.na(regime),
    is.character(type),
    is.numeric(lon_range),
    is.numeric(lat_range),
    is.numeric(depth_range)
  )

  reg_flt <- regime
  meth_flt <- type

  lat <- lon <- azi <- quality.numeric <- unc <- id <- depth <- numeric()
  type <- locality <- iso <- type <- date <- time <- isc_id <- character()
  eq_mag <- dist <- NULL

  qlt_flt <- quality
  lat_range <- sort(lat_range)
  lon_range <- sort(lon_range)
  depth_range <- sort(depth_range)

  out <- utils::read.csv(file) |>
    dplyr::rename_with(tolower) |>
    dplyr::filter(
      !is.na(lat), !is.na(lon),
      dplyr::between(lat, lat_range[1], lat_range[2]),
      dplyr::between(lon, lon_range[1], lon_range[2]),
      dplyr::between(depth, depth_range[1], depth_range[2]),
      type %in% meth_flt
    ) |>
    dplyr::mutate(
      azi = ifelse(azi == 999, NA, azi),
      quality = gsub("^X.*", "X", quality),
      quality.numeric = parse_wsm_quality(quality),
      quality = factor(quality),
      regime = dplyr::case_when(
        regime == "NF" ~ "N",
        regime == "NS" ~ "NS",
        regime == "TF" ~ "T",
        regime == "TS" ~ "TS",
        regime == "SS" ~ "S",
        regime == "U" ~ NA,
        TRUE ~ NA_character_,
        .default = regime
      ),
      unc = ifelse(is.na(sd), quality.numeric, sd),
      unc = ifelse(unc > quality.numeric, quality.numeric, unc),
      unc = ifelse(unc == 0, 15, unc),
      unc = as.numeric(unc),
      dplyr::across(dplyr::where(is.character), function(x) {
        ifelse(x == "", NA, x)
      }),
      eq_mag = as.numeric(eq_mag),
      # date = ymd(date),
      # time = hms(time)
    ) |>
    dplyr::filter(quality %in% qlt_flt, regime %in% reg_flt) |>
    sf::st_as_sf(coords = c("lon", "lat"), crs = sf::st_crs("WGS84"), remove = FALSE)

  if ("isc_id" %in% colnames(out)) {
    dplyr::select(out, id, lat, lon, azi, unc, type, depth, quality, regime, isc_id, locality:dist)
  } else {
    dplyr::select(out, id, lat, lon, azi, unc, type, depth, quality, regime, locality:iso) # 2016
  }
}

#' @rdname import_WSM
#' @export
download_WSM2016 <- function(destdir = tempdir(), load = TRUE, ...) {
  lifecycle::deprecate_warn("0.4.7", "download_WSM2016()", "download_WSM2016()")

  stopifnot(is.character(destdir))

  options(timeout = 999)

  address <- "https://datapub.gfz-potsdam.de/download/10.5880.WSM.2016.001/wsm2016.csv"
  dat_file <- tempfile("wsm2016_", destdir, fileext = ".csv")
  tryCatch(
    utils::download.file(address, destfile = dat_file)
  )

  options(timeout = 60)

  if (load) {
    load_WSM2016(dat_file, ...)
  } else {
    return(dat_file)
  }
}

#' @rdname import_WSM
#' @export
load_WSM2016 <- function(file, quality = c("A", "B", "C", "D", "E"), lat_range = c(-90, 90), lon_range = c(-180, 180),
                         depth_range = c(-Inf, Inf),
                         type = c("BO", "BOC", "BOT", "BS", "DIF", "FMA", "FMF", "FMS", "GFI", "GFM", "GFS", "GVA", "HF", "HFG", "HFM", "HFP", "OC", "PC", "SWB", "SWL", "SWS"),
                         regime = c("N", "NS", "T", "TS", "S", NA)) {
  lifecycle::deprecate_warn("0.4.7", "load_WSM2016()", "load_WSM()")

  stopifnot(
    is.character(quality),
    is.character(file),
    is.character(regime) | is.na(regime),
    is.character(type),
    is.numeric(lon_range),
    is.numeric(lat_range),
    is.numeric(depth_range)
  )

  reg_flt <- regime
  meth_flt <- type

  lat <- lon <- azi <- quality.numeric <- unc <- id <- depth <- numeric()
  type <- locality <- iso <- type <- date <- time <- character()
  eq_mag <- NULL

  qlt_flt <- quality
  lat_range <- sort(lat_range)
  lon_range <- sort(lon_range)
  depth_range <- sort(depth_range)

  utils::read.csv(file) |>
    dplyr::rename_with(tolower) |>
    dplyr::filter(
      dplyr::between(lat, lat_range[1], lat_range[2]),
      dplyr::between(lon, lon_range[1], lon_range[2]),
      dplyr::between(depth, depth_range[1], depth_range[2]),
      type %in% meth_flt
    ) |>
    dplyr::mutate(
      azi = ifelse(azi == 999, NA, azi),
      quality = factor(quality),
      regime = dplyr::case_when(
        regime == "NF" ~ "N",
        regime == "NS" ~ "NS",
        regime == "TF" ~ "T",
        regime == "TS" ~ "TS",
        regime == "SS" ~ "S",
        regime == "U" ~ NA,
        TRUE ~ NA_character_,
        .default = regime
      ),
      quality.numeric = parse_wsm_quality(quality),
      unc = ifelse(is.na(sd), quality.numeric, sd),
      unc = ifelse(unc > quality.numeric, quality.numeric, unc),
      unc = ifelse(unc == 0, 15, unc),
      unc = as.numeric(unc),
      dplyr::across(dplyr::where(is.character), function(x) {
        ifelse(x == "", NA, x)
      }),
      eq_mag = as.numeric(eq_mag),
      # date = ymd(date),
      # time = hms(time)
    ) |>
    dplyr::filter(quality %in% qlt_flt, regime %in% reg_flt) |>
    sf::st_as_sf(coords = c("lon", "lat"), crs = sf::st_crs("WGS84"), remove = FALSE) |>
    dplyr::select(id, lat, lon, azi, unc, type, depth, quality, regime, locality:iso)
}


#' Numerical values to World Stress Map Quality Ranking
#'
#' Assigns numeric values of the precision (sd.) of each measurement to the
#' categorical quality ranking of the World Stress Map (A, B, C, D, E, X).
#'
#' @param x Either a string or a character/factor vector of WSM quality ranking
#'
#' @returns `"numeric"`. the standard deviation of stress azimuth
#'
#' @references Heidbach, O., Barth, A., M<U+00FC>ller, B., Reinecker, J.,
#' Stephansson, O., Tingay, M., Zang, A. (2016). WSM quality
#' ranking scheme, database description and analysis guidelines for stress
#' indicator. *World Stress Map Technical Report* **16-01**, GFZ German Research
#' Centre for Geosciences. \doi{10.2312/wsm.2016.001}
#'
#' @name parse_wsm
#' @examples
#' parse_wsm_quality(c("A", "B", "C", "D", NA, "E", "X"))
#' data("san_andreas")
#' head(parse_wsm_quality(san_andreas$quality))
NULL

#' @rdname parse_wsm
#' @export
parse_wsm_quality <- function(x) {
  c(
    "A" = 15,
    "B" = 20,
    "C" = 25,
    "D" = 40,
    "E" = 90,
    "X" = 180,
    "Xmi" = 180,
    "Xne" = 180,
    "Xru" = 180
  )[x]
}

#' @rdname parse_wsm
#' @export
quantise_wsm_quality <- function(x) {
  lifecycle::deprecate_warn("0.3.7", "quantise_wsm_quality()", "parse_wsm_quality()")
  as.numeric(sapply(x, FUN = parse_wsm_quality))
}


#' Quick analysis of a stress data set
#'
#' Returns the converted azimuths, distances to the plate boundary,
#' statistics of the model, and some plots.
#'
#' @param x `data.frame` or `sf` object containing the coordinates of the point(s)
#' (`lat`, `lon`), the direction of
#' \eqn{\sigma_{Hmax}}{SHmax} `azi` and its standard deviation
#' \code{unc} (optional)
#' @param PoR Pole of Rotation. `data.frame` or object of class `"euler.pole"`
#' containing the geographical coordinates of the Euler pole
#' @param type Character. Type of plate boundary (optional). Can be
#' `"out"`, `"in"`, `"right"`, or
#' `"left"` for outward, inward, right-lateral, or left-lateral
#' moving plate boundaries, respectively. If `"none"` (the default), only
#' the PoR-equivalent azimuth is returned.
#' @param pb (optional) `sf` object of the plate boundary geometries in the
#' geographical coordinate system
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
#' \item{`stats`}{array with circular (weighted) mean, circular standard
#' deviation, circular variance, circular median, skewness, kurtosis, the 95%
#' confidence angle, circular dispersion, and the normalized Chi-squared test
#'  statistic}
#' \item{`test`}{list containing the test results of the (weighted) Rayleigh
#' test against the uniform distribution about the predicted orientation.}
#' }
#'
#' @export
#'
#' @seealso [PoR_shmax()], [distance_from_pb()], [norm_chisq()], [quick_plot()], [circular_summary()]
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
  tangential <- type %in% c("right", "left")
  res <- PoR_shmax(x, PoR, type)
  res <- cbind(res, PoR_coordinates(x, PoR))
  # if (!missing(pb)) {
  res$distance <- distance_from_pb(x, PoR, pb, tangential, ...)
  # }
  prd <- res$prd

  stats <- circular_summary(res$azi.PoR, 1 / x$unc)

  median <- circular_sample_median(res$azi.PoR)
  # mean <- circular_mean(res$azi.PoR, 1 / x$unc)
  # sd <- circular_sd(res$azi.PoR, 1 / x$unc)
  # var <- circular_var(res$azi.PoR, 1 / x$unc)
  disp <- circular_dispersion(res$azi.PoR, prd, 1 / x$unc)
  # conf <- confidence_angle(res$azi.PoR, w = 1 / x$unc)
  nchisq <- norm_chisq(res$azi.PoR, prd, unc = x$unc)
  rayleigh <- weighted_rayleigh(res$azi.PoR, prd, w = 1 / x$unc, quiet = TRUE)

  if (plot) {
    PoR_map(x, PoR, pb, type = type, deviation = TRUE)
    grDevices::dev.new()
    quick_plot(azi = res$azi.PoR, distance = res$distance, unc = x$unc, regime = x$regime, prd = prd)
  }

  list(
    result = res,
    stats =
      rbind(
        mean = stats["mean"], sd = stats["sd"], var = stats["var"],
        median = stats["median"],
        quasi_median = stats["quasi-median"], skewness = stats["skewness"], kurtosis = stats["kurtosis"],
        conf95 = stats["CI"],
        dispersion = disp, norm_chisq = nchisq
      ),
    test = rayleigh
  )
}

#' Extract azimuths of line segments
#'
#' @param x sf object of type `"LINESTRING"` or `"MULTILINESTRING"`
#' @param warn logical; if `TRUE`, warn if `"MULTILINESTRING"` (default).
#'
#' @return sf object of type `"POINT"` with the columns and entries of the first row of `x`
#'
#' @details
#' It is recommended to perform `line_azimuth()` on single line objects, i.e.
#' type `"LINESTRING"`, instead of `"MULTILINESTRING"`. This is because the azimuth
#' of the last point of a line will be calculated to the first point of the
#' next line otherwise. This will cause a warning message (if `warn = TRUE`).
#' For `"MULTILINESTRING"` objects, use `lines_azimuths()`.
#'
#' @importFrom sf st_cast st_coordinates st_as_sf st_crs st_drop_geometry
#' @export
#'
#' @name line_azimuth
#' @examples
#' data("plates")
#'
#' # one line:
#' subset(plates, pair == "af-eu") |>
#'   smoothr::densify() |>
#'   line_azimuth()
#'
#' # multiple lines:
#' lines_azimuths(plates[1:5, ])
NULL

#' @rdname line_azimuth
#' @export
line_azimuth <- function(x, warn = TRUE) {
  if (warn & (nrow(x) > 1 | any(sf::st_geometry_type(x) == "MULTILINESTRING"))) {
    warning("MULTILINESTRING object is not recommended")
  }
  mat <- sf::st_cast(x, "POINT", warn = FALSE) |>
    sf::st_coordinates()

  n <- nrow(mat)
  # a <- numeric(n - 1)
  # for (i in 1:(n - 1)) {
  #   a[i] <- get_azimuth(mat[i, 2], mat[i, 1], mat[i + 1, 2], mat[i + 1, 1])
  # }
  a <- vapply(
    1:(n - 1), function(i) {
      get_azimuth(mat[i, 2], mat[i, 1], mat[i + 1, 2], mat[i + 1, 1])
    },
    numeric(1)
  )

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
  la <- split(x, 1:nrow(x)) |>
    lapply(line_azimuth, warn = FALSE)
  do.call(rbind, la)
}
