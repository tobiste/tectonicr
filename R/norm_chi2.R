#' Normalized chi-square test
#'
#' A quantitative comparison between the predicted and observed directions of
#' SHmax is obtained by the calculation of the average azimuth and by a
#' normalized chi-square test.
#'
#' @author 2021 Tobias Stephan
#' @references Wdowinski, S., 1998, A theory of intraplate
#'   tectonics: Journal of Geophysical Research: Solid Earth, v. 103, p.
#'   5037-5059, http://dx.doi.org/10.1029/97JB03390.
#' @param obs observed SHmax, numeric vector
#' @param prd predicted SHmax, numeric vector of length of \code{obs}
#' @param unc uncertainty of observed SHmax, either numeric vector of length of
#' \code{obs} or a number
#' @return numeric vector
#' @details Test result are values are between 0 and 1 indicating the quality of
#' the predicted SHmax directions. Low values (<= 0.15) indicate good agreement,
#' high values (>0.7) indicate a systematic misfit between predicted and
#' observed SHmax directions
#' @export
#' @examples
#' data("nuvel1")
#' euler <- subset(nuvel1, nuvel1$ID == "na") # North America relative to
#' # Pacific plate
#' point <- data.frame(lat = 45, lon = 20)
#' prd <- model_shmax(point, euler)
#' norm_chi2(obs = 90, prd$sc, unc = 10)
norm_chi2 <- function(obs, prd, unc) {
  if (length(prd) != length(obs)) {
    stop("Observed and predicted values must have the same length!")
  }

  if (length(unc) != 1 & length(unc) != length(obs)) {
    stop("Uncertainties must be either a numeric value or a numeric value with length of observed values!")
  }

  if (length(unc) == 1) {
    unc <- rep(unc, length(obs))
  }

  w <- c()
  x <- c()
  y <- c()
  for (i in seq_along(obs)) {
    if (is.na(prd[i]) | is.na(obs[i])) {
      x[i] <- NA
      y[i] <- NA
    } else {
      w[i] <- deviation_norm(prd[i] - obs[i]%%180)
      x[i] <- (w[i] / unc[i])^2
      y[i] <- (90 / unc[i])^2
    }
  }
  nchi2 <- sum(x, na.rm = TRUE) / sum(y, na.rm = TRUE)

  return(nchi2)
}


#' @title Quasi Median on a Circle
#' @description Median of orientation data, i.e. pi-periodical data
#' @param x a numeric vector containing the orientations whose median is to be computed
#' @param na.rm logical; if true (default), any NA and NaN's are removed from x before the quantiles are computed.
#' @importFrom pracma atand cosd sind
#' @export
#' @seealso \code{\link[stats]{median}}
#' @references Ratanaruamkarn, S., Niewiadomska-Bugaj, M., Wang, J.-C. (2009). A New Estimator of a Circular Median. Communications in Statistics - Simulation and Computation, 38(6), 1269–1291. https://doi.org/10.1080/03610910902899950
#' @examples
#' circular_quasi_median(x = c(0, 45, 55, 40+180, 50+180))
circular_quasi_median <- function(x, na.rm=TRUE){
  if(na.rm==TRUE){
  } else {
    warning("NA values have been dropped\n")
  }
  x <- sort(x[!is.na(x)])
  n <- length(x)

  if(n%%2 != 0){
    m <- (n-1)/2
    qmed <- pracma::atand(
      pracma::sind(x[m+1])/pracma::cosd(x[m+1])
    )
  } else {
    m <- n/2
    qmed <- pracma::atand(
      (pracma::sind(x[m]) + pracma::sind(x[m+1])) / (pracma::cosd(x[m]) + pracma::cosd(x[m+1]))
    )
  }
  return(qmed)
}



#' @title Quasi Quartile on a Circle
#' @description The \code{"stats::quantile"} equivalent for circular orienation data
#' @param x a numeric vector containing the orientations whose sample quartile is to be computed
#' @param na.rm logical; if true (default), any NA and NaN's are removed from x before the quantiles are computed.
#' @importFrom pracma atand cosd sind
#' @export
#' @seealso \code{\link[stats]{stats}}
#' @references Ratanaruamkarn, S., Niewiadomska-Bugaj, M., Wang, J.-C. (2009). A New Estimator of a Circular Median. Communications in Statistics - Simulation and Computation, 38(6), 1269–1291. https://doi.org/10.1080/03610910902899950
#' @examples
#' circular_quasi_quartile(x = c(0, 45, 55, 40+180, 50+180))
circular_quasi_quartile <- function(x, na.rm=TRUE){
  if(na.rm==TRUE){
  } else {
    warning("NA values have been dropped\n")
  }
  x <- sort(x[!is.na(x)])
  n <- length(x)
  #ms <- 1:n

  med <- circular_quasi_median(x)

  if(n%%4 == 0){
    m <- n/4
    lq <- pracma::atand(
      pracma::sind(x[m+1]) / pracma::cosd(x[m+1])
    )
    uq <- pracma::atand(
      pracma::sind(x[3*m+1]) / pracma::cosd(x[3*m+1])
    )
  }
  if(n%%4 == 1){
    m <- (n-1)/4
    lq <- pracma::atand(
      (3*pracma::sind(x[m]) + pracma::sind(x[m+1])) / (3*pracma::cosd(x[m])+pracma::cosd(x[m+1]))
    )
    uq <- pracma::atand(
      (3*pracma::sind(x[3*m])+pracma::sind(x[3*m+1])) / (3*pracma::cosd(x[3*m])+pracma::cosd(x[3*m+1]))
    )
  }
  if(n%%4 == 2) {
    m <- (n-2)/4
    lq <- pracma::atand( (pracma::sind(x[m])+pracma::sind(x[m+1])) / (pracma::cosd(x[m])+pracma::cosd(x[m+1])) )
    uq <- pracma::atand( (pracma::sind(x[3*m])+pracma::sind(x[3*m+1])) / (pracma::cosd(x[3*m])+pracma::cosd(x[3*m+1])) )
  }
  if(n%%4 == 3) {
    m <- (n-2)/4
    lq <- pracma::atand( (pracma::sind(x[m]) + 3*pracma::sind(x[m+1])) / (pracma::cosd(x[m]) + 3*pracma::cosd(x[m+1])) )
    uq <- pracma::atand( (pracma::sind(x[3*m]) + 3*pracma::sind(x[3*m+1])) / (pracma::cosd(x[3*m]) + 3*pracma::cosd(x[3*m+1])) )
  }

  quantiles <- c(x[1], lq, med, uq, x[length(x)])
  names(quantiles) <- c("0%" ,  "25%" , "50%" , "75%" , "100%")
  return(quantiles)
}


#' @title Quasi Interquartile Range on a Circle
#' @description The \code{"stats::IQR"} equivalent for circular orienation data
#' @param x a numeric vector containing the orientations whose sample quartile is to be computed
#' @param na.rm logical; if true (default), any NA and NaN's are removed from x before the quantiles are computed.
#' @export
#' @references Reiter, K., Heidbach, O., Schmitt, D., Haug, K., Ziegler, M., Moeck, I. (2014). A revised crustal stress orientation database for Canada. Tectonophysics, 636, 111–124. https://doi.org/10.1016/j.tecto.2014.08.006
#' @seealso \code{\link[stats]{IQR}}
#' @examples
#' circular_quasi_interquartile_range(x = c(0, 45, 55, 40+180, 50+180))
circular_quasi_interquartile_range <- function(x, na.rm=TRUE){
  if(na.rm==TRUE){
  } else {
    warning("NA values have been dropped\n")
  }
  x <- sort(x[!is.na(x)])

  quantiles <- circular_quasi_quartile(x)
  qiroc <- as.numeric(quantiles[4] - quantiles[2])
  return(qiroc)
}
