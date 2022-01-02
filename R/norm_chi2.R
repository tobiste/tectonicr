#' Normalized chi-square test
#'
#' A quantitative comparison between the predicted and observed directions of
#' SHmax is obtained by the calculation of the average azimuth and by a
#' normalized chi-square test.
#'
#' @author Copyright (C) 2021 Tobias Stephan
#' @references Wdowinski, S., 1998, A theory of intraplate
#'   tectonics: Journal of Geophysical Research: Solid Earth, v. 103, p.
#'   5037-5059, http://dx.doi.org/10.1029/97JB03390.
#' @param obs observed SHmax, numeric vector
#' @param prd predicted SHmax, numeric vector of length of obs
#' @param unc uncertainty of observed SHmax, either numeric vector of length of
#' obs or a number
#' @return numeric vector
#' @details Test values are between 0 and 1 indicating the quality of the
#' predicted SHmax directions. Low values (<= 0.15) indicate good agreement,
#' high values (>0.7) indicate a systematic misfit between predicted and
#' observed SHmax directions
#' @export
#' @examples
#' data("nuvel1")
#' euler <- subset(nuvel1, nuvel1$ID=='na') #North America relative to Pacific plate
#' point <- data.frame(lat = 45, lon = 20)
#' prd <- model_shmax(point, euler)
#' norm_chi2(obs=90, prd, unc=10)
#'
norm_chi2 <- function(obs, prd, unc){
  if(length(prd) != length(obs)){
    stop("Observed and predicted values must have the same length!")
  }

  if(length(unc) != 1 & length(unc) != length(obs)){
    stop("Uncertainties must be either a numeric value or a numeric value with length of observed values!")
  }

  if(length(unc) == 1){
    unc <- rep(unc, length(obs))
  }

    w <- c()
    x <- c()
    y <- c()
    for(i in 1:length(obs)){
      w[i] <- deviation_norm(dev[i])
      x[i] <- ( w[i] / unc[i] )^2
      y[i] <- ( 90 / unc[i] )^2
    }
    nchi2 <- sum(x) / sum(y)

  return(nchi2)
}

