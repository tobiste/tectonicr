#' Global model of current plate motions (NUVEL-1)
#'
#' Global model of current plate motions (NUVEL-1) by DeMets et al. 1990
#'
#' @docType data
#'
#' @usage data(nuvel1)
#'
#' @format An object of class \code{"data.frame"}
#'
#' @references DeMets, C., Gordon, R. G., Argus, D. F., & Stein, S. (1990).
#' Current plate motions. Geophysical Journal International, 101(2), 425–478.
#' https://doi.org/10.1111/j.1365-246X.1990.tb06579.x
#' @keywords datasets
#' @examples
#'
#' data(nuvel1)
#' head(nuvel1)
"nuvel1"


#' World Stress Map
#'
#' The World Stress Map (WSM) is a global compilation of information on the
#' crustal present-day stress field. Version 2016
#'
#' @docType data
#'
#' @usage data(wsm2016)
#'
#' @format An object of class \code{"data.frame"}
#'
#' @references Heidbach, O., M. Rajabi, X. Cui, K. Fuchs, B. Müller, J.
#' Reinecker, K. Reiter, M. Tingay, F. Wenzel, F. Xie, M. O. Ziegler,
#' M.-L. Zoback, and M. D. Zoback (2018): The World Stress Map database
#' release 2016: Crustal stress pattern across scales. Tectonophysics, 744,
#' 484-498, doi:10.1016/j.tecto.2018.07.007
#' @keywords datasets
#' @examples
#'
#' data(wsm2016)
#' head(wsm2016)
"wsm2016"


#' Quantise World Stress Map quality ranking
#'
#' Quantise the World Stress Map A, B, C, D quality ranking
#'
#' @param x Either a string or a character vector of WSM quality ranking
#' @return \code{"integer"} or vector of type \code{"integer"}
#' @references Heidbach, O.; Barth, A.; Müller, B.; Reinecker, J.;
#' Stephansson, O.; Tingay, M.; Zang, A. (2016). WSM quality
#' ranking scheme, database description and analysis guidelines for stress
#' indicator. World Stress Map Technical Report 16-01, GFZ German Research
#' Centre for Geosciences. DOI: http://doi.org/10.2312/wsm.2016.001
#' @export
#' @examples
#' quantise_wsm_quality(c('A', 'B', 'C', 'D', NA))
#' data('wsm2016')
#' quantise_wsm_quality(wsm2016$quality)
quantise_wsm_quality <- function(x){
  azi.std <- c()
  for(i in seq_along(x)){
    azi.std[i] <- ifelse(x[i]=='A', 15,
                         ifelse(x[i]=='B', 20,
                                ifelse(x[i]=='C', 25,
                                       ifelse(x[i]=='D', 40, NA)
                                )
                         )
    )
  }
  return(azi.std)
}

