#' Example crustal stress dataset
#'
#' Subsets of the World Stress Map (WSM) compilation of information on the
#' crustal present-day stress field (Version 1.1. 2019).
#'
#' @details
#' \describe{
#'  \item{`'san_andreas"`}{contains 407 stress data adjacent to the San Andreas
#' Fault to be tested against a tangentially displaced plate boundary.}
#' \item{`"tibet"`}{contains 947 stress data from the Himalaya and Tibetan
#' plateau to be tested against an inward-moving displaced plate boundary.}
#' \item{`'iceland`}{contains 201 stress data from Iceland to be tested against a
#' outward-moving displaced plate boundary.}
#' }
#'
#' @docType data
#'
#' @format A `sf` object / `data.frame` with 10 columns. Each row represents a
#' different in-situ stress measurement:
#' \describe{
#'  \item{id}{Measurement identifier}
#'  \item{lat}{Latitude in degrees}
#'  \item{lon}{Longitude in degrees}
#'  \item{azi}{SHmax azimuth in degrees}
#'  \item{unc}{Measurement standard deviation (in degrees)}
#'  \item{type}{Type of measurement}
#'  \item{depth}{Depth in km}
#'  \item{quality}{WSM quality rank}
#'  \item{regime}{Stress regime}
#' }
#'
#' @references
#' Heidbach, O., Barth, A., Müller, B., Reinecker, J., Stephansson, O., Tingay,
#' M., & Zang, A. (2016). WSM quality ranking scheme, database description and
#' analysis guidelines for stress indicator. WSM Technical Report; 16-01.
#' GFZ German Research Centre for Geosciences. \doi{10.2312/WSM.2016.001}
#'
#' @source \url{https://www.world-stress-map.org/}
#'
#' @keywords datasets
#'
#' @seealso [download_WSM()] for description of columns and stress regime
#' acronyms
#'
#' @name stress_data
#'
#' @examples
#' data("san_andreas")
#' head(san_andreas)
#'
#' data("tibet")
#' head(tibet)
#'
#' data("iceland")
#' head(iceland)
NULL


#' @usage data('san_andreas')
#' @rdname stress_data
"san_andreas"

#' @usage data('tibet')
#' @rdname stress_data
"tibet"

#' @usage data('iceland')
#' @rdname stress_data
"iceland"
