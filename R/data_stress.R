#' Example crustal stress dataset
#'
#' Subsets of the World Stress Map (WSM) compilation of information on the
#' crustal present-day stress field (Version 1.1. 2019).
#'
#' @details
#' \describe{
#'  \item{`'san_andreas"`}{contains 407 stress data adjacent to the San Andreas
#' Fault to be tested against a tangentially displaced plate boundary.}
#' \item{`"tibet"`}{ contains 947 stress data from the Himalaya and Tibetan
#' plateau to be tested against an inward-moving displaced plate boundary.}
#' \item{`'iceland`}{contains 201 stress data from Iceland to be tested against a
#' outward-moving displaced plate boundary.}
#' }
#'
#'
#' @docType data
#'
#' @format A `sf` object / `data.frame` with 10 columns. Each row represents a different in-situ stress measurement:
#' \describe{
#'  \item{id}{Measurement identifier}
#'  \item{lat}{latitude in degrees}
#'  \item{lon}{longitude in degrees}
#'  \item{azi}{SHmax azimuth in degrees}
#'  \item{unc}{MMeasruement standard eviation (in degrees)}
#'  \item{type}{Type of measurement}
#'  \item{depth}{Depth in km}
#'  \item{quality}{WSM quality rank}
#'  \item{regime}{Stress regime}
#' }
#'
#' @references
#' Heidbach, Oliver; Rajabi, Mojtaba; Di Giacomo, Domenico; Harris, James;
#' Lammers, Steffi; Morawietz, Sophia; Pierdominici, Simona; Reiter, Karsten;
#' von Specht, Sebastian; Storchak, Dmitry; Ziegler, Moritz O. (2025): World
#' Stress Map Database Release 2025. GFZ Data Services.
#' \doi{10.5880/WSM.2025.001}
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
