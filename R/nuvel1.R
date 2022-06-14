#' Global model of current plate motions
#'
#' NUVEL-1 global model of current plate motions by DeMets et al. 1990
#'
#' @docType data
#'
#' @usage data('nuvel1')
#'
#' @format An object of class \code{"data.frame"}
#' \describe{
#'   \item{plate.name}{The rotating plate}
#'   \item{plate.rot}{The abbreviation of the plate's name}
#'   \item{lat,lon }{Coordinates of the Euler pole}
#'   \item{angle}{The amount of rotation (angle in 1 Myr)}
#'   \item{plate.fix}{The anchored plate, i.e. \code{plate.rot} moves relative
#'   to \code{plate.fix}}
#'   \item{source}{Reference to underlying study}
#' }
#' @references DeMets, C., Gordon, R. G., Argus, D. F., & Stein, S. (1990).
#' Current plate motions. *Geophysical Journal International*, 101(2), 425--478.
#' \doi{10.1111/j.1365-246X.1990.tb06579.x}.
#' @keywords datasets
#' @examples
#' data("nuvel1")
#' head("nuvel1")
"nuvel1"
