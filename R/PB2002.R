#' Global model of current plate motions
#'
#' PB2002 global model of current plate motions by Bird 2003
#'
#' @docType data
#'
#' @usage data('pb2002')
#'
#' @format An object of class `data.frame`
#' \describe{
#'   \item{plate.name}{The rotating plate}
#'   \item{plate.rot}{The abbreviation of the plate's name}
#'   \item{lat,lon }{Coordinates of the Euler pole}
#'   \item{angle}{The amount of rotation (angle in 1 Myr)}
#'   \item{plate.fix}{The anchored plate, i.e. `plate.rot` moves relative
#'   to `plate.fix`}
#'   \item{source}{Reference to underlying study}
#'   \item{area}{Size of the plate}
#' }
#' @references  Bird, P. (2003), An updated digital model of plate boundaries,
#' *Geochem. Geophys. Geosyst.*, **4**, 1027, doi: 10.1029/2001GC000252, 3.
#' @keywords datasets
#' @examples
#' data("pb2002")
#' head("pb2002")
"pb2002"
