#' @title Normalize angle between two directions
#'
#' @description Normalizes the angle between two directions to the acute angle
#' in between, i.e. angles between 0 and 90 degrees.
#'
#' @author Tobias Stephan
#' @param x numeric vector containing angles in degrees
#' @return numeric vector, acute angles between two directions, i.e. values
#' between 0° and 90°
#' @export
#' @examples
#'
#' deviation_norm(91)
deviation_norm <- function(x) {
  # deviation is between 0 and 90
  if (length(x) > 1) {
    for (i in seq_along(x)) {
      while (abs(x[i]) > 90) {
        x[i] <- 180 - abs(x[i])
      }
    }
  } else {
    while (abs(x) > 90) {
      x <- 180 - abs(x)
    }
  }
  return(abs(x))
}
