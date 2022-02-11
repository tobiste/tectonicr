position_center_spoke <- function() PositionCenterSpoke # position subclass
# "center_spoke" to center ggplot::geom_spoke() marker at its origin

#' @title  Centrically aligned geom_spoke marker
#'
#' @description \code{position} subclass "center_spoke" to center
#' \code{ggplot::geom_spoke()} marker at its origin
#' @export
#' @importFrom ggplot2 ggproto Position
PositionCenterSpoke <- ggplot2::ggproto("PositionCenterSpoke", ggplot2::Position,
  compute_panel = function(self, data, params, scales) {
    data$x <- 2 * data$x - data$xend
    data$y <- 2 * data$y - data$yend
    data$radius <- 2 * data$radius
    data
  }
)

#' Quantize World Stress Map quality ranking
#'
#' Quantize the World Stress Map A, B, C, D quality ranking
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
#' quantise_wsm_quality(c("A", "B", "C", "D", NA))
#' data("wsm2016")
#' quantise_wsm_quality(wsm2016$quality)
quantise_wsm_quality <- function(x) {
  azi.std <- c()
  for (i in seq_along(x)) {
    azi.std[i] <- ifelse(x[i] == "A", 15,
      ifelse(x[i] == "B", 20,
        ifelse(x[i] == "C", 25,
          ifelse(x[i] == "D", 40, NA)
        )
      )
    )
  }
  return(azi.std)
}
