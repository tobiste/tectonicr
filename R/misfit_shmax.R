#'  Deviation of observed and predicted directions of SHmax
#'
#' Calculate the angular difference between the observed and modeled direction
#' of maximum horizontal stresses along great circles, small circles, and
#' loxodromes of the relative plate motion´s Euler pole
#'
#' @author Copyright (C) 2021 Tobias Stephan
#' @param prd data.frame containing the modeled azimuths of SHmax, i.e. the
#' return object from model_shmax()
#' @param obs numeric vector containing the observed azimuth of SHmax,
#' same length as prd
#' @return data.frame containing the deviation angles between the  observed and
#' the modeled maximum horizontal stresses along great circles (dev.gc),
#' small circles (dev.sc), and loxodromes (dev.ld.cw, dev.ld.ccw)
#' @export
#' @examples

misfit_shmax <- function(prd, obs){
  if(length(obs)!=seq_along(prd$gc)){
    stop('prd and obs must have have the same length')
  }

  dev.gc <- c()
  dev.sc <- c()
  dev.ld.cw <- c()
  dev.ld.ccw <- c()

  for(i in seq_along(prd$gc)){
    if(is.na(obs[i])){
      dev.gc[i] <- NA
      dev.sc[i] <- NA
      dev.ld.cw[i] <- NA
      dev.ld.ccw[i] <- NA
      }
    else {

    # deviation from great circle direction
    dev.gc[i] <- deviation.norm(prd$gc[i] - obs[i])
    if(prd$gc[i] <= obs[i]){
      dev.gc[i] <- -1*dev.gc[i]
    }

    # deviation from small circle direction
    dev.sc[i] <- deviation.norm(prd$sc[i] - obs[i])
    if(prd$sc[i] <= obs[i]){
      dev.sc[i] <- -1*dev.sc[i]
    }

    # Deviation from counterclockwise (ccw) loxodrome
    dev.ld.ccw[i] <- deviation.norm(prd$ld.ccw - obs[i])
    if(prd$ld.ccw[i] <= obs[i]){
      dev.ld.ccw[i] <- -1*dev.ld.ccw[i]
    }

    # Deviation from clockwise (cw) loxodrome
    dev.ld.cw[i] <- deviation.norm(prd$ld.cw[i] - obs[i])
    if(prd$ld.cw[i] <= obs[i]){
      dev.ld.cw[i] <- -1*dev.ld.cw[i]
    }
    }
  }

  dev.df <- data.frame(dev.gc, dev.sc, dev.ld.cw, dev.ld.ccw)
  return(dev.df)
}
