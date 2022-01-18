## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(PlateTectonicStressR)

## ----direction_of_plate_motion------------------------------------------------
# Example:
point <- data.frame(lat = 45, lon = 20)
euler <- data.frame(lat = 90, lon = 0)
model <- model_shmax(point, euler)
print(model)

## ----deviation_of_plate_motion------------------------------------------------
deviation <- misfit_shmax(model, 90)
print(deviation)

## ----shmax_test---------------------------------------------------------------
data("nuvel1")
euler <- subset(nuvel1, nuvel1$ID == "na") # North America relative to Pacific plate
point <- data.frame(lat = 45, lon = 20)
prd <- model_shmax(point, euler)
norm_chi2(obs = 90, prd$sc, unc = 10)

## ----load_wsm2016-------------------------------------------------------------
data('wsm2016')
head(wsm2016)

## ----wsm_subset---------------------------------------------------------------
san.andreas <- subset(wsm2016, 
                      wsm2016$lat>=23 & wsm2016$lat<= 40 &
                      wsm2016$lon>=-126 & wsm2016$lon<=-108)

## ----san_andreas--------------------------------------------------------------
data("nuvel1")
euler <- subset(nuvel1, nuvel1$ID == "na")

san.andreas.prd <- model_shmax(san.andreas, euler)
san.andreas.deviation <- misfit_shmax(san.andreas.prd, san.andreas$azi)


