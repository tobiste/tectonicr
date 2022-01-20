## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(PlateTectonicStressR)
library(ggplot2) # load ggplot library

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


## ----san_andreas2-------------------------------------------------------------
san.andreas.res <- cbind(
  lat=san.andreas$lat, 
  lon=san.andreas$lon, 
  azi=san.andreas$azi,
  san.andreas.prd, 
  san.andreas.deviation
  )

## ----plates-------------------------------------------------------------------
data('PB2002') # load plate boundary data set

## ----plot, warning=FALSE------------------------------------------------------
ggplot(san.andreas.res) +
  coord_fixed(xlim = range(san.andreas.res$lon), ylim = range(san.andreas.res$lat)) +
  borders() +
  geom_path(data = broom::tidy(PB2002), aes(long, lat, group=group), color = 'red', lwd=3, alpha = .5) +
  geom_spoke(aes(x=lon, y = lat, radius= 1, angle = pracma::deg2rad(-azi+90), color=deviation_norm(dev.ld.ccw)), position = 'center_spoke') +
  scale_color_continuous(type = 'viridis', limits=c(0, 90), name = '|Deviation| in °', breaks = seq(0, 90, 22.5))

## ----san_andreas_nchi2--------------------------------------------------------
san.andreas$unc <- quantise_wsm_quality(san.andreas$quality) # convert the WSM quality ranking into the standard deviation of the measurement
norm_chi2(obs = san.andreas$azi, san.andreas.prd$ld.ccw, unc = san.andreas$unc)

## ----load_nuvel---------------------------------------------------------------
data("nuvel1")
euler <- subset(nuvel1, nuvel1$ID == "na") # North America relative to Pacific plate
euler$angle <- euler$rate 

## ----small_circles_around_ep, warning=FALSE-----------------------------------
euler.sm <- eulerpole_smallcircles(euler)

ggplot() +
 coord_fixed() +
 geom_path(data = broom::tidy(PB2002), aes(long, lat, group=group), color = 'red', lwd=1, alpha = .5) +
 borders() +
 labs(title='North America relative to Pacific plate', subtitle='source: NUVEL1') +
 geom_path(data=broom::tidy(euler.sm), aes(long, lat, group=group, lty='small circles'), color="blue") + 
 geom_point(data=euler, aes(lon, lat), color="blue") +
 geom_point(data=euler, aes(lon+180, -lat), color="blue", size = 2)

## ----great_circles_around_ep, warning=FALSE-----------------------------------
euler.gm <- eulerpole_greatcircles(euler)

ggplot() +
  coord_fixed() +
  geom_path(data = broom::tidy(PB2002), aes(long, lat, group=group), color = 'red', lwd=1, alpha = .5) +
  borders() +
  labs(title='North America relative to Pacific plate', subtitle='source: NUVEL1') +
  geom_path(data=broom::tidy(euler.sm), aes(long, lat, group=group, lty='small circles'), color="blue") +
  geom_path(data=broom::tidy(euler.gm), aes(long, lat, group=group, lty='great circles'), color="blue") +
  geom_point(data=euler, aes(lon, lat), color="blue", size = 2) +
  geom_point(data=euler, aes(lon+180, -lat), color="blue", size = 2)

