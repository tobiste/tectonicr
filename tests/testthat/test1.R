# test get_azimuth()
p1 <- c(35, 45) # Baghdad
p2 <- c(35, 135) # Osaka
get_azimuth(p1, p2)

p3 <- c(35, NA) # add NA values
get_azimuth(p3, p2)

# test model_shmax
data("nuvel1")
euler <- subset(nuvel1, nuvel1$ID == "na") # North America relative to Pacific
# plate
point <- data.frame(lat = 45, lon = 20)
prd <- model_shmax(point, euler)

# test mistfit_shmax
misfit1 <- misfit_shmax(prd, obs = 90)

# test with data.frames
data("wsm2016")
set.seed(12)
points <- dplyr::slice_sample(wsm2016, n = 10)
prd2 <- model_shmax(points, euler)
misfits2 <- misfit_shmax(prd2, points$azi)

test2 <- norm_chi2(obs = points$azi, prd = prd2$sc, unc = 10)

#misfit_shmax(c(1, 2), c(1))

rotation_matrix(c(0, 1, 0), 90)

longitude_modulo(-361)
abs_vel(0.21, 0, r=1)

quantise_wsm_quality(c('A', 'E', 'F', 'G', 5))

circular_quasi_median(NA)
circular_quasi_median(c(15, 16))
circular_quasi_median(c(15, 15, 16))
circular_quasi_quartile(NA)
circular_quasi_quartile(c(15, 16))
circular_quasi_quartile(c(15, 15, 16))


ep1 <- data.frame(lat = 91, lon = 0, angle = 1)
eulerpole_smallcircles(ep1)
eulerpole_greatcircles(ep1)
#eulerpole_loxodromes(ep1)
eulerpole_loxodromes(ep1, sense = 'dextral')
eulerpole_loxodromes(ep1, sense = 'sinistral')
eulerpole_loxodromes(ep1, angle = 0, sense = 'sinistral')
eulerpole_loxodromes(ep1, angle = 90, sense = 'sinistral')
eulerpole_loxodromes(ep1, angle = -95, sense = 'sinistral')

ep2 <- data.frame(lat = 90, lon = -185, angle = 1)
eulerpole_smallcircles(ep2)
eulerpole_greatcircles(ep2)
eulerpole_loxodromes(ep2, sense = 'dextral')

ep3 <- data.frame()
eulerpole_smallcircles(ep3)
eulerpole_greatcircles(ep3)
eulerpole_loxodromes(ep3, sense = 'dextral')
