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
