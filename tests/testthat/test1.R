data("PB2002")
data("wsm2016")
data("nuvel1")

# test model_shmax
euler <- subset(nuvel1, nuvel1$ID == "na") # North America relative to Pacific
point <- data.frame(lat = 45, lon = 20)
prd <- model_shmax(point, euler)

# test mistfit_shmax
misfit1 <- misfit_shmax(prd, obs = 90)

# test with data.frames
set.seed(12)
points <- dplyr::slice_sample(wsm2016, n = 10)
prd2 <- model_shmax(points, euler)
misfits2 <- misfit_shmax(prd2, points$azi)
test2 <- norm_chi2(obs = points$azi, prd = prd2$sc, unc = 10)

ep1 <- data.frame(lat = 91, lon = 0, angle = 1)
sm.sf <- eulerpole_smallcircles(ep1)
gc.sf <- eulerpole_greatcircles(ep1)
ld.sf <- eulerpole_loxodromes(ep1, sense = "sinistral")

sm.sp <- eulerpole_smallcircles(ep1, returnclass = "sp")
gc.sp <- eulerpole_greatcircles(ep1, returnclass = "sp")
ld.sp <- eulerpole_loxodromes(ep1, sense = "sinistral", returnclass = "sp")

eulerpole_loxodromes(ep1, sense = "dextral")
eulerpole_loxodromes(ep1, sense = "sinistral")
eulerpole_loxodromes(ep1, angle = 0, sense = "sinistral")
eulerpole_loxodromes(ep1, angle = 90, sense = "sinistral")
eulerpole_loxodromes(ep1, angle = -95, sense = "sinistral")

ep2 <- data.frame(lat = 90, lon = -185, angle = 1)
eulerpole_smallcircles(ep2)
eulerpole_greatcircles(ep2)
eulerpole_loxodromes(ep2, sense = "dextral")

ep3 <- data.frame()
eulerpole_smallcircles(ep3)
eulerpole_greatcircles(ep3)
eulerpole_loxodromes(ep3, sense = "dextral")

p1 <- c(35, 45) # Baghdad
p2 <- c(35, 135) # Osaka
p3 <- c(35, NA) # add NA values
get_azimuth(p3, p2)


test_that("Output of functions is as expected", {
  expect_equal(norm_chi2(NA, NA, NA), NaN)
  expect_equal(norm_chi2(1, NA, NA), NaN)
  expect_equal(longitude_modulo(-361), -1)
  expect_equal(abs_vel(0.21, 0, r = 1), 0)
  expect_equal(quantise_wsm_quality(c("A", "E", "F", "G", 5)), c(15, NA, NA, NA, NA))
  expect_equal(circular_quasi_median(c(15, 16)), 15.5)
  expect_equal(circular_quasi_median(c(15, 15, 16)), 15)
  expect_equal(circular_quasi_interquartile_range(c(15, 16, 15, 15)), 1)
  expect_equal(deviation_norm(91), 89)
})


test_that("Statistics return NULL when too few numbers", {
  expect_null(circular_quasi_quartile(c(15, 16)))
  expect_null(circular_quasi_quartile(c(15, 15, 16)))
})

test_that("type of object returned is as expected", {
  expect_vector(get_azimuth(p1, p2), ptype = double(), size = 1)
  expect_type(rotation_matrix(c(0, 1, 0), 90), "double")
  expect_s3_class(sm.sf, "sf")
  expect_s3_class(gc.sf, "sf")
  expect_s3_class(ld.sf, "sf")
  expect_s4_class(sm.sp, "SpatialLinesDataFrame")
  expect_s4_class(gc.sp, "SpatialLinesDataFrame")
  expect_s4_class(ld.sp, "SpatialLinesDataFrame")
})

test_that("Message expected", {
  expect_message(circular_quasi_median(c(12, NA)))
  expect_message(circular_quasi_quartile(c(12, NA)))
  expect_message(circular_quasi_interquartile_range(c(12, NA)))
})


test_that("Error message if incorrect type argument", {
  expect_error(misfit_shmax(c(1, 2), c(1)))
  expect_error(angle_vectors(NA, NA))
  expect_error(angle_vectors(c(1, 0), NA))
  expect_error(cartesian_to_geographical(1))
  expect_error(geographical_to_cartesian(1))
  expect_error(cartesian_to_geographical(NA))
  expect_error(geographical_to_cartesian(NA))
  expect_error(norm_chi2(obs = c(1, 2), prd = 1, unc = 1))
  expect_error(norm_chi2(obs = c(1, 2), prd = 1, unc = c(1, 2, 3)))
})
