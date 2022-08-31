data("pb2002")
data("plates")
# data("wsm2016")
data("nuvel1")
data("nuvel1_plates")
data("san_andreas")
# test rotations



# test model_shmax
euler <- subset(nuvel1, nuvel1$ID == "na") # North America relative to Pacific
point <- data.frame(lat = 45, lon = 20)
prd <- model_shmax(point, euler)

# test mistfit_shmax
misfit1 <- misfit_shmax(prd, obs = 90)

# test with data.frames
set.seed(12)
points <- dplyr::slice_sample(san_andreas, n = 10)
prd2 <- model_shmax(points, euler)
misfits2 <- misfit_shmax(prd2, points$azi)
test2 <- norm_chisq(obs = points$azi, prd = prd2$sc, unc = 10)

ep1 <- data.frame(lat = 91, lon = 0, angle = 1)
sm.sf <- eulerpole_smallcircles(ep1)
gc.sf <- eulerpole_greatcircles(ep1)
ld.sf <- eulerpole_loxodromes(ep1, cw = FALSE)

# sm.sp <- eulerpole_smallcircles(ep1, returnclass = "sp")
# gc.sp <- eulerpole_greatcircles(ep1, returnclass = "sp")
# ld.sp <- eulerpole_loxodromes(ep1, cw = FALSE, returnclass = "sp")

eulerpole_loxodromes(ep1, cw = TRUE)
eulerpole_loxodromes(ep1, cw = FALSE)
eulerpole_loxodromes(ep1, angle = 0, cw = FALSE)
eulerpole_loxodromes(ep1, angle = -95, cw = FALSE)

eulerpole_paths(ep1)

ep2 <- data.frame(lat = 90, lon = -185, angle = 1)
eulerpole_smallcircles(ep2)
eulerpole_greatcircles(ep2)
eulerpole_loxodromes(ep2, cw = TRUE)

ep3 <- data.frame()
# eulerpole_smallcircles(ep3)
# eulerpole_greatcircles(ep3)
# eulerpole_loxodromes(ep3, cw = TRUE)
# eulerpole_paths(ep3)

# euler_rot(ep1, 45)

p1 <- c(35, 45) # Baghdad
p2 <- c(35, 135) # Osaka
p3 <- c(35, NA) # add NA values
get_azimuth(p3, p2)



euler <- subset(nuvel1, nuvel1$plate.rot == "na")


plate_boundary <- subset(plates, plates$plateA %in% c("na", "pa") &
  plates$plateB %in% c("na", "pa"))

distance_from_pb(
  x = san_andreas, ep = euler, pb = plate_boundary, tangential = TRUE
)

distance_from_pb(
  x = san_andreas, ep = euler, pb = plate_boundary, tangential = TRUE, km = TRUE
)

test.vals <- c(175, 179, 2, 4)
test.weights <- 1 / c(5, 1, 2, 4)

# test expected output values --------------------------------------------------
test_that("Output of functions is as expected", {
  expect_equal(longitude_modulo(-361), -1)
  expect_equal(abs_vel(0.21, 0, r = 1), 0)
  expect_equal(quantise_wsm_quality(c("A", "E", "F", "G", 5)), c(15, NA, NA, NA, NA))
  expect_equal(circular_quasi_median(c(15, 16)), 15.5)
  expect_equal(circular_quasi_median(c(15, 15, 16)), 15)
  expect_equal(circular_quasi_IQR(c(15, 16, 15, 15)), 1)
  expect_equal(deviation_norm(91), 89)
  # expect_equal(circular_quasi_median(c(12, NA)), 12)
  # expect_equal(circular_quasi_median(c(180, 0)), 0)
  # expect_equal(circular_weighted_mean(c(180, 0)), 0)
  expect_equal(cartesian_to_geographical(c(10, 0, 0)), c(0, 0))
  expect_equal(geographical_to_cartesian(c(90, 0)), c(0, 0, 1))
  # expect_equal(circular_quasi_median(test.vals), circular_weighted_median(test.vals))
  # expect_equal(circular_quasi_IQR(test.vals), circular_weighted_IQR(test.vals))
  # expect_equal(circular_mean(test.vals), circular_weighted_mean(test.vals))
})

# test output is NULL ----------------------------------------------------------

test_that("Statistics return NULL when too few numbers", {
  expect_null(circular_quasi_quantile(c(15, 16)))
  expect_null(circular_quasi_quantile(c(15, 15, 16)))
})

# test type --------------------------------------------------------------------

test_that("type of object returned is as expected", {
  expect_vector(get_azimuth(p1, p2), ptype = double(), size = 1)
  # expect_type(rotation_matrix(c(0, 1, 0), 90), "double")
  expect_s3_class(sm.sf, "sf")
  expect_s3_class(gc.sf, "sf")
  expect_s3_class(ld.sf, "sf")
  # expect_s4_class(sm.sp, "SpatialLinesDataFrame")
  # expect_s4_class(gc.sp, "SpatialLinesDataFrame")
  # expect_s4_class(ld.sp, "SpatialLinesDataFrame")
})

# test message -----------------------------------------------------------------
test_that("Message expected", {
  expect_message(circular_quasi_quantile(c(12, NA, 10, 11, 9)))
  expect_message(norm_chisq(c(12, NA), 1, 1))
})

# test warning -----------------------------------------------------------------
# test_that("Warning expected", {
#   expect_warning(rotation_angle(rotation_matrix(c(0, 0, 1), 0.000001)))
# })

# test error -------------------------------------------------------------------
test_that("Error message if incorrect type argument", {
  expect_error(misfit_shmax(c(1, 2), c(1)))
  expect_error(cartesian_to_geographical(1))
  expect_error(geographical_to_cartesian(1))
  expect_error(cartesian_to_geographical(NA))
  expect_error(geographical_to_cartesian(NA))
  expect_error(norm_chisq(1, NA, NA))
  expect_error(norm_chisq(NA, NA, NA))
  expect_error(norm_chisq(2, 3, 3, na.rm = "typo"))
  expect_error(norm_chisq(obs = c(1, 2), prd = 1, unc = c(1, 2, 3)))
  # expect_error(rotation_angle(as.character(rotation_matrix(c(0, 0, 1))), 1))
  # expect_error(as.character(rotation_axis(c(0, 0, 1)), 1))
  # expect_error(as.character(rotation_matrix(c(0, 0, 1)), 1))
  expect_error(euler_pole(90, 0, NA, "test"))
  # expect_error(euler_from_rot(C(1, 2, 3)))
  expect_error(circular_quasi_IQR(c(12, NA, 10, 9, "Inf", 7)))
  expect_error(PoR_shmax(stress, 10))
  # expect_error(euler_rot(c(90, 0), "test"))
  expect_error(distance_from_pb(san_andreas, euler, plate_boundary, tangential = "typo"))
  expect_error(distance_from_pb(x = stress, ep = euler, pb = san_andreas, tangential = TRUE))
  expect_error(equivalent_rotation(nuvel1, fixed = "test"))
  expect_error(eulerpole_smallcircles(ep3))
  expect_error(eulerpole_greatcircles(ep3))
  expect_error(eulerpole_loxodromes(ep3))
  expect_error(eulerpole_loxodromes(ep1, angle = 90, cw = FALSE))
  expect_error(eulerpole_paths(ep3))
})
