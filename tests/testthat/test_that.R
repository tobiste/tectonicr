data("pb2002")
data("nuvel1_plates")

# test the examples from manuscript
# Load the San Andreas-Gulf of California dataset:
data("san_andreas")
# Load plate boundary geometries
data("plates")
# Extract boundary between Pacific (pa) and North American (na) plates
na_pa_boundary <- subset(plates, plates$pair == "na-pa")
# Load the current plate motion (cpm) models:
data("cpm_models")
morvel <- cpm_models[["NNR-MORVEL56"]] # select MORVEL model
# Relative plate motion between Pacific and North American plates:
na_pa <- equivalent_rotation(morvel, fixed = "na", rot = "pa")
# Transform stress data set and test against predicted left-lateral tangential plate boundary (left):
stress_analysis(san_andreas, PoR = na_pa, type = "right", pb = na_pa_boundary, plot = FALSE)
# Interpolate the stress field in the PoR coordinate system:
PoR_stress2grid(san_andreas, na_pa)

data("tibet")
eu_in_boundary <- subset(plates, plates$pair == "eu-in")
eu_in <- equivalent_rotation(morvel, fixed = "eu", rot = "in")
stress_analysis(tibet, PoR = eu_in, type = "in", pb = eu_in_boundary, plot = FALSE)
PoR_stress2grid(tibet, eu_in)

data("iceland")
eu_na_boundary <- subset(plates, plates$pair == "eu-na")
eu_na <- equivalent_rotation(morvel, fixed = "na", rot = "eu")
stress_analysis(iceland, PoR = eu_na, type = "out", pb = eu_na_boundary, plot = FALSE)
PoR_stress2grid(iceland, eu_na)


# test model_shmax
euler <- subset(nuvel1, nuvel1$plate.rot == "na") # North America relative to Pacific
point <- data.frame(lat = 45, lon = 20)
prd <- model_shmax(point, euler)

# test mistfit_shmax
misfit1 <- deviation_shmax(prd, obs = 90)

# test with data.frames
set.seed(12)
points <- dplyr::slice_sample(san_andreas, n = 10)
prd2 <- model_shmax(points, euler)
misfits2 <- deviation_shmax(prd2, points$azi)
test2 <- norm_chisq(obs = points$azi, prd = prd2$sc, unc = 10)

ep1 <- data.frame(lat = 91, lon = 0, angle = 1)
sm.sf <- eulerpole_smallcircles(ep1)
gc.sf <- eulerpole_greatcircles(ep1)
ld.sf <- eulerpole_loxodromes(ep1, cw = FALSE)

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

p1 <- c(35, 45) # Baghdad
p2 <- c(35, 135) # Osaka
p3 <- c(35, NA) # add NA values
get_azimuth(p3[1], p3[2], p2[1], p2[2])

ep.geo <- c(20, 33)
q.geo <- c(10, 45)
q.por <- geographical_to_PoR_quat(q.geo, ep.geo)
q.por
PoR_to_geographical_quat(q.por, ep.geo)


euler <- subset(nuvel1, nuvel1$plate.rot == "na")


plate_boundary <- subset(plates, plates$plateA %in% c("na", "pa") &
  plates$plateB %in% c("na", "pa"))

distance_from_pb(
  x = san_andreas, PoR = euler, pb = plate_boundary, tangential = TRUE
)

distance_from_pb(
  x = san_andreas, PoR = euler, pb = plate_boundary, tangential = TRUE, km = TRUE
)

test.vals <- c(175, 179, 2, 4)
test.weights <- 1 / c(5, 1, 2, 4)

# test expected output values --------------------------------------------------
test_that("Output of functions is as expected", {
  expect_equal(longitude_modulo(-361), -1)
  expect_equal(abs_vel(0.21, 0, r = 1), 0)
  expect_equal(as.numeric(parse_wsm_quality(c("A", "E", "F", "G", 5))), c(15, 90, NA, NA, NA))
  expect_equal(circular_median(c(15, 16)), 15.5)
  expect_equal(circular_median(c(15, 15, 16)), 15)
  expect_equal(circular_IQR(c(15, 16, 15, 15)), 1)
  expect_equal(deviation_norm(175, 5), 10)
  expect_equal(deviation_norm(c(-10, 100, 170, 190, 260, 280)), c(10, 80, 10, 10, 80, 80))
  expect_equal(cartesian_to_geographical(c(10, 0, 0)), c(0, 0))
  expect_equal(geographical_to_cartesian(c(90, 0)), c(0, 0, 1))
})

sa.por <- PoR_shmax(san_andreas, na_pa, "right")
test_that("Compe to {circular} package", {
  expect_equal(
    circular_mean(sa.por$azi.PoR),
    ((circular::mean.circular(circular::circular(2 * sa.por$azi.PoR, units = "degrees", modulo = "asis")) |> as.numeric() / 2) %% 180)
  )
})


# test output is NULL ----------------------------------------------------------
# test_that("Statistics return NULL when too few numbers", {
#   expect_null(circular_quantiles(c(15, 16)))
# })

# test type --------------------------------------------------------------------

test_that("type of object returned is as expected", {
  expect_vector(get_azimuth(p1[1], p1[2], p2[1], p2[2]), ptype = double(), size = 1)
  expect_s3_class(sm.sf, "sf")
  expect_s3_class(gc.sf, "sf")
  expect_s3_class(ld.sf, "sf")
})

# test message -----------------------------------------------------------------
test_that("Message expected", {
  # expect_message(norm_chisq(c(12, NA), 1, 1))
  expect_message(circular_quantiles(c(15, 16)))
  expect_message(circular_quantiles(c(15, 15, 16)))
})

# test warning -----------------------------------------------------------------
# test_that("Warning expected", {
#   expect_warning(rotation_angle(rotation_matrix(c(0, 0, 1), 0.000001)))
# })

# test error -------------------------------------------------------------------
test_that("Error message if incorrect type argument", {
  expect_error(deviation_shmax(c(1, 2), c(1)))
  expect_error(cartesian_to_geographical(1))
  expect_error(geographical_to_cartesian(1))
  expect_error(cartesian_to_geographical(NA))
  expect_error(geographical_to_cartesian(NA))
  expect_error(norm_chisq(1, NA, NA))
  expect_error(norm_chisq(NA, NA, NA))
  expect_error(norm_chisq(2, 3, 3, na.rm = "typo"))
  expect_error(norm_chisq(obs = c(1, 2), prd = 1, unc = c(1, 2, 3)))
  expect_error(euler_pole(90, 0, NA, "test"))
  expect_error(circular_quasi_IQR(c(12, NA, 10, 9, "Inf", 7)))
  expect_error(PoR_shmax(stress, 10))
  expect_error(distance_from_pb(san_andreas, euler, plate_boundary, tangential = "typo"))
  expect_error(distance_from_pb(x = stress, PoR = euler, pb = san_andreas, tangential = TRUE))
  expect_error(equivalent_rotation(nuvel1, fixed = "test"))
  expect_error(eulerpole_smallcircles(ep3))
  expect_error(eulerpole_greatcircles(ep3))
  expect_error(eulerpole_loxodromes(ep3))
  expect_error(eulerpole_loxodromes(ep1, angle = 90, cw = FALSE))
  expect_error(eulerpole_paths(ep3))
})


## test azimuth conversion -----------------------------------------------------

## this test works but seems to fail lately on older mac releases :-(
# test_that("Azimuth back conversion", {
#   na_pa <- subset(nuvel1, nuvel1$plate.rot == "na")
#   san_andreas$azi.PoR <- PoR_shmax(san_andreas, na_pa)
#
#   eu_na <- equivalent_rotation(nuvel1, "eu", "na")
#   iceland$azi.PoR <- PoR_shmax(iceland, eu_na)
#
#   eu_in <- equivalent_rotation(nuvel1, "eu", "in")
#   tibet$azi.PoR <- PoR_shmax(tibet, eu_in)
#
#   expect_equal(PoR2Geo_azimuth(san_andreas, na_pa), san_andreas$azi %% 180)
#   expect_equal(PoR2Geo_azimuth(iceland, eu_na), iceland$azi %% 180)
#   expect_equal(PoR2Geo_azimuth(tibet, eu_in), tibet$azi %% 180)
#
#   san_andreas_por <- geographical_to_PoR(san_andreas, na_pa)
#   por_crds <- sf::st_coordinates(san_andreas_por) |> as.data.frame()
#   san_andreas_por$lat.PoR <- por_crds$Y
#   san_andreas_por$lon.PoR <- por_crds$X
#
#   expect_equal(round(PoR2Geo_azimuth(san_andreas_por, na_pa), 12), san_andreas$azi %% 180)
# })

## tes coordinates -------------------------------------------------------------

test_that("Coordinate conversion sf", {
  san_andreas_por <- geographical_to_PoR(san_andreas, na_pa) |>
    PoR_to_geographical_sf(na_pa)
  por_crds <- sf::st_coordinates(san_andreas_por)
  geo_crds <- sf::st_coordinates(san_andreas)

  expect_equal(por_crds, geo_crds)
})

test_that("Coordinate conversion df", {
  geo_crds <- sf::st_coordinates(san_andreas)

  san_andreas_por2 <- geographical_to_PoR(san_andreas, na_pa) |>
    PoR_to_geographical(na_pa)
  por_crds2 <- cbind(X = san_andreas_por2$lon, Y = san_andreas_por2$lat)

  expect_equal(por_crds2, geo_crds)
})


## test dispersion -------------------------------------------------------------

test_that("Max Circular dispersion is 1!", {
  expect_equal(circular_dispersion(c(0, 180), 90, axial = TRUE), 1)
  expect_equal(circular_dispersion(c(0, 180), 90, axial = FALSE), .5)
  expect_equal(circular_dispersion(rep(270, 2), 90, axial = TRUE), 0)
  expect_equal(circular_dispersion(rep(270, 2), 90, axial = FALSE), 1)
})

test_that("Max Circular distance is 1!", {
  expect_equal(circular_distance(0, 90, axial = TRUE), 1)
  expect_equal(circular_distance(0, 90, axial = FALSE), .5)
  expect_equal(circular_distance(270, 90, axial = TRUE), 0)
  expect_equal(circular_distance(270, 90, axial = FALSE), 1)
})
