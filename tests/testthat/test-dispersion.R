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
