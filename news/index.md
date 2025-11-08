# Changelog

## tectonicr (development version)

## tectonicr 0.4.7 *2025-05-22*

CRAN release: 2025-05-22

- NEW: compatibility for World Stress Map Database 2025.
- [`download_WSM()`](https://tobiste.github.io/tectonicr/reference/import_WSM.md)
  has option `"version"` for the user to decide which WSM version should
  be downloaded. The most recent 2025 version is the default.
- `parse_quality()` assigns an uncertainty of 90 degrees to X-ranked
  stress data in the 2025 WSM version.
- NEW: calculate shortest distance to plate boundary (which is not
  always the best choice! see function description for details.)
- major bug fix in
  [`est.kappa()`](https://tobiste.github.io/tectonicr/reference/est.kappa.md):
  converts to directional data by doubling the angles
- calculating circular mode in
  [`circular_summary()`](https://tobiste.github.io/tectonicr/reference/circular_summary.md)
  now optional due to performance issues
- some weighting parameters in
  [`stress2grid()`](https://tobiste.github.io/tectonicr/reference/stress2grid.md)
  have been renamed to me more alike with {gstat}
- choice of algorithm to calculate confidence interval now optional in
  [`circular_summary()`](https://tobiste.github.io/tectonicr/reference/circular_summary.md)
- minor bug fixes

## tectonicr 0.4.6 *2025-02-27*

CRAN release: 2025-03-01

- performance boost for
  [`stress2grid()`](https://tobiste.github.io/tectonicr/reference/stress2grid.md)
  and friends

## tectonicr 0.4.5 *2025-02-26*

CRAN release: 2025-02-26

- new function
  [`distance_binned_stats()`](https://tobiste.github.io/tectonicr/reference/distance_binned_stats.md)
  to calculate summary statistics along distance intervals: The function
  supersedes
  [`distroll_circstats()`](https://tobiste.github.io/tectonicr/reference/rolling_test_dist.md)
  and friends because it is faster and more flexible
- [`PoR_to_geographical()`](https://tobiste.github.io/tectonicr/reference/por_transformation.md)
  and
  [`geographical_to_PoR()`](https://tobiste.github.io/tectonicr/reference/por_transformation.md)
  now accept `data.frame`, `sf` or raster as input
- new
  [`PoR_azimuth()`](https://tobiste.github.io/tectonicr/reference/PoR_azi.md)
  function to transform directions in the PoR coordinate system /
  doesn’t need to be Shmax data
- bug fix in
  [`circular_distance()`](https://tobiste.github.io/tectonicr/reference/dispersion.md)
  and
  [`circular_dispersion()`](https://tobiste.github.io/tectonicr/reference/dispersion.md)
  when using `axial=FALSE`
- bug fix in
  [`dvm()`](https://tobiste.github.io/tectonicr/reference/vonmises.md)
  when using `axial=FALSE`

## tectonicr 0.4.4

CRAN release: 2024-12-11

- additional plate motion models available in `cpm_models`
- Use of {circular} package in
  [`dvm()`](https://tobiste.github.io/tectonicr/reference/vonmises.md),
  [`pvm()`](https://tobiste.github.io/tectonicr/reference/vonmises.md),
  [`qvm()`](https://tobiste.github.io/tectonicr/reference/vonmises.md)
  and
  [`rvm()`](https://tobiste.github.io/tectonicr/reference/vonmises.md)

## tectonicr 0.4.3.92

- [`PoR_distance()`](https://tobiste.github.io/tectonicr/reference/PoR_distance.md):
  function to calculate distance to PoR

## tectonicr 0.4.3

- performance upgrade on spatial interpolation
- spatial interpolation of more summary statistics
- new vignette on spatial analysis
- optional grid lines added to rose

## tectonicr 0.4.0.9001 *2024-09-12*

- bug fix in `load_wsm()`
- CI now double when axial data

## tectonicr 0.4.0 *2024-08-08*

CRAN release: 2024-08-21

- minor performance upgrade
- CRAN submission

## tectonicr 0.3.11 *2024-07-26*

- bug fixes
- minor performance upgrade
- QQ plot for circular data

## tectonicr 0.3.10 *2024-07-13*

- more statistical estimators
- jittered circular dot plot

## tectonicr 0.3.9 *2024-06-20*

- stacked dots for
  [`rose()`](https://tobiste.github.io/tectonicr/reference/rose.md)
  diagram
- density as multiples of a von Mises distribution added for circular
  plots

## tectonicr 0.3.8 *2024-06-17*

- [`superimposed_shmax()`](https://tobiste.github.io/tectonicr/reference/superimposed_shmax.md)
  and
  [`superimposed_shmax_PB()`](https://tobiste.github.io/tectonicr/reference/superimposed_shmax_PB.md)
  to model the stress orientation using multiple plate boundaries

## tectonicr 0.3.7 *2024-06-09*

- download WSM2016 data from GFZ server using
  [`download_WSM2016()`](https://tobiste.github.io/tectonicr/reference/import_WSM.md)

## tectonicr 0.3.6

- bug fixes in
  [`weighted_rayleigh()`](https://tobiste.github.io/tectonicr/reference/weighted_rayleigh.md)

## tectonicr 0.3.2 *2024-05-27*

CRAN release: 2024-05-27

- adjusted due to functions’ move from `spatstat.geom` to
  `spatstat.univar`

## tectonicr 0.3.0 *2024-05-14*

- weighting powers added to spatial interpolation
  [`stress2grid()`](https://tobiste.github.io/tectonicr/reference/stress2grid.md)
- [`deviation_norm()`](https://tobiste.github.io/tectonicr/reference/deviation_norm.md)
  accepts two arguments

## tectonicr 0.2.98 *2024-04-07*

- minor fixes

## tectonicr 0.2.97 *2024-04-07*

- [`dispersion_grid()`](https://tobiste.github.io/tectonicr/reference/kernel_dispersion.md)
  deprecated and replaced by
  [`kernel_dispersion()`](https://tobiste.github.io/tectonicr/reference/kernel_dispersion.md)
  (details in vignette E)

## tectonicr 0.2.96 *2023-10-15*

- cran update

## tectonicr 0.2.95 *2023-10-15*

CRAN release: 2023-11-02

- bug fixes in
  [`rose()`](https://tobiste.github.io/tectonicr/reference/rose.md)
  (e.g. symmetrical fans when axial data is plotted)
- weighted rose diagram enabled
- add single line and fans to
  [`rose()`](https://tobiste.github.io/tectonicr/reference/rose.md)
- add mean and confidence interval to
  [`rose()`](https://tobiste.github.io/tectonicr/reference/rose.md)
- corrected typos in manual
- estimate kappa of a von Mises distribution:`est.kapp()`
- no doubling of angles when testing Watson distribution
  [`watson_test()`](https://tobiste.github.io/tectonicr/reference/watson_test.md)

## tectonicr 0.2.94 *2023-09-25*

- mean direction and spread on
  [`rose()`](https://tobiste.github.io/tectonicr/reference/rose.md)
  diagram

## tectonicr 0.2.93 *2023-09-10*

CRAN release: 2023-09-22

- prepared for CRAN submission
- bootstrap statistics of circular dispersion
- multiple angles as input for circular dispersion and daughter
  functions
- spatial distribution of the dispersion

## tectonicr 0.2.92 *2023-05-16*

- statistical tests for circular uniformity and goodness-of-fit, e.g. 
  [`rayleigh_test()`](https://tobiste.github.io/tectonicr/reference/rayleigh_test.md)
- bug fixes

## tectonicr 0.2.8 *2023-03-01*

- estimator for the error of predictions:
  [`prd_err()`](https://tobiste.github.io/tectonicr/reference/prd_err.md)
- area-weighted
  [`rose()`](https://tobiste.github.io/tectonicr/reference/rose.md)
  diagrams

## tectonicr 0.2.7 *2023-01-25*

- optimized some functions for better performance
- bug fixes

## tectonicr 0.2.6 *2023-01-06*

- distance to plate boundary in km
- coordinate transformation using quaternions
- plot the transformed azimuth vs. distance to plate boundary
- quick plotting:
  [`quick_plot()`](https://tobiste.github.io/tectonicr/reference/quick_plot.md)

## tectonicr 0.2.5

- Calculation of rotation replaced by quaternions to boost performance
- `sp` class output for small circles, great circles and loxodromes
  deprecated
- new functions added to calculate mean/median stress direction, e.g. 
  [`circular_mean()`](https://tobiste.github.io/tectonicr/reference/circle_stats.md),
  [`circular_median()`](https://tobiste.github.io/tectonicr/reference/circle_stats.md)

## tectonicr 0.1

- new functions to rotate stress directions and data points into PoR
  coordinate system
- calculate distance of data point from plate boundaries
- Added functions to calculate relative plate motions from a set of
  absolute plate motions or different relative plate motions
