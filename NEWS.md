# tectonicr (development version)

# tectonicr 0.4.9

# tectonicr 0.4.9 _2026-08-28_

* NEW: Added additional distribution functions (Circular Uniform, Wrapped Normal, Wrapped Cauchy) 
for density estimation. The functions are also available for kernel density estimation and calculating the circular mode. 
* NEW: Statistical tests: Angular Randomisation Test of Homogeneity (`ar_test()`), 
Watson-Wheeler permutation test (`waton_wheeler_test_perm()`), Watson 2-sample permutation test `watson_two_test_perm()`, 
and permutation Rayleigh test (`rayleigh_test_perm()`.
* NEW: Kernel density estimation in circular plot can be plotted as filled polygon through option `fill`
* NEW: 2 example files datasets (`homing` and `striae`) added for testing directional data
* Many functions have been optimized and debugged using Claude.AI (including `stress2grid()`, `rayleigh_test()`, `kuiper_test()` and `watson_test()`
* North label "N" in `circular_plot()` (and thus in `rose()`) now a placed `text` instead of `mtext`.

# tectonicr 0.4.8 _2025-12-11_

* NEW: `ortensor2d()` to calculate summary statistics for axial data
* NEW: `data2PoR()`: Convenience function to add PoR coordinates and PoR azimuths to data
* NEW: `geom_azimuth()` and `geom_azimuthpoint()`: convenience functions to plot directions as bars in ggplot2
* NEW: `weighting()` convenience function to assign weightings for angles based on optional algorithms
* "N" in rose diagram now in bold
* stacked dot and density function in rose plots are visualized in expended plotting windows 
* Deviation is now positive for counterclockwise deviation of observed azimuth wrt. predicted azimuth
* `est.kappa()` now uses approximation (as precise as original but faster). Original algorithm now in `est.kappa.MLE()`
* minor bug fix in `relative_rotation()`
* safe condition checks with `isTRUE` or `isFALSE`

# tectonicr 0.4.7 _2025-05-22_

* NEW: compatibility for World Stress Map Database 2025. 
* `download_WSM()` has option `"version"` for the user to decide which WSM version should be downloaded. The most recent 2025 version is the default.
* `parse_quality()` assigns an uncertainty of 90 degrees to X-ranked stress data in the 2025 WSM version.
* NEW: calculate shortest distance to plate boundary (which is not always the best choice! see function description for details.)
* major bug fix in `est.kappa()`: converts to directional data by doubling the angles
* calculating circular mode in `circular_summary()` now optional due to performance issues
* some weighting parameters in `stress2grid()` have been renamed to me more alike with {gstat}
* choice of algorithm to calculate confidence interval now optional in `circular_summary()`
* minor bug fixes

# tectonicr 0.4.6 _2025-02-27_

* performance boost for `stress2grid()` and friends

# tectonicr 0.4.5 _2025-02-26_

* new function `distance_binned_stats()` to calculate summary statistics along 
distance intervals: The function supersedes `distroll_circstats()` and friends 
because it is faster and more flexible
* `PoR_to_geographical()` and `geographical_to_PoR()` now accept `data.frame`, `sf` or raster as input
* new `PoR_azimuth()` function to transform directions in the PoR coordinate system / doesn't need to be Shmax data
* bug fix in `circular_distance()` and `circular_dispersion()` when using `axial=FALSE`
* bug fix in `dvm()` when using `axial=FALSE`

# tectonicr 0.4.4

* additional plate motion models available in `cpm_models`
* Use of {circular} package in `dvm()`, `pvm()`, `qvm()` and `rvm()`

# tectonicr 0.4.3.92
* `PoR_distance()`: function to calculate distance to PoR

# tectonicr 0.4.3

* performance upgrade on spatial interpolation
* spatial interpolation of more summary statistics
* new vignette on spatial analysis
* optional grid lines added to rose

# tectonicr 0.4.0.9001  _2024-09-12_

* bug fix in `load_wsm()`
* CI now double when axial data

# tectonicr 0.4.0  _2024-08-08_

* minor performance upgrade
* CRAN submission

# tectonicr 0.3.11 _2024-07-26_

* bug fixes
* minor performance upgrade
* QQ plot for circular data

# tectonicr 0.3.10 _2024-07-13_

* more statistical estimators
* jittered circular dot plot

# tectonicr 0.3.9 _2024-06-20_

* stacked dots for `rose()` diagram
* density as multiples of a von Mises distribution added for circular plots

# tectonicr 0.3.8 _2024-06-17_

* `superimposed_shmax()` and `superimposed_shmax_PB()` to model the stress 
orientation using multiple plate boundaries

# tectonicr 0.3.7 _2024-06-09_

* download WSM2016 data from GFZ server using `download_WSM2016()`

# tectonicr 0.3.6

* bug fixes in `weighted_rayleigh()`

# tectonicr 0.3.2 _2024-05-27_

* adjusted due to functions' move from `spatstat.geom` to `spatstat.univar`

# tectonicr 0.3.0 _2024-05-14_

* weighting powers added to spatial interpolation `stress2grid()`
* `deviation_norm()` accepts two arguments

# tectonicr 0.2.98 _2024-04-07_

* minor fixes

# tectonicr 0.2.97 _2024-04-07_

* `dispersion_grid()` deprecated and replaced by `kernel_dispersion()` 
(details in vignette E)

# tectonicr 0.2.96 _2023-10-15_

* cran update

# tectonicr 0.2.95 _2023-10-15_

* bug fixes in `rose()` (e.g. symmetrical fans when axial data is plotted)
* weighted rose diagram enabled
* add single line and fans to `rose()`
* add mean and confidence interval to `rose()`
* corrected typos in manual
* estimate kappa of a von Mises distribution:` est.kapp()`
* no doubling of angles when testing Watson distribution `watson_test()`

# tectonicr 0.2.94 _2023-09-25_

* mean direction and spread on `rose()` diagram

# tectonicr 0.2.93 _2023-09-10_

* prepared for CRAN submission
* bootstrap statistics of circular dispersion
* multiple angles as input for circular dispersion and daughter functions
* spatial distribution of the dispersion

# tectonicr 0.2.92 _2023-05-16_

* statistical tests for circular uniformity and goodness-of-fit, e.g. 
`rayleigh_test()`
* bug fixes

# tectonicr 0.2.8 _2023-03-01_

* estimator for the error of predictions: `prd_err()`
* area-weighted `rose()` diagrams

# tectonicr 0.2.7 _2023-01-25_

* optimized some functions for better performance
* bug fixes

# tectonicr 0.2.6 _2023-01-06_

* distance to plate boundary in km
* coordinate transformation using quaternions
* plot the transformed azimuth vs. distance to plate boundary
* quick plotting: `quick_plot()`

# tectonicr 0.2.5

* Calculation of rotation replaced by quaternions to boost performance
* `sp` class output for small circles, great circles and loxodromes deprecated
* new functions added to calculate mean/median stress direction, e.g. 
`circular_mean()`, `circular_median()`

# tectonicr 0.1

* new functions to rotate stress directions and data points into PoR coordinate
system
* calculate distance of data point from plate boundaries
* Added functions to calculate relative plate motions from a set of absolute 
plate motions or different relative plate motions

# PlateTectonicStressR 0.0.1

* New `euler_loxodrome()` function to construct loxodromes directing towards an 
given point or Euler pole.

# PlateTectonicStressR 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.
