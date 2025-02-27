# tectonicr 0.4.6

# tectonicr 0.4.5 _2025-02-26_

* new function `distance_binned_stats()` to calculate summary statistics along distance intervals: The function superseds which superseded `distroll_circstats()` and friends because it is faster and more flexible
* `PoR_to_geographical()` and `geographical_to_PoR()` now accept `data.frame`, `sf` or raster as input
* new `PoR_azimuth()` function to transform directions in the PoR coordinate system / doesn't need to be Shmax data
* bug fix in `circular_distance()` and `circular_dispersion()` when using `axial=FALSE`
* bug fix in `dvm()` when using `axial=FALSE`

# tectonicr 0.4.4

* additional plate motion models available in `cpm_models`
* Use of {circular} package in `dvm()`, `pvm()`, `qvm()` and `rvm()`

# tectonicr 0.4.3.92
* PoR_distance(): function to calculate distance to PoR

# tectonicr 0.4.3

* performance upgrade on spatial interpolation
* spatial interpolation of more summary statistics
* new vignette on spatial analysis
* optional grid lines added to rose

# tectonicr 0.4.0.9001  _2024-09-12_

* bug fix in load_wsm()
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
