# tectonicr (development version)

# tectonicr 0.2.95 _2023-10-15_

* bug fixes in `rose()` (e.g. symmetrical fans when axial data is plotted)
* weighted rose diagram enabled
* add single line and fans to `rose()`
* add mean and confidence interval to `rose()`
* corrected typos in manual
* estimate kappa of von Mises distribution:` est.kapp()`
* no doubling of angles when testing Watson distribution `watson_test()`

# tectonicr 0.2.94 _2023-09-25_

* mean direction and spread on rose diagram (`rose()`)

# tectonicr 0.2.93 _2023-09-10_

* prepared for CRAN submission
* bootstrap statistics of circular dispersion
* multiple angles as input for circular dispersion and daughter functions
* spatial distribution of the dispersion

# tectonicr 0.2.92 _2023-05-16_

* statistical tests for circular uniformity and goodness-of-fit
* bug fixes

# tectonicr 0.2.8 _2023-03-01_

* estimator for the error of predictions
* area-weighted rose diagrams

# tectonicr 0.2.7 _2023-01-25_

* optimized some functions for better performance
* bug fixes

# tectonicr 0.2.6 _2023-01-06_

* distance to plate boundary in km
* coordinate transformation using quaternions
* plot the transformed azimuth vs. distance to plate boundary
* quick plotting

# tectonicr 0.2.5

* Calculation of rotation replaced by quaternions to boost performance
* "sp" class output for small circles, great circles and loxodromes deprecated
* new functions added to calculate mean/median stress direction, e.g. weighted mean

# tectonicr 0.1

* new functions to rotate stress directions and data points into PoR coordinate
system
* calculate distance of data point from plate boundaries
* Added functions to calculate relative plate motions from a set of absolute 
plate motions or different relative plate motions

# PlateTectonicStressR 0.0.1

* New `euler_loxodrome()` function to construct loxodromes directing towards an given point or Euler pole.

# PlateTectonicStressR 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.
