# Plotting Stress Analysis Results

Creates a set of plots including the azimuth as a function of the
distance to the plate boundary, the Norm Chi-squared as a function of
the distance to the plate boundary, the circular distance (and
dispersion) a function of the distance to the plate boundary, a
von-Mises QQ plot, and a rose diagram of the quality-weighted frequency
distribution of the azimuths.

## Usage

``` r
quick_plot(azi, distance, prd, unc = NULL, regime, width = 51)
```

## Arguments

- azi:

  numeric. Azimuth of \\\sigma\_{Hmax}\\

- distance:

  numeric. Distance to plate boundary

- prd:

  numeric. the predicted direction of \\\sigma\_{Hmax}\\

- unc:

  numeric. Uncertainty of observed \\\sigma\_{Hmax}\\, either a numeric
  vector or a number

- regime:

  character vector. The stress regime (following the classification of
  the World Stress Map)

- width:

  integer. window width (in number of observations) for moving average
  of the azimuths, circular dispersion, and Norm Chi-square statistics.
  If `NULL`, an optimal width will be estimated.

## Value

four R base plots

## Details

Plot 1 shows the transformed azimuths as a function of the distance to
the plate boundary. The red line indicates the rolling circular mean,
stippled red lines indicate the 95% confidence interval about the mean.

Plot 2 shows the normalized \\\chi^2\\ statistics as a function of the
distance to the plate boundary. The red line shows the rolling
\\\chi^2\\ statistic.

Plot 3 shows the circular distance of the transformed azimuths to the
predicted azimuth, as a function of the distance to the plate boundary.
The red line shows the rolling circular dispersion about the prediction.

Plot 4 give the rose diagram of the transformed azimuths.

## See also

[`PoR_shmax()`](https://tobiste.github.io/tectonicr/reference/PoR_azi.md),
[`distance_from_pb()`](https://tobiste.github.io/tectonicr/reference/distance_from_pb.md),
[`circular_mean()`](https://tobiste.github.io/tectonicr/reference/circle_stats.md),
[`circular_dispersion()`](https://tobiste.github.io/tectonicr/reference/dispersion.md),
[`confidence_interval_fisher()`](https://tobiste.github.io/tectonicr/reference/confidence_interval_fisher.md),
[`norm_chisq()`](https://tobiste.github.io/tectonicr/reference/norm_chisq.md),
[`weighted_rayleigh()`](https://tobiste.github.io/tectonicr/reference/weighted_rayleigh.md),
[`vm_qqplot()`](https://tobiste.github.io/tectonicr/reference/vm_qqplot.md)

## Examples

``` r
data("nuvel1")
na_pa <- subset(nuvel1, nuvel1$plate.rot == "na")

data("plates")
plate_boundary <- subset(plates, plates$pair == "na-pa")

data("san_andreas")
res <- PoR_shmax(san_andreas, na_pa, "right")
d <- distance_from_pb(san_andreas, na_pa, plate_boundary, tangential = TRUE)
quick_plot(res$azi.PoR,
  distance = d, prd = res$prd, unc = san_andreas$unc,
  regime = san_andreas$regime
)
```
