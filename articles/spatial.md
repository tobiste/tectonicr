# Spatial Analysis

This vignette demonstrates some additional spatially interpolated
statistics of a stress field.

``` r

library(tectonicr)
library(ggplot2) # load ggplot library
```

``` r

data("san_andreas")

data("cpm_models")
por <- cpm_models[["NNR-MORVEL56"]] |>
  equivalent_rotation("na", "pa")

plate_boundary <- subset(plates, plates$pair == "na-pa")
```

[`circular_summary()`](https://tobiste.github.io/tectonicr/reference/circular_summary.md)
yields several statistics estimates for a given vector of angles, such
as mean, median, standard deviation, quasi-quantiles, mode, 95%
confidence angle, as wells as the moments (, i.e. 2nd moment = variance,
3rd = skewness, 4th = kurtosis):

``` r

circular_summary(san_andreas$azi, w = weighting(san_andreas$unc))
#>            n         mean           sd          var          25% quasi-median 
#> 1126.0000000   10.8538228   23.8438493    0.2927477   15.0000000   35.5357179 
#>          75%       median           CI     skewness     kurtosis            R 
#>  160.0000000    9.0000000    5.4651431   -0.2710770    1.1239222    0.7072523
```

## Spatial analysis

Spatial analysis and interpolation of stress data using
[`stress2grid_stats()`](https://tobiste.github.io/tectonicr/reference/stress2grid.md)
or
[`PoR_stress2grid_stats()`](https://tobiste.github.io/tectonicr/reference/PoR_stress2grid.md)
(analysis in the PoR coordinate system) uses a moving window with a user
defined cell-size (im km) and calculates the summary statistics within
each cell:

``` r

spatial_stats_R <- PoR_stress2grid_stats(san_andreas, PoR = por, gridsize = 1, R_range = 100)
subset(spatial_stats_R, !is.na(mean)) |> head()
#> Simple feature collection with 6 features and 23 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: -110.8569 ymin: 34.99224 xmax: -107.9222 ymax: 39.97615
#> Geodetic CRS:  WGS 84
#>     lon.PoR   lat.PoR mean.PoR       sd        var  25%.PoR quasi-median.PoR
#> 24 84.83339 -61.15394 164.8719 13.04123 0.09842750 173.0119         152.1271
#> 25 85.83339 -61.15394 167.1085 13.82763 0.10995901 147.5408         150.0343
#> 26 86.83339 -61.15394 164.6122 12.53953 0.09135068 172.9765         152.1271
#> 35 75.83339 -60.15394 170.5340 10.09080 0.06014989 147.5408         155.3596
#> 42 82.83339 -60.15394 169.2616 10.85242 0.06923893 147.5408         155.2982
#> 43 83.83339 -60.15394 167.6901 11.32732 0.07519259 173.3593         152.1271
#>     75%.PoR median.PoR       CI  skewness   kurtosis     meanR   R N       mdr
#> 24 155.0990   163.6913 45.41541 0.6473004  -6.344281 0.9015725 100 5 0.6719937
#> 25 173.2871   168.4892 54.35222 2.1805378 -11.516448 0.8900410 100 4 0.1758576
#> 26 155.9182   163.6913 43.64028 0.4693911  -3.152742 0.9086493 100 5 0.6370203
#> 35 173.2871   168.4892 39.39468 1.9447050  12.268748 0.9398501 100 4 0.6239140
#> 42 173.2871   168.4892 42.39134 1.9469927   7.779855 0.9307611 100 4 0.6932675
#> 43 155.3578   163.6913 39.38947 1.6286915   4.699328 0.9248074 100 5 0.6433647
#>                      geometry      lat       lon     mean      25% quasi-median
#> 24 POINT (-110.1257 39.15542) 39.15542 -110.1257 42.48196 50.62196     29.73716
#> 25 POINT (-110.4569 39.56431) 39.56431 -110.4569 45.38443 25.81673     28.31031
#> 26 POINT (-110.7839 39.97615) 39.97615 -110.7839 43.55496 51.91922     31.06982
#> 35 POINT (-107.9222 34.99224) 34.99224 -107.9222 41.62418 18.63096     26.44979
#> 42 POINT (-110.5034 37.78685) 37.78685 -110.5034 44.89325 23.17242     30.92983
#> 43 POINT (-110.8569 38.19922) 38.19922 -110.8569 43.97151 49.64076     28.40849
#>         75%   median
#> 24 32.70908 41.30134
#> 25 51.56308 46.76515
#> 26 34.86091 42.63400
#> 35 44.37730 39.57938
#> 42 48.91877 44.12084
#> 43 31.63922 39.97267
```

One can also specify a range of cell-sizes for a wavelength analysis:

``` r

spatial_stats <- PoR_stress2grid_stats(san_andreas, PoR = por, gridsize = 1, R_range = seq(50, 350, 100), mode = TRUE)
```

The mean azimuth for each grid cell:

``` r

trajectories <- eulerpole_loxodromes(x = por, n = 40, cw = FALSE)

ggplot(spatial_stats) +
  geom_sf(data = plate_boundary, color = "red") +
  geom_sf(data = trajectories, lty = 2) +
  geom_azimuth(data = san_andreas, aes(lon, lat, angle = azi), radius = .17, linewidth = .5, color = "grey30") +
  geom_azimuth(aes(lon, lat, angle = mean, alpha = sd, color = mdr), radius = .5, lwd = 1) +
  coord_sf(xlim = range(san_andreas$lon), ylim = range(san_andreas$lat)) +
  scale_alpha(name = "Standard deviation", range = c(1, .25)) +
  scale_color_viridis_c(
    "Wavelength\n(R-normalized mean distance)",
    limits = c(0, 1),
    breaks = seq(0, 1, .25)
  ) +
  facet_wrap(~R)
```

![](spatial_files/figure-html/plot-1.png)

To filter the range of search windows to only keep the shortest
wavelength (R) with the least variance for each grid cell, use
[`compact_grid2()`](https://tobiste.github.io/tectonicr/reference/compact-grid.md).

``` r

spatial_stats_comp <- spatial_stats |>
  compact_grid2(var)
```

Interpolated median stress field color-coded by the skewness within each
search window:

``` r

ggplot(spatial_stats_comp) +
  geom_sf(data = plate_boundary, color = "red") +
  geom_sf(data = trajectories, lty = 2) +
  geom_azimuth(data = san_andreas, aes(lon, lat, angle = azi), radius = .15, color = "grey30") +
  geom_azimuth(aes(lon, lat, angle = median, alpha = CI, color = skewness), radius = .25, lwd = 1) +
  coord_sf(xlim = range(san_andreas$lon), ylim = range(san_andreas$lat)) +
  scale_alpha(name = "95% CI", range = c(1, .25)) +
  scale_color_viridis_c(
    "Skewness"
  )
```

![](spatial_files/figure-html/unnamed-chunk-4-1.png)

Interpolated mode of the stress field color-coded by the absolute
kurtosis within each search window:

``` r

ggplot(spatial_stats_comp) +
  geom_sf(data = plate_boundary, color = "red") +
  geom_sf(data = trajectories, lty = 2) +
  geom_azimuth(data = san_andreas, aes(lon, lat, angle = azi), radius = .15, color = "grey30") +
  geom_azimuth(aes(lon, lat, angle = mode, alpha = CI, color = abs(kurtosis)), radius = .25, lwd = 1) +
  coord_sf(xlim = range(san_andreas$lon), ylim = range(san_andreas$lat)) +
  scale_alpha(name = "95% CI", range = c(1, .25)) +
  scale_color_viridis_c(
    "|Kurtosis|"
  )
```

![](spatial_files/figure-html/unnamed-chunk-5-1.png)

## Heat maps for the spatial statistics

[`PoR_stress2grid_stats()`](https://tobiste.github.io/tectonicr/reference/PoR_stress2grid.md)
and
[`stress2grid_stats()`](https://tobiste.github.io/tectonicr/reference/stress2grid.md)
allow to create heatmaps showing the spatial patterns of any desired
statistical estimate (from
[`circular_summary()`](https://tobiste.github.io/tectonicr/reference/circular_summary.md)).
Some examples:

### Spatial central moments

#### Spatial variance

``` r

ggplot(spatial_stats_comp) +
  ggforce::geom_voronoi_tile(
    aes(lon, lat, fill = var),
    max.radius = .7, normalize = FALSE
  ) +
  scale_fill_viridis_c(limits = c(0, 1)) +
  geom_sf(data = plate_boundary, color = "red") +
  geom_sf(data = trajectories, lty = 2) +
  geom_azimuth(
    aes(lon, lat, angle = mean),
    radius = .25, lwd = .2, colour = "white"
  ) +
  coord_sf(xlim = range(san_andreas$lon), ylim = range(san_andreas$lat))
```

![](spatial_files/figure-html/variance-1.png)

#### Skewness:

Skewness is a measure for the asymmetry of the probability distribution.
It can be either counterclockwise or clockwise skewed, hence values can
range between negative and positive numbers, respectively. This can be
best visualized in a diverging color-sequence:

``` r

ggplot(spatial_stats_comp) +
  ggforce::geom_voronoi_tile(
    aes(lon, lat, fill = skewness),
    max.radius = .7, normalize = FALSE
  ) +
  scale_fill_gradient2("|Skewness|", low = "#001260", mid = "#EBE5E0", high = "#590007") +
  geom_sf(data = plate_boundary, color = "red") +
  geom_sf(data = trajectories, lty = 2) +
  geom_azimuth(
    aes(lon, lat, angle = median, alpha = var),
    radius = .2, lwd = .2, colour = "white"
  ) +
  scale_alpha("variance", range = c(1, 0)) +
  coord_sf(xlim = range(san_andreas$lon), ylim = range(san_andreas$lat))
```

![](spatial_files/figure-html/skew-1.png)

#### Kurtosis

Kurtosis is a measure of the “tailedness” of the probability
distribution. Here, colors are in a square-root scale:

``` r

ggplot(spatial_stats_comp) +
  ggforce::geom_voronoi_tile(
    aes(lon, lat, fill = abs(kurtosis)),
    max.radius = .7, normalize = FALSE
  ) +
  scale_fill_viridis_c("|Kurtosis|", transform = "sqrt") +
  geom_sf(data = plate_boundary, color = "red") +
  geom_sf(data = trajectories, lty = 2) +
  geom_azimuth(
    aes(lon, lat, angle = mode, alpha = var),
    radius = .25, lwd = .2, colour = "white"
  ) +
  scale_alpha("variance", range = c(1, 0)) +
  coord_sf(xlim = range(san_andreas$lon), ylim = range(san_andreas$lat))
```

![](spatial_files/figure-html/kurtosis-1.png)

### Kernel dispersion

Another way to analyse spatial misfits is the kernel dispersion,
i.e. the local dispersion within a user-defined window (kernel). The
kernel´s half width can be a single number (km) or a range of widths.
The latter requires to compact the grid result (`x`) to find the
smallest kernel size containing the the least dispersion
(`compact_grid(x, 'dispersion')`).

> It is recommended to calculate the kernel dispersion on PoR
> transformed data to avoid angle distortions due to projections.

``` r

san_andreas_por <- san_andreas
san_andreas_por$azi <- PoR_shmax(san_andreas, por, "right")$azi.PoR # transform to PoR azimuth
san_andreas_por$prd <- 135 # test direction
san_andreas_kdisp <- kernel_dispersion(san_andreas_por, gridsize = 1, R_range = seq(50, 350, 100))
san_andreas_kdisp <- compact_grid(san_andreas_kdisp, "dispersion")

ggplot(san_andreas_kdisp) +
  ggforce::geom_voronoi_tile(
    aes(lon, lat, fill = stat),
    max.radius = .7, normalize = FALSE
  ) +
  scale_fill_viridis_c("Dispersion", limits = c(0, 1)) +
  geom_sf(data = trajectories, lty = 2) +
  geom_azimuth(
    data = san_andreas,
    aes(lon, lat, angle = azi, alpha = unc),
    radius = .25, lwd = .2, colour = "white"
  ) +
  scale_alpha("Standard deviation", range = c(1, .25)) +
  coord_sf(xlim = range(san_andreas$lon), ylim = range(san_andreas$lat))
```

![](spatial_files/figure-html/kernel_disp-1.png)
