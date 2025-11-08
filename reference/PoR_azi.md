# Azimuth Conversion from Geographical to PoR Coordinate Reference System

Transforms azimuths and models the direction of maximum horizontal
stress \\\sigma\_{Hmax}\\ in the Euler pole (Pole of Rotation)
coordinate reference system. When type of plate boundary is given, it
also gives the deviation from the theoretically predicted azimuth of
\\\sigma\_{Hmax}\\, the circular distance, and the normalized \\\chi^2\\
statistics.

## Usage

``` r
PoR_azimuth(x, PoR, axial = TRUE)

PoR_shmax(x, PoR, type = c("none", "in", "out", "right", "left"), axial = TRUE)
```

## Arguments

- x:

  `sf` object or a `data.frame` containing the coordinates of the
  point(s) (`lat`, `lon` columns). `x` must contain the direction of
  \\\sigma\_{Hmax}\\ as column `azi`, its standard deviation (column
  `unc`) is optional).

- PoR:

  `data.frame` or object of class `euler.pole` containing the
  geographical coordinates of the Eule pole.

- axial:

  logical. Whether the azimuth is axial (0-180) or directional (0-360).

- type:

  Character. Type of plate boundary (optional). Can be `"out"`, `"in"`,
  `"right"`, or `"left"` for outward, inward, right-lateral, or
  left-lateral moving plate boundaries, respectively. If `"none"` (the
  default), only the PoR-equivalent azimuth is returned.

## Value

`PoR_azimuth` returns numeric vector of the transformed azimuth in
degrees. `PoR_shmax` returns either a numeric vector of the azimuths in
the transformed coordinate system (in degrees), or a `"data.frame"` with

- `azi.PoR`:

  the transformed azimuths (in degrees),

- `prd`:

  the predicted azimuths (in degrees),

- `dev`:

  the deviation between the transformed and the predicted azimuth in
  degrees (positive for counterclockwise deviation of observed azimuth
  wrt. predicted azimuth),

- `nchisq`:

  the Norm \\\chi^2\\ test statistic, and

- `cdist`:

  the angular distance between the transformed and the predicted
  azimuth.

## Details

The theoretical azimuth of \\\sigma\_{Hmax}\\ in the pole of rotation
reference system is 0 (or 180), 45, 90, 135 degrees if the stress is
sourced by an outward, sinistral, inward, or dextral moving plate
boundary, respectively. directions of \\\sigma\_{Hmax}\\ with respect to
the four plate boundary types.

## References

Stephan, T., Enkelmann, E., and Kroner, U. "Analyzing the horizontal
orientation of the crustal stress adjacent to plate boundaries". *Sci
Rep* 13. 15590 (2023).
[doi:10.1038/s41598-023-42433-2](https://doi.org/10.1038/s41598-023-42433-2)
.

## See also

[`model_shmax()`](https://tobiste.github.io/tectonicr/reference/model_shmax.md)
to compute the theoretical direction of \\\sigma\_{Hmax}\\ in the
geographical reference system.
[`deviation_shmax()`](https://tobiste.github.io/tectonicr/reference/deviation_shmax.md)
to compute the deviation of the modeled direction from the observed
direction of \\\sigma\_{Hmax}\\.
[`norm_chisq()`](https://tobiste.github.io/tectonicr/reference/norm_chisq.md)
to calculate the normalized \\\chi^2\\ statistics.
[`circular_distance()`](https://tobiste.github.io/tectonicr/reference/dispersion.md)
to calculate the angular distance.

## Examples

``` r
data("nuvel1")
# North America relative to Pacific plate:
PoR <- subset(nuvel1, nuvel1$plate.rot == "na")

data("san_andreas")
res <- PoR_shmax(san_andreas, PoR, type = "right")
head(res)
#>      azi.PoR prd       dev     nchisq      cdist
#> 1 173.240166 135 -38.24017 0.18053214 0.38311043
#> 2   1.058195 135 -46.05820 2.21486507 0.51846479
#> 3 147.609558 135 -12.60956 0.01962975 0.04765752
#> 4 163.607388 135 -28.60739 0.10103489 0.22925429
#> 5 152.233017 135 -17.23302 0.03666381 0.08776909
#> 6 150.611386 135 -15.61139 0.03008832 0.07242085
```
