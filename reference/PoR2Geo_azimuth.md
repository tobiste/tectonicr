# Azimuth Conversion From PoR to Geographical Coordinate Reference System

Conversion of PoR azimuths into geographical azimuths

## Usage

``` r
PoR2Geo_azimuth(x, PoR, axial = TRUE)
```

## Arguments

- x:

  `data.frame` containing the PoR equivalent azimuths (`azi.PoR`), and
  either the geographical coordinates of the point(s) or the
  PoR-equivalent coordinates.

- PoR:

  `data.frame` containing the geographical location of the Euler pole
  (`lat`, `lon`)

- axial:

  logical. Whether the azimuth is axial (0-180) or directional (0-360).

## Value

numeric vector of transformed azimuths (in degrees)

## References

Stephan, T., Enkelmann, E., and Kroner, U. "Analyzing the horizontal
orientation of the crustal stress adjacent to plate boundaries". *Sci
Rep* 13. 15590 (2023).
[doi:10.1038/s41598-023-42433-2](https://doi.org/10.1038/s41598-023-42433-2)
.

## See also

[`PoR_shmax()`](https://tobiste.github.io/tectonicr/reference/PoR_azi.md)

## Examples

``` r
data("nuvel1")
# North America relative to Pacific plate:
PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
data("san_andreas")
san_andreas$azi.PoR <- PoR_shmax(san_andreas, PoR)

# convert back to geo CRS
PoR2Geo_azimuth(san_andreas, PoR) |> head()
#> [1] 50 54 24 41 30 27
```
