# Deviation of Observed and Predicted Directions of Maximum Horizontal Stress

Calculate the angular difference between the observed and modeled
direction of maximum horizontal stresses (\\\sigma\_{Hmax}\\) along
great circles, small circles, and loxodromes of the relative plate
motion's Euler pole

## Usage

``` r
deviation_shmax(prd, obs)
```

## Arguments

- prd:

  `data.frame` containing the modeled azimuths of \\\sigma\_{Hmax}\\,
  i.e. the return object from
  [`model_shmax()`](https://tobiste.github.io/tectonicr/reference/model_shmax.md)

- obs:

  Numeric vector containing the observed azimuth of \\\sigma\_{Hmax}\\,
  same length as `prd`

## Value

An object of class `data.frame`

- dev.gc:

  Deviation of observed stress from modeled \\\sigma\_{Hmax}\\ following
  great circles

- dev.sc:

  Small circles

- dev.ld.cw:

  Clockwise loxodromes

- dev.ld.ccw:

  Counter-clockwise loxodromes

## Details

Deviation is positive for counterclockwise deviation of observed azimuth
wrt. predicted azimuth.

## References

Stephan, T., Enkelmann, E., and Kroner, U. "Analyzing the horizontal
orientation of the crustal stress adjacent to plate boundaries". *Sci
Rep* 13. 15590 (2023).
[doi:10.1038/s41598-023-42433-2](https://doi.org/10.1038/s41598-023-42433-2)
.

## See also

[`model_shmax()`](https://tobiste.github.io/tectonicr/reference/model_shmax.md)
to calculate the theoretical direction of \\\sigma\_{Hmax}\\.

## Author

Tobias Stephan

## Examples

``` r
data("nuvel1")
# North America relative to Pacific plate:
PoR <- subset(nuvel1, nuvel1$plate.rot == "na")

# the point where we want to model the SHmax direction:
point <- data.frame(lat = 45, lon = 20)

prd <- model_shmax(point, PoR)
deviation_shmax(prd, obs = 90)
#>     dev.gc    dev.sc dev.ld.cw dev.ld.ccw
#> 1 42.45436 -47.54564  87.45436  -2.545636
```
