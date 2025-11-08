# Strike of the plate boundary projected on data point

The fault's strike in the PoR CRS projected on the data point along the
predicted stress trajectories.

## Usage

``` r
projected_pb_strike(x, PoR, pb, tangential = FALSE, ...)
```

## Arguments

- x, pb:

  `sf` objects of the data points and the plate boundary geometries in
  the geographical coordinate system

- PoR:

  Pole of rotation. `"data.frame"` or object of class `"euler.pole"`
  containing the geographical coordinates of the Euler pole

- tangential:

  Logical. Whether the plate boundary is a tangential boundary (`TRUE`)
  or an inward and outward boundary (`FALSE`, the default).

- ...:

  optional arguments passed to
  [`smoothr::densify()`](https://strimas.com/smoothr/reference/densify.html)

## Value

Numeric vector of the strike direction of the plate boundary (in degree)

## Details

Useful to calculate the beta angle, i.e. the angle between SHmax
direction (in PoR CRS!) and the fault's strike (in PoR CRS). The beta
angle is the same in geographical and PoR coordinates.

## Note

The algorithm calculates the great circle bearing between line vertices.
Since transform plate boundaries represent small circle lines in the PoR
system, this great-circle azimuth is only a approximation of the true
(small-circle) azimuth.

## Examples

``` r
data("nuvel1")
na_pa <- subset(nuvel1, nuvel1$plate.rot == "na")

data("plates")
plate_boundary <- subset(plates, plates$pair == "na-pa")

data("san_andreas")
res <- projected_pb_strike(
  x = san_andreas, PoR = na_pa, pb = plate_boundary, tangential = TRUE
)
head(res)
#> [1] 76.74072 69.04668 82.89240 76.75248 84.79126 60.25133
head(san_andreas$azi - res) # beta angle
#> [1] -26.74072 -15.04668 -58.89240 -35.75248 -54.79126 -33.25133
```
