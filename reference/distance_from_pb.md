# Distance from plate boundary

Absolute distance of data points from the nearest plate boundary

## Usage

``` r
distance_from_pb(x, PoR, pb, tangential = FALSE, km = FALSE, ...)
```

## Arguments

- x:

  `sf` or `data.frame` objects of the data points in geographical
  coordinate system

- PoR:

  Pole of Rotation. `"data.frame"` or object of class `"euler.pole"`
  containing the geographical coordinates of the Pole of Rotation

- pb:

  `sf` objects of the plate boundary geometries in the geographical
  coordinate system

- tangential:

  Logical. Whether the plate boundary is a tangential boundary (`TRUE`)
  or an inward and outward boundary (`FALSE`, the default).

- km:

  Logical. Whether the distance is expressed in kilometers (`TRUE`) or
  in degrees (`FALSE`, the default).

- ...:

  optional arguments passed to
  [`smoothr::densify()`](https://strimas.com/smoothr/reference/densify.html)

## Value

Numeric vector of the great circle distances in units defined by `km`.

## Details

The distance to the plate boundary is the longitudinal or latitudinal
difference between the data point and the plate boundary (along the
closest latitude or longitude) for inward/outward or tangential plate
boundaries, respectively.

## Note

Stresses emanate from the plate boundary along great circles, small
circles or loxodromes associated with the pole of rotation. Hence the
emanation distance is not necessarily the shortest distance to the plate
boundary, which is measured along a great circle unrelated to the pole
of rotation. The differences are particularly notable when the plate
boundary is kinked or for convergent and divergent plate boundaries.

## References

Wdowinski, S. (1998). A theory of intraplate tectonics. Journal of
Geophysical Research: Solid Earth, 103(3), 5037\<U+2013\>5059.
http://dx.doi.org/10.1029/97JB03390

## Examples

``` r
data("nuvel1")
na_pa <- subset(nuvel1, nuvel1$plate.rot == "na")

data("plates")
plate_boundary <- subset(plates, plates$pair == "na-pa")

data("san_andreas")
res <- distance_from_pb(
  x = san_andreas, PoR = na_pa, pb = plate_boundary, tangential = TRUE
)
head(res)
#> [1] -2.446542 -2.751949 -2.705483 -3.162091 -4.988116 -6.358818

res.km <- distance_from_pb(
  x = san_andreas, PoR = na_pa, pb = plate_boundary, tangential = TRUE, km = TRUE
)
range(res.km)
#> [1] -1025.9419   597.3071
```
