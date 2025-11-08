# Degrees to Radians

Helper functions to transform between angles in degrees and radians.

## Usage

``` r
rad2deg(rad)

deg2rad(deg)
```

## Arguments

- rad:

  (array of) angles in radians.

- deg:

  (array of) angles in degrees.

## Value

numeric. angle in degrees or radians.

## Examples

``` r
deg2rad(seq(-90, 90, 15))
#>  [1] -1.5707963 -1.3089969 -1.0471976 -0.7853982 -0.5235988 -0.2617994
#>  [7]  0.0000000  0.2617994  0.5235988  0.7853982  1.0471976  1.3089969
#> [13]  1.5707963
rad2deg(seq(-pi / 2, pi / 2, length = 13))
#>  [1] -90 -75 -60 -45 -30 -15   0  15  30  45  60  75  90
```
