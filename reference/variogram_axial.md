# Empirical trace variogram of axial data, for choosing the model

Bins the empirical variogram gamma(h) = mean(1 - cos(2 \* (a_i - a_j)))
over pairs of observations. Use it to choose `nugget` and `range_km` for
[`stress2grid_kriging()`](https://tobiste.github.io/tectonicr/reference/stress2grid_kriging.md),
e.g. by overlaying the fitted model.

## Usage

``` r
variogram_axial(x, breaks = seq(0, 500, 50))
```

## Arguments

- x:

  `sf` object with `azi` (degrees).

- breaks:

  numeric. Distance bin edges in km.

## Value

data.frame with `lag` (km), `gamma` and `n` (pairs per bin).
