# SHmax direction resulting from multiple plate boundaries

Calculates a \\\sigma\_{Hmax}\\ direction at given coordinates, sourced
by multiple plate boundaries. This first-order approximation is the
circular mean of the superimposed theoretical directions, weighted by
the rotation rates of the underlying PoRs.

## Usage

``` r
superimposed_shmax(df, PoRs, types, absolute = TRUE, PoR_weighting = NULL)
```

## Arguments

- df:

  `data.frame` containing the coordinates of the point(s) (`lat`,
  `lon`), and the direction of \\\sigma\_{Hmax}\\ `azi` (in degrees)

- PoRs:

  multirow `data.frame` or `"euler.pole"` object that must contain
  `lat`, `lon` and `angle`

- types:

  character vector with length equal to number of rows in `PoRs`. Type
  of plate boundary. Must be `"out"`, `"in"`, `"right"`, or `"left"` for
  outward, inward, right-lateral, or left-lateral moving plate
  boundaries, respectively.

- absolute:

  logical. Whether the resultant azimuth should be weighted using the
  absolute rotation at the points or the angular rotation of the PoRs.

- PoR_weighting:

  (optional) numeric vector with length equal to number of rows in
  `PoRs`. Extra weightings for the used PoRs.

## Value

two column vector. `azi` is the resultant azimuth in degrees /
geographical CRS), `R` is the resultant length.

## See also

[`model_shmax()`](https://tobiste.github.io/tectonicr/reference/model_shmax.md)

[`superimposed_shmax_PB()`](https://tobiste.github.io/tectonicr/reference/superimposed_shmax_PB.md)
for considering distances to plate boundaries

## Examples

``` r
data(san_andreas)
data(nuvel1)
pors <- subset(nuvel1, plate.rot %in% c("eu", "na"))
res <- superimposed_shmax(san_andreas, pors, types = c("in", "right"), PoR_weighting = c(2, 1))
head(res)
#>           azi         R
#> [1,] 156.6593 0.6609390
#> [2,] 150.5433 0.6532568
#> [3,] 155.9269 0.6555497
#> [4,] 157.0336 0.6434240
#> [5,] 155.9564 0.5960833
#> [6,] 152.3433 0.5636785
```
