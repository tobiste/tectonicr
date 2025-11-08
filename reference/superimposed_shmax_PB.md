# SHmax direction resulting from multiple plate boundaries considering distance to plate boundaries

Calculates a \\\sigma\_{Hmax}\\ direction at given coordinates, sourced
by multiple plate boundaries. This first-order approximation is the
circular mean of the superimposed theoretical directions, weighted by
the rotation rates of the underlying PoRs, the inverse distance to the
plate boundaries, and the type of plate boundary.

## Usage

``` r
superimposed_shmax_PB(
  x,
  pbs,
  model,
  rotation_weighting = TRUE,
  type_weights = c(divergent = 1, convergent = 3, transform_L = 2, transform_R = 2),
  idp = 1
)
```

## Arguments

- x:

  grid. An object of `sf`, `sfc` or 2-column matrix

- pbs:

  plate boundaries. `sf` object

- model:

  `data.frame` containing the Euler pole parameters. See
  [`equivalent_rotation()`](https://tobiste.github.io/tectonicr/reference/equivalent_rotation.md)
  for details.

- rotation_weighting:

  logical.

- type_weights:

  named vector.

- idp:

  numeric. Weighting power of inverse distance. The higher the number,
  the less impact far-distant boundaries have. When set to `0`, no
  weighting is applied.

## Value

two-column matrix. `azi` is the resultant azimuth (in degrees), `R` is
the resultant length.

## See also

[`superimposed_shmax()`](https://tobiste.github.io/tectonicr/reference/superimposed_shmax.md)

## Examples

``` r
na_grid <- sf::st_make_grid(san_andreas, what = "centers", cellsize = 1)
na_plate <- subset(plates, plateA == "na" | plateB == "na")
cpm <- cpm_models[["NNR-MORVEL56"]]

# make divergent to ridge-push:
na_plate <- transform(na_plate, type = ifelse(na_plate$pair == "eu-na", "convergent", type))

res <- superimposed_shmax_PB(na_grid, na_plate, model = cpm, idp = 2)
head(res)
#>           azi        R
#> [1,] 36.78752 167.7019
#> [2,] 34.10380 178.7247
#> [3,] 31.98870 189.8442
#> [4,] 29.85757 201.7330
#> [5,] 28.37406 213.9229
#> [6,] 27.33043 226.3990
```
