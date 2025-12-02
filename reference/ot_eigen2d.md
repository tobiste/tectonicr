# Decompositon of orientation tensor

Decompositon of orientation tensor

## Usage

``` r
ot_eigen2d(x, w = NULL, axial = TRUE)
```

## Arguments

- x:

  numeric vector. Values in degrees.

- w:

  (optional) Weights. A vector of positive numbers and of the same
  length as `x`.

- axial:

  logical. Whether the data are axial, i.e. pi-periodical (`TRUE`, the
  default) or directional, i.e. \\2 \pi\\-periodical (`FALSE`).

## Examples

``` r
test <- rvm(100, mean = 0, k = 19)
ot_eigen2d(test, axial = FALSE)
#> Error in ot_eigen2d(test, axial = FALSE): could not find function "ot_eigen2d"

data("nuvel1")
PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
sa.por <- PoR_shmax(san_andreas, PoR, "right")
sa_eig <- ot_eigen2d(sa.por$azi.PoR, w = weighting(san_andreas$unc))
#> Error in ot_eigen2d(sa.por$azi.PoR, w = weighting(san_andreas$unc)): could not find function "ot_eigen2d"
print(sa_eig)
#> Error: object 'sa_eig' not found

rose(sa.por$azi.PoR)

rose_line(sa_eig$vectors, col = c('blue', 'green'))
#> Error: object 'sa_eig' not found
```
