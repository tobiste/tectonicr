# Decomposition of orientation tensor

Eigenvector decomposition of the orientations.

## Usage

``` r
ot_eigen2d(x, w = NULL, axial = TRUE, scale = FALSE)
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

- scale:

  logical. Whether the Eigenvalues should be scaled so they sum up to 1.

## Value

list of Eigenvalues and the angles corresponding to the Eigenvectors.

## Details

Eigenvalues can be interpreted as the fraction of the data explained by
the orientation of the associated Eigenvector.

## See also

[`ortensor2d()`](https://tobiste.github.io/tectonicr/reference/ortensor2d.md)

## Examples

``` r
test <- rvm(100, mean = 0, k = 10)
ot_eigen2d(test, axial = FALSE)
#> eigen() decomposition
#> $values
#> [1] 0.90093653 0.09906347
#> 
#> $vectors
#> [1]  1.008136 91.008136
#> 

data("nuvel1")
PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
sa.por <- PoR_shmax(san_andreas, PoR, "right")
sa_eig <- ot_eigen2d(sa.por$azi.PoR, w = weighting(san_andreas$unc), scale = TRUE)
print(sa_eig)
#> eigen() decomposition
#> $values
#> [1] 0.7051428 0.2948572
#> 
#> $vectors
#> [1] 139.17907  49.17907
#> 

rose(sa.por$azi.PoR, muci = FALSE)
rose_line(sa_eig$vectors, col = c('red', 'green'),
  radius = sa_eig$values, lwd = 2)
graphics::legend("topright",
  legend = round(sa_eig$values, 2),
  col = c('red', 'green'), lty = 1)
```
