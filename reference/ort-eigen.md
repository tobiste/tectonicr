# Decomposition of Orientation Tensor

Eigenvector decomposition of the orientation tensor

## Usage

``` r
ot_eigen2d(x, w = NULL, axial = TRUE, scale = FALSE)

principal_direction(x, w = NULL, axial = TRUE)

orientation_strength(x, w = NULL, axial = TRUE)
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
  Only applicable when weighting are specified, otherwise the
  eigenvalues are always scaled.

## Value

`ot_eigen2d` returns a list of the Eigenvalues and the angles
corresponding to the Eigenvectors. `principal_direction()` and
`orientation_strength()` are convenience functions to return the
orientation of the largest eigenvalue, and the orientation strength,
respectively.

## Details

The two perpendicular **Eigenvectors** are the "principal directions"
towards the highest and the lowest concentration of orientation data.

The **Eigenvalues** can be interpreted as the fractions of the data
explained by the orientation of the associated principal direction.
Thus, the strength of the orientation is the largest eigenvalue
normalized by the sum of the eigenvalues (`scale=TRUE`).

## Note

Eigenvalues and eigenvectors of the orientation tensor (inertia tensor)
are also called "principle moments of inertia" and "principle axes of
inertia", respectively.

## See also

[`ortensor2d()`](https://tobiste.github.io/tectonicr/reference/ortensor2d.md)

## Examples

``` r
test <- rvm(100, mean = 0, k = 10)
ot_eigen2d(test, axial = FALSE)
#> eigen() decomposition
#> $values
#> [1] 0.90790834 0.09209166
#> 
#> $vectors
#> [1]  3.588165 93.588165
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
rose_line(sa_eig$vectors,
  col = c("red", "green"),
  radius = sa_eig$values, lwd = 2
)
graphics::legend("topright",
  legend = round(sa_eig$values, 2),
  col = c("red", "green"), lty = 1
)


principal_direction(sa.por$azi.PoR)
#> [1] 139.8762

orientation_strength(sa.por$azi.PoR)
#> [1] 0.6904934
```
