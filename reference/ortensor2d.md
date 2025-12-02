# Orientation tensor

Orientation tensor

## Usage

``` r
ortensor2d(x, w = NULL, norm = FALSE)
```

## Arguments

- x:

  numeric vector. Values in degrees.

- w:

  (optional) Weights. A vector of positive numbers and of the same
  length as `x`.

- norm:

  logical. Whether the tensor should be normalized.

## Value

2x2 matrix

## Note

\$\$E = x \cdot x^{T}\$\$

## See also

[`ot_eigen2d()`](https://tobiste.github.io/tectonicr/reference/ot_eigen2d.md)

## Examples

``` r
test <- rvm(100, mean = 0, k = 10)
ortensor2d(test)
#>            [,1]       [,2]
#> [1,] 0.90471294 0.05095728
#> [2,] 0.05095728 0.09528706
```
