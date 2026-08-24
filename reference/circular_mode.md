# Circular Mode

MLE angle (maximum density) using a circular distribution kernel with
specified concentration

## Usage

``` r
circular_mode(x, ...)
```

## Arguments

- x:

  numeric. A vector of angles (in degrees) from which the estimate is to
  be computed.

- ...:

  parameters passed to
  [`circular_density()`](https://tobiste.github.io/tectonicr/reference/circular_density.md)

## Value

numeric

## References

N.I. Fisher (1993) Statistical Analysis of Circular Data, Cambridge
University Press.

## Examples

``` r
set.seed(20250411)
x <- rvm(10, 0, 100)

# Mode of von Mises kernel density (the default)
circular_mode(x)
#> [1] 176.8297

# Mode of wrapped Cauchy kernel density
circular_mode(x, kernel = "wrappedcauchy")
#> [1] 176.1252
```
