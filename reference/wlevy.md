# Wrapped Lévy Distribution

Probability density function of the wrapped Lévy distribution

## Usage

``` r
dwlevy(theta, mean, c, axial = FALSE, tol = 1e-08, P0 = 100, P.max = 5000)
```

## Arguments

- theta:

  numeric. Angular value in degrees

- mean:

  numeric. The mean vector in degrees.

- c:

  numeric. Scale factor of the Lévy distribution. Small values indicate
  concentrated data near the mean.

- axial:

  logical. Whether the data are axial, i.e. \\\pi\\-periodical (`TRUE`,
  the default) or directional, i.e. \\2 \pi\\-periodical (`FALSE`).

- tol:

  numeric. The precision in evaluating the distribution function or the
  quantile.

- P0, P.max:

  numeric. Arguments Control how the truncated series (the
  Fourier/characteristic-function sum) is grown and when it gives up. In
  practice you'll almost never touch either one, as for typical `c` (\\c
  \le 0.1\\) convergence happens within the first one or two doublings,
  well under `P0`.

## See also

[cunif](https://tobiste.github.io/tectonicr/reference/cunif.md),
[wnorm](https://tobiste.github.io/tectonicr/reference/wnorm.md),
[wcauchy](https://tobiste.github.io/tectonicr/reference/wcauchy.md), and
[vonmises](https://tobiste.github.io/tectonicr/reference/vonmises.md)

## Examples

``` r
set.seed(1)
x <- rwcauchy(5, mean = 90, rho = exp(-1))

dwlevy(x, mean = 80, c = 10)
#> Error in dwlevy(x, mean = 80, c = 10): could not find function "dwlevy"
```
