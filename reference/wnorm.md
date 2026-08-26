# The Wrapped Normal Distribution

Density, probability distribution function, quantiles, and random
generation for the circular normal distribution with mean \\\mu\\ and
standard deviation \\\sigma\\.

## Usage

``` r
rwnorm(n, mean = 0, sd = 1)

dwnorm(theta, mean = 0, sd = 1, axial = FALSE, log = FALSE)

pwnorm(theta, mean = 0, sd = 1, axial = FALSE, from = NULL, ...)

qwnorm(
  p,
  mean = 0,
  sd = 1,
  axial = FALSE,
  from = NULL,
  tol = .Machine$double.eps^(0.6),
  ...
)
```

## Arguments

- n:

  number of observations. If `length(n) > 1`, the length is taken to be
  the number required.

- mean:

  numeric. The mean vector in degrees.

- sd:

  numeric. standard deviation of the (unwrapped) normal distribution in
  degrees.

- theta:

  numeric. Angular value in degrees

- axial:

  logical. Whether the data are axial, i.e. \\\pi\\-periodical (`TRUE`,
  the default) or directional, i.e. \\2 \pi\\-periodical (`FALSE`).

- log:

  logical. If `TRUE`, probabilities p are given as \\\log(p)\\.

- from:

  if `NULL` is set to \\\mu-\pi\\. This is the value from which the
  `pvm` and `qvm` are evaluated. in degrees.

- ...:

  optional parameters passed to underlying circular functions
  [`circular::pwrappednormal()`](https://rdrr.io/pkg/circular/man/wrappednormal.html)
  and
  [`circular::qwrappednormal()`](https://rdrr.io/pkg/circular/man/wrappednormal.html)

- p:

  numeric. Vector of probabilities with values in \\\[0,1\]\\.

- tol:

  numeric. The precision in evaluating the distribution function or the
  quantile.

## Value

`dwnorm` gives the density, `pwnorm` gives the probability of the
wrapped normal distribution function, `rwnorm` generates random deviates
(in degrees), and `qwnorm` provides quantiles (in degrees).

## See also

[cunif](https://tobiste.github.io/tectonicr/reference/cunif.md), wnorm,
[wcauchy](https://tobiste.github.io/tectonicr/reference/wcauchy.md), and
[vonmises](https://tobiste.github.io/tectonicr/reference/vonmises.md)

## Examples

``` r
set.seed(1)
x <- rwnorm(5, mean = 90, sd = 5)

dwnorm(x, mean = 90, sd = 5, axial = FALSE)
#> [1] 0.06557252 0.07845431 0.05627449 0.02235206 0.07557240
dwnorm(x, mean = 90, sd = 5, axial = TRUE)
#> [1] 0.06557252 0.07845431 0.05627449 0.02235206 0.07557240

pwnorm(x, mean = 90, sd = 5)
#> [1] 0.2655087 0.5728534 0.2016819 0.9446753 0.6291140
qwnorm(c(.25, .5, .75), mean = 90, sd = 5)
#> [1] 1.511936 1.570796 1.629657
```
