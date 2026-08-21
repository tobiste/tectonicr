# The Wrapped Cauchy Distribution

Density, probability distribution function, quantiles, and random
generation for the circular wrapped Cauchy distribution with mean
\\\mu\\ and rho \\\rho\\

## Usage

``` r
rwcauchy(n, mean, rho)

dwcauchy(theta, mean, rho, axial = FALSE, log = FALSE)

pwcauchy(
  theta,
  mean,
  rho,
  axial = FALSE,
  from = NULL,
  lower.tail = TRUE,
  log.p = FALSE
)

qwcauchy(
  p,
  mean,
  rho,
  axial = FALSE,
  from = NULL,
  lower.tail = TRUE,
  log.p = FALSE
)
```

## Arguments

- n:

  integer. Number of observations in degrees

- mean:

  numeric. The mean vector in degrees.

- rho:

  numeric. Concentration parameter in the range (0, 1)

- theta:

  numeric. Angular value in degrees

- axial:

  logical. Whether the data are axial, i.e. \\\pi\\-periodical (`TRUE`,
  the default) or directional, i.e. \\2 \pi\\-periodical (`FALSE`).

- log:

  logical. If `TRUE`, probabilities p are given as log(p).

- from:

  if `NULL` is set to \\\mu-\pi\\. This is the value from which the pvm
  and qvm are evaluated. in degrees.

- lower.tail:

  logical. If `TRUE` (default), probabilities are \\P(\Theta \le
  \theta)\\, otherwise \\P(\Theta \> \theta)\\.

- log.p:

  logical. If `TRUE`, probabilities p are given as log(p).

- p:

  numeric. Vector of probabilities with values in \\\[0,1\]\\.

## Value

`dwcauchy` gives the density, `pwcauchy` gives the probability of the
wrapped Cauchy distribution function, `rwcauchy` generates random
deviates (in degrees), and `qrwcauchy` provides quantiles (in degrees).

## See also

[vonmises](https://tobiste.github.io/tectonicr/reference/vonmises.md)
and [cunif](https://tobiste.github.io/tectonicr/reference/cunif.md)

## Examples

``` r
set.seed(1)
x <- rwcauchy(5, mean = 90, rho = exp(-1))

dwcauchy(x, mean = 90, rho = exp(-1))
#> [1] 0.002990156 0.001451825 0.001673540 0.005563540 0.004075399
dwcauchy(x, mean = 90, rho = exp(-1), axial = TRUE)
#> [1] 0.003057101 0.004218155 0.002953053 0.009144468 0.004528445

pwcauchy(x, mean = 90, rho = exp(-1))
#> [1] 0.7948404 0.9396023 0.9072813 0.4004565 0.7210158
qwcauchy(c(.25, .5, .75), mean = 90, rho = exp(-1))
#> [1]  40.39506  90.00000 139.60494
```
