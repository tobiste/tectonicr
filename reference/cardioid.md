# The Cardiod Distribution

Density, probability distribution function, quantiles, and random
generation for the Cardiod distribution with mean \\\mu\\ and
concentration \\\rho\\.

## Usage

``` r
rcard(n, mean, rho)

dcard(theta, mean, rho, axial = FALSE)

pcard(theta, mean, rho, lower.tail = TRUE, log.p = FALSE, axial = FALSE)

qcard(
  p,
  mean,
  rho,
  lower.tail = TRUE,
  log.p = FALSE,
  axial = FALSE,
  tol = .Machine$double.eps^0.6
)
```

## Arguments

- n:

  number of observations. If `length(n) > 1`, the length is taken to be
  the number required.

- mean:

  numeric. The mean vector in degrees.

- rho:

  numeric. The concentration parameter \\\rho\\, where \\\|\rho \| \le
  \frac{1}{2}\\

- theta:

  numeric. Angular value in degrees

- axial:

  logical. Whether the data are axial, i.e. \\\pi\\-periodical (`TRUE`,
  the default) or directional, i.e. \\2 \pi\\-periodical (`FALSE`).

- lower.tail:

  logical; if `TRUE` (default), probabilities are \\P\[X \le x\]\\,
  otherwise, \\P\[X \> x\]\\.

- p:

  numeric. Vector of probabilities with values in \\\[0,1\]\\.

- tol:

  numeric. The precision in evaluating the distribution function or the
  quantile.

- log, log.p:

  logical; if `TRUE`, probabilities/densities are given as logarithms.

## Value

`dcard` gives the density, `pcard` gives the probability of the Cardiod
distribution function, `rcard` generates random deviates (in degrees),
and `qcard` provides quantiles (in degrees).

## References

Jeffreys, H., 1948. Theory of Probability, second ed. Oxford University
Press, Oxford

Pewsey, A. (2026). On Jeffreys's cardioid distribution. Computational
Statistics and Data Analysis. https://doi.org/10.1016/j.csda.2025.108248

## See also

[wnorm](https://tobiste.github.io/tectonicr/reference/wnorm.md),
[wcauchy](https://tobiste.github.io/tectonicr/reference/wcauchy.md),
[cunif](https://tobiste.github.io/tectonicr/reference/cunif.md), and
[vonmises](https://tobiste.github.io/tectonicr/reference/vonmises.md)

## Examples

``` r
set.seed(1)
x <- rcard(5, mean = 90, rho = 0.1)

dcard(x, mean = 90, rho = 0.2)
#> [1] 0.2225845 0.2157478 0.1554590 0.1225191 0.2172391
dcard(x, mean = 90, rho = 0.2, axial = TRUE)
#> [1] 0.4437777 0.3922209 0.1918442 0.2753177 0.4029664

pcard(x, mean = 90, rho = 0.2)
#> [1] 0.2946220 0.4185331 0.6364617 0.9140068 0.2204834

qcard(c(.25, .5, .75), mean = 90, rho = 0.2)
#> [1]  73.56566 139.62265 233.42537
```
