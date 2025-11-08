# The von Mises Distribution

Density, probability distribution function, quantiles, and random
generation for the circular normal distribution with mean and kappa.

## Usage

``` r
rvm(n, mean, kappa)

dvm(theta, mean, kappa, log = FALSE, axial = FALSE)

pvm(theta, mean, kappa, from = NULL, tol = 1e-20)

qvm(p, mean = 0, kappa, from = NULL, tol = .Machine$double.eps^(0.6), ...)
```

## Arguments

- n:

  integer. Number of observations in degrees

- mean:

  numeric. Mean angle in degrees

- kappa:

  numeric. Concentration parameter in the range (0, Inf\]

- theta:

  numeric. Angular value in degrees

- log:

  logical. If `TRUE`, probabilities p are given as log(p).

- axial:

  logical. Whether the data are axial, i.e. \\\pi\\-periodical (`TRUE`,
  the default) or directional, i.e. \\2 \pi\\-periodical (`FALSE`).

- from:

  if `NULL` is set to \\\mu-\pi\\. This is the value from which the pvm
  and qvm are evaluated. in degrees.

- tol:

  numeric. The precision in evaluating the distribution function or the
  quantile.

- p:

  numeric. Vector of probabilities with values in \\\[0,1\]\\.

- ...:

  parameters passed to
  [`stats::integrate()`](https://rdrr.io/r/stats/integrate.html).

## Value

`dvm` gives the density, `pvm` gives the probability of the von Mises
distribution function, `rvm` generates random deviates (in degrees), and
`qvm` provides quantiles (in degrees).

## Examples

``` r
set.seed(1)
x <- rvm(5, mean = 90, kappa = 2)

dvm(x, mean = 90, kappa = 2)
#> [1] 0.46942367 0.01695767 0.21318638 0.49589993 0.08107754
dvm(x, mean = 90, kappa = 2, axial = TRUE)
#> [1] 0.71367193 0.14001117 0.06570123 0.88231652 0.01932479

pvm(x, mean = 90, kappa = 2)
#> [1] 0.6542335 0.9908071 0.1148932 0.3986252 0.9568411
qvm(c(.25, .5, .75), mean = 90, kappa = 2)
#> [1]  59.65254  90.00000 120.34746
```
