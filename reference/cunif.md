# The Circular Uniform Distribution

Density, probability distribution function, quantiles, and random
generation for the circular uniform distribution.

## Usage

``` r
rcunif(n, axial = FALSE)

dcunif(theta, axial = FALSE, log = FALSE)

pcunif(theta, axial = FALSE, lower.tail = TRUE, log.p = FALSE)

qcunif(p, axial = FALSE, lower.tail = TRUE, log.p = FALSE)
```

## Arguments

- n:

  integer. Number of observations in degrees

- axial:

  logical. Whether the data are axial, i.e. \\\pi\\-periodical (`TRUE`,
  the default) or directional, i.e. \\2 \pi\\-periodical (`FALSE`).

- theta:

  numeric. Angular value in degrees

- log, log.p:

  logical. If `TRUE`, probabilities p are given as log(p).

- lower.tail:

  logical. If `TRUE` (default), probabilities are \\P(\Theta \le
  \theta)\\, otherwise \\P(\Theta \> \theta)\\.

- p:

  numeric. Vector of probabilities with values in \\\[0,1\]\\.

## Value

`dcunif` gives the density, `pcunif` gives the probability of circular
uniform distribution function, `rcunif` generates random deviates (in
degrees), and `qcunif` provides quantiles (in degrees).

## See also

[wcauchy](https://tobiste.github.io/tectonicr/reference/wcauchy.md) and
[vonmises](https://tobiste.github.io/tectonicr/reference/vonmises.md)

## Examples

``` r
set.seed(1)
x <- rcunif(5)

dcunif(x)
#> [1] 0.002777778 0.002777778 0.002777778 0.002777778 0.002777778
dcunif(x, axial = TRUE)
#> [1] 0.005555556 0.005555556 0.005555556 0.005555556 0.005555556

pcunif(x)
#> [1] 0.2655087 0.3721239 0.5728534 0.9082078 0.2016819
qcunif(c(.25, .5, .75))
#> [1]  90 180 270
```
