# Weighted Goodness-of-fit Test for Circular Data

Weighted version of the Rayleigh test (or V0-test) for uniformity
against a distribution with a priori expected von Mises concentration.

## Usage

``` r
weighted_rayleigh(x, mu = NULL, w = NULL, axial = TRUE, quiet = FALSE)
```

## Arguments

- x:

  numeric vector. Values in degrees

- mu:

  The *a priori* expected direction (in degrees) for the alternative
  hypothesis.

- w:

  numeric vector weights of length `length(x)`. If `NULL`, the
  non-weighted Rayleigh test is performed.

- axial:

  logical. Whether the data are axial, i.e. \\\pi\\-periodical (`TRUE`,
  the default) or directional, i.e. \\2 \pi\\-periodical (`FALSE`).

- quiet:

  logical. Prints the test's decision.

## Value

a list with the components:

- `R` or `C`:

  mean resultant length or the dispersion (if `mu` is specified). Small
  values of `R` (large values of `C`) will reject uniformity. Negative
  values of `C` indicate that vectors point in opposite directions (also
  lead to rejection).

- `statistic`:

  Test statistic

- `p.value`:

  significance level of the test statistic

## Details

The Null hypothesis is uniformity (randomness). The alternative is a
distribution with a (specified) mean direction (`mu`). If
`statistic >= p.value`, the null hypothesis of randomness is rejected
and angles derive from a distribution with a (or the specified) mean
direction.

## See also

[`rayleigh_test()`](https://tobiste.github.io/tectonicr/reference/rayleigh_test.md)

## Examples

``` r
# Load data
data("cpm_models")
data(san_andreas)
PoR <- equivalent_rotation(cpm_models[["NNR-MORVEL56"]], "na", "pa")
sa.por <- PoR_shmax(san_andreas, PoR, "right")
data("iceland")
PoR.ice <- equivalent_rotation(cpm_models[["NNR-MORVEL56"]], "eu", "na")
ice.por <- PoR_shmax(iceland, PoR.ice, "out")
data("tibet")
PoR.tib <- equivalent_rotation(cpm_models[["NNR-MORVEL56"]], "eu", "in")
tibet.por <- PoR_shmax(tibet, PoR.tib, "in")

# GOF test:
weighted_rayleigh(tibet.por$azi.PoR, mu = 90, w = 1 / tibet$unc)
#> Reject Null Hypothesis
#> $C
#> [1] 0.5321474
#> 
#> $statistic
#> [1] 25.68679
#> 
#> $p.value
#> [1] 2.004315e-143
#> 
weighted_rayleigh(ice.por$azi.PoR, mu = 0, w = 1 / iceland$unc)
#> Reject Null Hypothesis
#> $C
#> [1] 0.3728874
#> 
#> $statistic
#> [1] 11.67322
#> 
#> $p.value
#> [1] 1.826316e-32
#> 
weighted_rayleigh(sa.por$azi.PoR, mu = 135, w = 1 / san_andreas$unc)
#> Reject Null Hypothesis
#> $C
#> [1] 0.8042366
#> 
#> $statistic
#> [1] 38.16524
#> 
#> $p.value
#> [1] 3.589557e-315
#> 
```
