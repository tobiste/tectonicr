# Rayleigh Test of Circular Uniformity

Performs a Rayleigh test for uniformity of circular/directional data by
assessing the significance of the mean resultant length.

## Usage

``` r
rayleigh_test(x, mu = NULL, axial = TRUE, quiet = FALSE)
```

## Arguments

- x:

  numeric vector. Values in degrees

- mu:

  (optional) The specified or known mean direction (in degrees) in
  alternative hypothesis

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

  test statistic

- `p.value`:

  significance level of the test statistic

## Details

- \\H_0\\::

  angles are randomly distributed around the circle.

- \\H_1\\::

  angles are from non-uniformly distribution with unknown mean direction
  and mean resultant length (when `mu` is `NULL`. Alternatively (when
  `mu` is specified), angles are non-uniformly distributed around a
  specified direction.

If `statistic > p.value`, the null hypothesis is rejected, i.e. the
length of the mean resultant differs significantly from zero, and the
angles are not randomly distributed.

## Note

Although the Rayleigh test is consistent against (non-uniform) von Mises
alternatives, it is not consistent against alternatives with `p = 0` (in
particular, distributions with antipodal symmetry, i.e. axial data).
Tests of non-uniformity which are consistent against all alternatives
include Kuiper's test
([`kuiper_test()`](https://tobiste.github.io/tectonicr/reference/kuiper_test.md))
and Watson's \\U^2\\ test
([`watson_test()`](https://tobiste.github.io/tectonicr/reference/watson_test.md)).

## References

Fisher, N. I. (1993) Statistical Analysis of Circular Data, Cambridge
University Press.

## See also

[`mean_resultant_length()`](https://tobiste.github.io/tectonicr/reference/mean_resultant_length.md),
[`circular_mean()`](https://tobiste.github.io/tectonicr/reference/circle_stats.md),
[`norm_chisq()`](https://tobiste.github.io/tectonicr/reference/norm_chisq.md),
[`kuiper_test()`](https://tobiste.github.io/tectonicr/reference/kuiper_test.md),
[`watson_test()`](https://tobiste.github.io/tectonicr/reference/watson_test.md),
[`weighted_rayleigh()`](https://tobiste.github.io/tectonicr/reference/weighted_rayleigh.md)

## Examples

``` r
# Example data from Mardia and Jupp (1999), pp. 93
pidgeon_homing <- c(55, 60, 65, 95, 100, 110, 260, 275, 285, 295)
rayleigh_test(pidgeon_homing, axial = FALSE) # Do not reject null hypothesis.
#> Do Not Reject Null Hypothesis
#> $R
#> [1] 0.2228717
#> 
#> $statistic
#> [1] 0.4967179
#> 
#> $p.value
#> [1] 0.6201354
#> 
# R = 0.22; stat = 0.497, p = 0.62

# Example data from Davis (1986), pp. 316
finland_striae <- c(
  23, 27, 53, 58, 64, 83, 85, 88, 93, 99, 100, 105, 113,
  113, 114, 117, 121, 123, 125, 126, 126, 126, 127, 127, 128, 128, 129, 132,
  132, 132, 134, 135, 137, 144, 145, 145, 146, 153, 155, 155, 155, 157, 163,
  165, 171, 172, 179, 181, 186, 190, 212
)
rayleigh_test(finland_striae, axial = FALSE) # reject null hypothesis
#> Reject Null Hypothesis
#> $R
#> [1] 0.8003694
#> 
#> $statistic
#> [1] 32.67015
#> 
#> $p.value
#> [1] 6.479397e-15
#> 
rayleigh_test(finland_striae, mu = 105, axial = FALSE) # reject null hypothesis
#> Reject Null Hypothesis
#> $C
#> [1] 0.7300887
#> 
#> $statistic
#> [1] 7.373534
#> 
#> $p.value
#> [1] 2.130845e-13
#> 

# Example data from Mardia and Jupp (1999), pp. 99
atomic_weight <- c(
  rep(0, 12), rep(3.6, 1), rep(36, 6), rep(72, 1),
  rep(108, 2), rep(169.2, 1), rep(324, 1)
)
rayleigh_test(atomic_weight, 0, axial = FALSE) # reject null hypothesis
#> Reject Null Hypothesis
#> $C
#> [1] 0.7237434
#> 
#> $statistic
#> [1] 5.014241
#> 
#> $p.value
#> [1] 5.348331e-08
#> 

# San Andreas Fault Data:
data(san_andreas)
rayleigh_test(san_andreas$azi) # reject null hypothesis
#> Reject Null Hypothesis
#> $R
#> [1] 0.6822906
#> 
#> $statistic
#> [1] 524.1761
#> 
#> $p.value
#> [1] 2.255452e-228
#> 
data("nuvel1")
PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
sa.por <- PoR_shmax(san_andreas, PoR, "right")
rayleigh_test(sa.por$azi.PoR, mu = 135) # reject null hypothesis
#> Reject Null Hypothesis
#> $C
#> [1] 0.7009544
#> 
#> $statistic
#> [1] 33.26396
#> 
#> $p.value
#> [1] 1.420092e-239
#> 
```
