# Watson's \\U^2\\ Test of Circular Uniformity

Watson's test statistic is a rotation-invariant Cramer - von Mises test.
non-parametric, rank-based alternative to one-sample

## Usage

``` r
watson_test(
  x,
  alpha = NULL,
  dist = c("uniform", "vonmises"),
  axial = TRUE,
  quiet = FALSE
)
```

## Arguments

- x:

  numeric vector. Values in degrees

- alpha:

  Significance level of the test. Valid levels are `0.01`, `0.05`, and
  `0.1`. This argument may be omitted (`NULL`, the default), in which
  case, a range for the p-value will be returned.

- dist:

  Distribution to test for. The default, `"uniform"`, is the uniform
  distribution. `"vonmises"` tests the von Mises distribution.

- axial:

  logical. Whether the data are axial, i.e. \\\pi\\-periodical (`TRUE`,
  the default) or circular, i.e. \\2 \pi\\-periodical (`FALSE`).

- quiet:

  logical. Prints the test's decision.

## Value

list containing the test statistic `statistic` and the significance
level `p.value` (the critical value).

## Details

If `statistic > p.value`, the null hypothesis is rejected. If not,
randomness (uniform distribution) cannot be excluded.

## References

Mardia and Jupp (1999). Directional Statistics. John Wiley and Sons.

## See also

Other Tests:
[`ar_test()`](https://tobiste.github.io/tectonicr/reference/ar_test.md),
[`kuiper_test()`](https://tobiste.github.io/tectonicr/reference/kuiper_test.md),
[`norm_chisq()`](https://tobiste.github.io/tectonicr/reference/norm_chisq.md),
[`rayleigh-test`](https://tobiste.github.io/tectonicr/reference/rayleigh-test.md),
[`watson_two_sample`](https://tobiste.github.io/tectonicr/reference/watson_two_sample.md),
[`watson_wheeler_test_perm()`](https://tobiste.github.io/tectonicr/reference/watson_wheeler_test_perm.md),
[`weighted-rayleigh-test`](https://tobiste.github.io/tectonicr/reference/weighted-rayleigh-test.md)

## Examples

``` r
# Example data from Mardia and Jupp (1999), pp. 93
pidgeon_homing <- c(55, 60, 65, 95, 100, 110, 260, 275, 285, 295)
watson_test(pidgeon_homing, axial = FALSE, alpha = .05)
#> Do Not Reject Null Hypothesis
#> $statistic
#> [1] 0.1153633
#> 
#> $p.value
#> [1] 0.187
#> 

# San Andreas Fault Data:
data(san_andreas)
data("nuvel1")
PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
sa.por <- PoR_shmax(san_andreas, PoR, "right")
watson_test(sa.por$azi.PoR, alpha = .05)
#> Reject Null Hypothesis
#> $statistic
#> [1] 31.68687
#> 
#> $p.value
#> [1] 0.187
#> 
watson_test(sa.por$azi.PoR, alpha = .05, dist = "vonmises")
#> Reject Null Hypothesis
#> $statistic
#> [1] 77.28103
#> 
#> $p.value
#> [1] 0.101
#> 
watson_test(sa.por$azi.PoR, alpha = .05, dist = "vonmises")
#> Reject Null Hypothesis
#> $statistic
#> [1] 77.28103
#> 
#> $p.value
#> [1] 0.101
#> 
```
