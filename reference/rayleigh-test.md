# Rayleigh Test of Circular Uniformity

A test to determine whether a sample of circular or directional data is
evenly spread out or clustered around a single specific direction. The
test assesses the significance of the mean resultant length.
`rayleigh_test_perm()` estimates p-values by Monte Carlo simulation of
uniform samples.

## Usage

``` r
rayleigh_test(x, mu = NULL, axial = TRUE, alpha = 0.05, quiet = FALSE)

rayleigh_test_perm(x, mu = NULL, axial = TRUE, n_perm = 1000L)
```

## Arguments

- x:

  numeric vector. Values in degrees

- mu:

  (optional) The specified or known mean direction (in degrees) in
  alternative hypothesis.

- axial:

  logical. Whether the data are axial, i.e. \\\pi\\-periodical (`TRUE`,
  the default) or directional, i.e. \\2 \pi\\-periodical (`FALSE`). In
  case of axial data, the angles will be doubled for the test.

- alpha:

  Significance level of the test ranging between 0 and 1.

- quiet:

  logical. Prints the test's decision.

- n_perm:

  integer. Number of simulated uniform samples.

## Value

a list with the components:

- `R` or `C`:

  mean resultant length or, if `mu` is specified, the mean cosine of the
  deviations from `mu`. Large values of `R` (large positive values of
  `C`) reject uniformity. Negative values of `C` indicate that the
  vectors point away from `mu`; the one-sided V-test does not reject
  them.

- `statistic`:

  test statistic

- `p.value`:

  significance level of the test statistic

## Details

### Hypotheses

**Null Hypothesis** \\H_0\\: The population is distributed uniformly
(randomly) around the circle with no preferred direction.

**Alternative Hypothesis** \\H_1\\: The population is not uniform and
has a unimodal (single-peaked) concentration in a preferred direction.
When `mu` is specified (V-test), the alternative is a unimodal
distribution centred on the specified direction; the test is one-sided.

*Mean Resultant Length* (\\\bar{R}\\ or R): A value between 0 and 1 that
measures how concentrated the data points are.

- \\\bar{R}\\ = 0: The data is completely spread out around the circle

- \\\bar{R}\\ = 1: All data points point in the exact same direction.

*p-value* The probability of seeing data this clustered purely by chance
under the assumption of uniformity.

### Interpretation

- Small p-value (p \< 0.05): Reject the null hypothesis. The length of
  the mean resultant differs significantly from zero, and the angles are
  not randomly distributed. You have strong evidence that the data
  points share a significant preferred or mean direction (unimodal
  clustering).

- Large p-value (p \\\ge\\ 0.05): Fail to reject the null hypothesis.
  There is not enough evidence to claim a preferred direction, meaning
  the data looks random or uniform around the circle.

## Note

Although the Rayleigh test is consistent against (non-uniform) von Mises
alternatives, it is not consistent against alternatives with a mean
resultant length \\\rho = 0\\ (in particular, distributions with
antipodal symmetry, i.e. undoubled axial data). Tests of non-uniformity
which are consistent against all alternatives include Kuiper's test
([`kuiper_test()`](https://tobiste.github.io/tectonicr/reference/kuiper_test.md))
and Watson's \\U^2\\ test
([`watson_test()`](https://tobiste.github.io/tectonicr/reference/watson_test.md)).

### Limitations

- The test assumes a unimodal alternative (one main peak).

- If your data has two opposite clusters (bimodal or axial data, like a
  bi-directional line trend), the Rayleigh test can yield a
  high/non-significant p-value because the opposing vectors cancel each
  other out.

## References

Fisher, N. I. (1993) Statistical Analysis of Circular Data, Cambridge
University Press.

## See also

[`mean_resultant_length()`](https://tobiste.github.io/tectonicr/reference/mean_resultant_length.md),
[`circular_mean()`](https://tobiste.github.io/tectonicr/reference/circle_stats.md)

Other One-Sample-Tests:
[`hodges_ajne_test()`](https://tobiste.github.io/tectonicr/reference/hodges_ajne_test.md),
[`kuiper_test()`](https://tobiste.github.io/tectonicr/reference/kuiper_test.md),
[`norm_chisq()`](https://tobiste.github.io/tectonicr/reference/norm_chisq.md),
[`rao_spacing_test()`](https://tobiste.github.io/tectonicr/reference/rao_spacing_test.md),
[`watson_test()`](https://tobiste.github.io/tectonicr/reference/watson_test.md),
[`weighted-rayleigh-test`](https://tobiste.github.io/tectonicr/reference/weighted-rayleigh-test.md)

## Examples

``` r
# Example data from Mardia and Jupp (1999), pp. 93
rayleigh_test(homing, axial = FALSE) # Do not reject null hypothesis.
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
rayleigh_test_perm(homing, axial = FALSE)
#> 
#>  Rayleigh test of uniformity (Monte Carlo p-value based on 1000 uniform
#>  samples)
#> 
#> data:  homing
#> Rbar = 0.22287, n = 10, n_perm = 1000, p-value = 0.6653
#> alternative hypothesis: unimodal distribution
#> sample estimates:
#>         R      Rbar 
#> 2.2287169 0.2228717 
#> 

# Example data from Davis (1986), pp. 316
rayleigh_test(striae, axial = FALSE) # reject null hypothesis
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
rayleigh_test_perm(striae, axial = FALSE) # reject null hypothesis
#> 
#>  Rayleigh test of uniformity (Monte Carlo p-value based on 1000 uniform
#>  samples)
#> 
#> data:  striae
#> Rbar = 0.80037, n = 51, n_perm = 1000, p-value = 0.000999
#> alternative hypothesis: unimodal distribution
#> sample estimates:
#>          R       Rbar 
#> 40.8188388  0.8003694 
#> 

rayleigh_test(striae, mu = 105, axial = FALSE) # reject null hypothesis
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
rayleigh_test_perm(striae, mu = 105, axial = FALSE) # reject null hypothesis
#> 
#>  Rayleigh V-test of uniformity (Monte Carlo p-value based on 1000
#>  uniform samples)
#> 
#> data:  striae
#> C = 0.73009, n = 51, n_perm = 1000, p-value = 0.000999
#> alternative hypothesis: unimodal distribution with mean direction 105
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
rayleigh_test_perm(sa.por$azi.PoR, mu = 135, n_perm = 1e3) # reject null hypothesis
#> 
#>  Rayleigh V-test of uniformity (Monte Carlo p-value based on 1000
#>  uniform samples)
#> 
#> data:  sa.por$azi.PoR
#> C = 0.70095, n = 1126, n_perm = 1000, p-value = 0.000999
#> alternative hypothesis: unimodal distribution with mean direction 135
#> 
```
