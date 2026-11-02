# Watson-Wheeler Test of Homogeneity

A non-parametric statistical test used to determine whether two
independent samples of circular data (angles, directions, or periodic
times) come from the same underlying population distribution. The
difference between the samples can be in either the mean direction or
the dispersion. P-values are estimated using permutation.

## Usage

``` r
watson_wheeler_test_perm(x, y, axial = TRUE, n_perm = 1000L, alpha = NULL)
```

## Arguments

- x, y:

  numeric vectors. Angles in degrees

- axial:

  logical. Whether the data are axial, i.e. \\\pi\\-periodical (`TRUE`,
  the default) or directional, i.e. \\2 \pi\\-periodical (`FALSE`). In
  case of axial data, the angles will be doubled for the test.

- n_perm:

  integer. Number of permutations

- alpha:

  Significance level of the test ranging between 0 and 1. For
  `dist = "vonmises"` it must be one of `0.1`, `0.05`, `0.01` or `0`
  (only the tabulated p-value range is returned).

## Value

A list with class `"htest"`

## Details

The statistic is \\R^\* = 2 (n - 1) (C_1^2 + S_1^2) / (n_1 n_2)\\ on the
uniform scores \\2 \pi r_i / n\\ of the ranks \\r_i\\ in the pooled
sample (Mardia & Jupp, 2000, eq. 8.3.4). Tied values get their average
rank.

### Hypotheses

**Null Hypothesis (\\H_0\\)** The samples come from identical
populations (meaning both the mean direction and the dispersion/variance
are homogeneous across groups).

**Alternative Hypothesis (\\H_1\\)**: The samples come from different
population distributions, which can be due to a difference in the mean
direction, a difference in variance/concentration, or both.

### Interpretation

**Test Statistic**: Higher values indicate larger discrepancies between
the angular distributions of your groups.

- If the p-value is less than your significance level (commonly
  \\\alpha\\ = 0.05), you **reject** the null hypothesis. This means you
  have strong evidence that the groups differ significantly in their
  central direction or spread around the circle.

- If the p-value is greater than 0.05, you **fail to reject** the null
  hypothesis, meaning there is no statistically significant evidence of
  difference among the groups.

## Note

Important Considerations & Limitations:

- **Sensitivity to both mean and variance**: Because it detects
  differences in either mean or variance, a significant result does not
  automatically mean the mean angles are different; it could be driven
  entirely by differences in concentration (variance).

- For more than two samples, or the asymptotic \\\chi^2\\ version, see
  [`uniform_scores_test()`](https://tobiste.github.io/tectonicr/reference/uniform_scores_test.md).

## See also

Other Two-Sample-Tests:
[`ar_test()`](https://tobiste.github.io/tectonicr/reference/ar_test.md),
[`equal_kappa_test()`](https://tobiste.github.io/tectonicr/reference/equal_kappa_test.md),
[`uniform_scores_test()`](https://tobiste.github.io/tectonicr/reference/uniform_scores_test.md),
[`walraff_test()`](https://tobiste.github.io/tectonicr/reference/walraff_test.md),
[`watson_two_sample`](https://tobiste.github.io/tectonicr/reference/watson_two_sample.md),
[`watson_williams_test()`](https://tobiste.github.io/tectonicr/reference/watson_williams_test.md)

## Examples

``` r
set.seed(20250411)
x1 <- c(35, 45, 50, 55, 60, 70, 85, 95, 105, 120)
x2 <- c(75, 80, 90, 100, 110, 130, 135, 140, 150, 160, 165)
watson_wheeler_test_perm(x1, x2, axial = FALSE)
#> 
#>  Watson-Wheeler test of homogeneity based on 1000 permutations
#> 
#> data:  x1 and x2
#> R* = 3.6783, n = 21, n_perm = 1000, p-value = 0.1489
#> alternative hypothesis: Samples come from different distributions
#> 

data1 <- rvm(n=20, mean = 0, kappa=3)
data2 <- rvm(n=20, mean = 90, kappa=2)
watson_wheeler_test_perm(data1, data2, axial = FALSE)
#> 
#>  Watson-Wheeler test of homogeneity based on 1000 permutations
#> 
#> data:  data1 and data2
#> R* = 28.652, n = 40, n_perm = 1000, p-value = 0.000999
#> alternative hypothesis: Samples come from different distributions
#> 

# San Andreas Fault Data:
data(san_andreas)
data("nuvel1")
PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
sa.por <- PoR_shmax(san_andreas, PoR, "right")
watson_wheeler_test_perm(sa.por$azi.PoR, rvm(100, 135, 10))
#> 
#>  Watson-Wheeler test of homogeneity based on 1000 permutations
#> 
#> data:  sa.por$azi.PoR and rvm(100, 135, 10)
#> R* = 1.9205, n = 1226, n_perm = 1000, p-value = 0.3836
#> alternative hypothesis: Samples come from different distributions
#> 
```
