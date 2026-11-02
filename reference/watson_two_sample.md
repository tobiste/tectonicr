# Watson's Two-Sample Test of Homogeneity

Performs Watson's test for homogeneity on two samples of circular data.
`watson_two_test_perm()` uses permutation to estimate p-values.

## Usage

``` r
watson_two_test(x, y, alpha = 0.05, axial = TRUE, quiet = FALSE)

watson_two_test_perm(x, y, axial = TRUE, n_perm = 1000L, alpha = NULL)
```

## Arguments

- x, y:

  numeric vectors. Angles in degrees

- alpha:

  Significance level of the test ranging between 0 and 1. For
  `dist = "vonmises"` it must be one of `0.1`, `0.05`, `0.01` or `0`
  (only the tabulated p-value range is returned).

- axial:

  logical. Whether the data are axial, i.e. \\\pi\\-periodical (`TRUE`,
  the default) or directional, i.e. \\2 \pi\\-periodical (`FALSE`). In
  case of axial data, the angles will be doubled for the test.

- quiet:

  logical. Prints the test's decision.

- n_perm:

  integer. Number of permutations

## Value

A list with class `"htest"`

## Details

A two-sample Watson's \\U^2\\ test determines whether two independent
groups of circular data (angles or directions) come from the same
underlying distribution. Ties (e.g. azimuths rounded to whole degrees)
are handled by evaluating the difference of the empirical distribution
functions only at distinct values, weighted by the number of tied
observations.

### Hypotheses

**Null Hypothesis** (\\H_0\\): The two samples come from the same
circular distribution (the two groups of angles share a common
distribution around the circle).

**Alternative Hypothesis** (\\H_1\\): The two samples come from
different circular distributions (the groups are oriented or dispersed
differently around the circle).

### Interpretation

The Test Statistic (\\U^2\\) measures the distance between the
cumulative distribution functions of the two circular samples. A larger
\\U^2\\ value means the two sets of angles look more different from each
other.

The p-value of `watson_two_test_perm()` represents the probability of
getting a U² value as large as (or larger than) your observed value
purely by chance, assuming the null hypothesis is true. It is calculated
by shuffling the group labels across your data many times to build a
reference permutation distribution.

#### Making a Decision

- Low p-value (\\p \le \alpha\\, usually 0.05): Reject the null
  hypothesis. Conclude that the two groups have significantly different
  circular distributions.

- High p-value (\\p \> \alpha\\): Fail to reject the null hypothesis.
  There is not enough evidence to say the two groups are distributed
  differently around the circle.

## Note

Critical values and p-values of `watson_two_test()` are obtained using
the asymptotic distribution of the test statistic. They should only be
used for combined sample sizes in excess of 17. Tables are available for
smaller sample sizes and can be found in Mardia (1972) for instance;
alternatively use `watson_two_test_perm()`.

## See also

Other Two-Sample-Tests:
[`ar_test()`](https://tobiste.github.io/tectonicr/reference/ar_test.md),
[`equal_kappa_test()`](https://tobiste.github.io/tectonicr/reference/equal_kappa_test.md),
[`uniform_scores_test()`](https://tobiste.github.io/tectonicr/reference/uniform_scores_test.md),
[`walraff_test()`](https://tobiste.github.io/tectonicr/reference/walraff_test.md),
[`watson_wheeler_test_perm()`](https://tobiste.github.io/tectonicr/reference/watson_wheeler_test_perm.md),
[`watson_williams_test()`](https://tobiste.github.io/tectonicr/reference/watson_williams_test.md)

## Examples

``` r
set.seed(20250411)
x1 <- c(35, 45, 50, 55, 60, 70, 85, 95, 105, 120)
x2 <- c(75, 80, 90, 100, 110, 130, 135, 140, 150, 160, 165)
watson_two_test(x1, x2, axial = FALSE)
#> 
#>  Watson two-sample U2 test
#> 
#> data:  x1 and x2
#> U2 = 0.14574, n = 21, p-value = 0.1126
#> alternative hypothesis: Samples come from different distributions
#> 
watson_two_test_perm(x1, x2, axial = FALSE)
#> 
#>  Watson two-sample U2 test based on 1000 permutations
#> 
#> data:  x1 and x2
#> U2 = 0.14574, n = 21, n_perm = 1000, p-value = 0.1199
#> alternative hypothesis: Samples come from different distributions
#> 

data1 <- rvm(n=20, mean = 0, kappa=3)
data2 <- rvm(n=20, mean = 90, kappa=2)
watson_two_test(data1, data2, axial = FALSE)
#> 
#>  Watson two-sample U2 test
#> 
#> data:  data1 and data2
#> U2 = 0.61775, n = 40, p-value = 1.012e-05
#> alternative hypothesis: Samples come from different distributions
#> 
watson_two_test_perm(data1, data2, axial = FALSE)
#> 
#>  Watson two-sample U2 test based on 1000 permutations
#> 
#> data:  data1 and data2
#> U2 = 0.61775, n = 40, n_perm = 1000, p-value = 0.000999
#> alternative hypothesis: Samples come from different distributions
#> 


# San Andreas Fault Data:
data(san_andreas)
data("nuvel1")
PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
sa.por <- PoR_shmax(san_andreas, PoR, "right")
watson_two_test(sa.por$azi, rvm(100, 135, 10), alpha = 0.05)
#> 
#>  Watson two-sample U2 test
#> 
#> data:  sa.por$azi and rvm(100, 135, 10)
#> U2 = 0.20134, n = 1226, p-value = 0.03759
#> alternative hypothesis: Samples come from different distributions
#> 
watson_two_test_perm(sa.por$azi, rvm(100, 135, 10), alpha = 0.05)
#> 
#>  Watson two-sample U2 test based on 1000 permutations
#> 
#> data:  sa.por$azi and rvm(100, 135, 10)
#> U2 = 0.21981, n = 1226, n_perm = 1000, p-value = 0.02498
#> alternative hypothesis: Samples come from different distributions
#> 
```
