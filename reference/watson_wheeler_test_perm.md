# Watson-Wheeler Test of Homogeneity of Means

A a non-parametric statistical test used to determine whether two or
more independent samples of circular data (angles, directions, or
periodic times) come from the same underlying population distribution.
The difference between the samples can be in either the mean or the
variance.

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

  Significance level of the test. Valid levels are `0.01`, `0.05`, and
  `0.1`. This argument may be omitted (`NULL`, the default), in which
  case, a range for the p-value will be returned.

## Value

list

## Details

### Hypotheses

**Null Hypothesis (H₀)** The samples come from identical populations
(meaning both the mean direction and the dispersion/variance are
homogeneous across groups).

**Alternative Hypothesis (H₁)**: At least one sample comes from a
different population distribution, which can be due to a difference in
the mean direction, a difference in variance/concentration, or both.

### Interpretation

**Test Statistic (W)**: This value follows an approximate chi-squared
(\\\Chi^2\\) distribution. Higher values of W indicate larger
discrepancies between the angular distributions of your groups.

- If the p-value is less than your significance level (commonly α =
  0.05), you **reject** the null hypothesis. This means you have strong
  evidence that the groups differ significantly in their central
  direction or spread around the circle.

- If the p-value is greater than 0.05, you **fail to reject** the null
  hypothesis, meaning there is no statistically significant evidence of
  difference among the groups.

## Note

Important Considerations & Limitations:

- **Sensitivity to both mean and variance**: Because it detects
  differences in either mean or variance, a significant result does not
  automatically mean the mean angles are different; it could be driven
  entirely by differences in concentration (variance).

- **Sample size requirement**: The chi-squared approximation requires
  each group to have a minimum sample size (typically at least 10
  elements per group) to remain valid.

## See also

Other Tests:
[`ar_test()`](https://tobiste.github.io/tectonicr/reference/ar_test.md),
[`kuiper_test()`](https://tobiste.github.io/tectonicr/reference/kuiper_test.md),
[`norm_chisq()`](https://tobiste.github.io/tectonicr/reference/norm_chisq.md),
[`rayleigh-test`](https://tobiste.github.io/tectonicr/reference/rayleigh-test.md),
[`watson_test()`](https://tobiste.github.io/tectonicr/reference/watson_test.md),
[`watson_two_sample`](https://tobiste.github.io/tectonicr/reference/watson_two_sample.md),
[`weighted-rayleigh-test`](https://tobiste.github.io/tectonicr/reference/weighted-rayleigh-test.md)

## Examples

``` r
set.seed(20250411)
x1 <- c(35, 45, 50, 55, 60, 70, 85, 95, 105, 120)
x2 <- c(75, 80, 90, 100, 110, 130, 135, 140, 150, 160, 165)
watson_wheeler_test_perm(x1, x2, axial = FALSE)
#> $statistic
#> [1] 3.67827
#> 
#> $p.value
#> [1] 0.1688312
#> 
#> $alpha
#> NULL
#> 
#> $reject
#> NULL
#> 

data1 <- rvm(n=20, mean = 0, kappa=3)
data2 <- rvm(n=20, mean = 90, kappa=2)
watson_wheeler_test_perm(data1, data2, axial = FALSE)
#> $statistic
#> [1] 15.96367
#> 
#> $p.value
#> [1] 0.000999001
#> 
#> $alpha
#> NULL
#> 
#> $reject
#> NULL
#> 

# San Andreas Fault Data:
data(san_andreas)
data("nuvel1")
PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
sa.por <- PoR_shmax(san_andreas, PoR, "right")
watson_wheeler_test_perm(sa.por$azi.PoR, rvm(100, 135, 10))
#> $statistic
#> [1] 2.934737
#> 
#> $p.value
#> [1] 0.2547453
#> 
#> $alpha
#> NULL
#> 
#> $reject
#> NULL
#> 
```
