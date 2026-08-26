# Watson's Two-Sample Test of Homogeneity

Performs Watson's test for homogeneity on two samples of circular data.
`watson_two_test_perm()` uses permutation to estimate p-values.

## Usage

``` r
watson_two_test(x, y, alpha = NULL, axial = TRUE, quiet = FALSE)

watson_two_test_perm(x, y, axial = TRUE, n_perm = 1000L)
```

## Arguments

- x, y:

  numeric vectors. Angles in degrees

- alpha:

  Significance level of the test. Valid levels are `0.01`, `0.05`, and
  `0.1`. This argument may be omitted (`NULL`, the default), in which
  case, a range for the p-value will be returned.

- axial:

  logical. Whether the data are axial, i.e. \\\pi\\-periodical (`TRUE`,
  the default) or directional, i.e. \\2 \pi\\-periodical (`FALSE`).

- quiet:

  logical. Prints the test's decision.

- n_perm:

  integer. Number of permutations

## Value

list

## Details

Watson's two-sample test of homogeneity is performed, and the results
are printed. If alpha is specified and non-zero, the test statistic is
printed along with the critical value and decision. If alpha is omitted,
the test statistic is printed and a range for the p-value of the test is
given.

Critical values for the test statistic are obtained using the asymptotic
distribution of the test statistic. It is recommended to use the
obtained critical values and ranges for p-values only for combined
sample sizes in excess of 17. Tables are available for smaller sample
sizes and can be found in Mardia (1972) for instance.

## See also

Other Tests:
[`ar_test()`](https://tobiste.github.io/tectonicr/reference/ar_test.md),
[`kuiper_test()`](https://tobiste.github.io/tectonicr/reference/kuiper_test.md),
[`norm_chisq()`](https://tobiste.github.io/tectonicr/reference/norm_chisq.md),
[`rayleigh-test`](https://tobiste.github.io/tectonicr/reference/rayleigh-test.md),
[`watson_test()`](https://tobiste.github.io/tectonicr/reference/watson_test.md),
[`watson_wheeler_test_perm()`](https://tobiste.github.io/tectonicr/reference/watson_wheeler_test_perm.md),
[`weighted-rayleigh-test`](https://tobiste.github.io/tectonicr/reference/weighted-rayleigh-test.md)

## Examples

``` r
set.seed(20250411)
x1 <- c(35, 45, 50, 55, 60, 70, 85, 95, 105, 120)
x2 <- c(75, 80, 90, 100, 110, 130, 135, 140, 150, 160, 165)
watson_two_test(x1, x2, axial = FALSE)
#> $statistic
#> [1] 0.1457431
#> 
#> $p.value
#> [1] "P-value > 0.10"
#> 
watson_two_test_perm(x1, x2, axial = FALSE)
#> $statistic
#> [1] 0.1457431
#> 
#> $p.value
#> [1] 0.1208791
#> 

data1 <- rvm(n=20, mean = 0, kappa=3)
data2 <- rvm(n=20, mean = 90, kappa=2)
watson_two_test(data1, data2, axial = FALSE)
#> $statistic
#> [1] 0.61775
#> 
#> $p.value
#> [1] "P-value < 0.001"
#> 
watson_two_test_perm(data1, data2, axial = FALSE)
#> $statistic
#> [1] 0.61775
#> 
#> $p.value
#> [1] 0.000999001
#> 


# San Andreas Fault Data:
data(san_andreas)
data("nuvel1")
PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
sa.por <- PoR_shmax(san_andreas, PoR, "right")
watson_two_test(sa.por$azi, 135, alpha = 0.05)
#> Do not reject null hypothesis
#> $statistic
#> [1] 0.08340728
#> 
#> $p.value
#> [1] 0.187
#> 
watson_two_test_perm(sa.por$azi, rvm(100, 135, 10))
#> $statistic
#> [1] 0.201153
#> 
#> $p.value
#> [1] 0.02897103
#> 
```
