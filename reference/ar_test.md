# Angular Randomisation Test of Homogeneity

Performs Angular Randomisation Test for homogeneity on two samples of
circular data after Ruxton et al. (2023). P-values are estimated using
permutation.

## Usage

``` r
ar_test(x, y, n_perm = 1000L, axial = TRUE)
```

## Arguments

- x, y:

  numeric vectors. Angles in degrees

- n_perm:

  integer. Number of permutations

- axial:

  logical. Whether the data are axial, i.e. \\\pi\\-periodical (`TRUE`,
  the default) or directional, i.e. \\2 \pi\\-periodical (`FALSE`).

## Value

list containing the test statistic and the p-value

## References

Ruxton, G.D., Malkemper, E.P. & Landler, L. Evaluating the power of a
recent method for comparing two circular distributions: an alternative
to the Watson U2 test. Sci Rep 13, 10007 (2023).
https://doi.org/10.1038/s41598-023-36960-1

## See also

Other Tests:
[`kuiper_test()`](https://tobiste.github.io/tectonicr/reference/kuiper_test.md),
[`norm_chisq()`](https://tobiste.github.io/tectonicr/reference/norm_chisq.md),
[`rayleigh-test`](https://tobiste.github.io/tectonicr/reference/rayleigh-test.md),
[`watson_test()`](https://tobiste.github.io/tectonicr/reference/watson_test.md),
[`watson_two_sample`](https://tobiste.github.io/tectonicr/reference/watson_two_sample.md),
[`watson_wheeler_test_perm()`](https://tobiste.github.io/tectonicr/reference/watson_wheeler_test_perm.md),
[`weighted-rayleigh-test`](https://tobiste.github.io/tectonicr/reference/weighted-rayleigh-test.md)

## Examples

``` r
set.seed(20250411)
x1 <- c(35, 45, 50, 55, 60, 70, 85, 95, 105, 120)
x2 <- c(75, 80, 90, 100, 110, 130, 135, 140, 150, 160, 165)
ar_test(x1, x2, axial = FALSE)
#> $statistic
#> [1] 104.3707
#> 
#> $p.value
#> [1] 0.003996004
#> 

# San Andreas Fault Data:
data(san_andreas)
data("nuvel1")
PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
sa.por <- PoR_shmax(san_andreas, PoR, "right")
ar_test(sa.por$azi.PoR, rvm(100, 135, 10))
#> $statistic
#> [1] 92652.2
#> 
#> $p.value
#> [1] 0.993007
#> 
```
