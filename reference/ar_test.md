# Angular Randomisation Test of Homogeneity

Performs Angular Randomisation Test for homogeneity on two samples of
circular data after Ruxton et al. (2023). P-values are estimated using
permutation.

## Usage

``` r
ar_test(x, y, n_perm = 1000L, axial = FALSE, alpha = NULL)
```

## Arguments

- x, y:

  numeric vectors. Angles in degrees

- n_perm:

  integer. Number of permutations

- axial:

  logical. Whether the data are axial, i.e. \\\pi\\-periodical (`TRUE`)
  or directional, i.e. \\2 \pi\\-periodical (`FALSE`, the default).

- alpha:

  (optional) numeric. Significance level of the test (values between 0
  and 1).

## Value

list containing the test statistic, the p-value, the significance value
`alpha` and a logical decision whether to reject the null hypothesis or
not.

## Details

**Null Hypothesis** (\\H_0\\): The two circular samples share an
identical underlying probability distribution.

**Alternative Hypothesis** (\\H\_{1}\\): The two samples come from
different distributions.

### Interpretation

- Small p-value (\\p \< \alpha\\, e.g., \<0.05): Reject the null
  hypothesis. This indicates strong evidence that the two samples come
  from different circular distributions (differing in central
  tendency/mean direction or shape).

- Large p-value (\\p \ge \alpha\\): Fail to reject the null hypothesis;
  there is insufficient evidence to claim the two circular samples
  differ.

## Note

**Concentration Differences**: The test can suffer from markedly lower
statistical power if the underlying unimodal distributions differ by
concentration (dispersion/spread) rather than location—especially with
small, uneven sample sizes where the smaller sample comes from the more
concentrated distribution.

**Axial/Multimodal Data**: ART performs poorly and loses power when
applied to axially symmetric or symmetrically multimodal distributions.

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
ar_test(x1, x2)
#> $statistic
#> [1] 104.3707
#> 
#> $p.value
#> [1] 0.003996004
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
ar_test(sa.por$azi.PoR, rvm(100, 135, 10), axial = TRUE, alpha = 0.05)
#> $statistic
#> [1] 92652.2
#> 
#> $p.value
#> [1] 0.993007
#> 
#> $alpha
#> [1] 0.05
#> 
#> $reject
#> [1] FALSE
#> 
```
