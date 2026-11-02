# Watson's \\U^2\\ Test for Goodness-of-Fit against known Distribution

A non-parametric statistical test used for circular data to determine
whether a sample fits a specified theoretical distribution

## Usage

``` r
watson_test(
  x,
  alpha = 0.05,
  dist = c("uniform", "vonmises"),
  axial = TRUE,
  quiet = FALSE
)
```

## Arguments

- x:

  numeric vector. Values in degrees

- alpha:

  Significance level of the test ranging between 0 and 1. For
  `dist = "vonmises"` it must be one of `0.1`, `0.05`, `0.01` or `0`
  (only the tabulated p-value range is returned).

- dist:

  Distribution to test for. The default, `"uniform"`, is the circular
  uniform distribution. `"vonmises"` tests the von Mises distribution.

- axial:

  logical. Whether the data are axial, i.e. \\\pi\\-periodical (`TRUE`,
  the default) or directional, i.e. \\2 \pi\\-periodical (`FALSE`). In
  case of axial data, the angles will be doubled for the test.

- quiet:

  logical. Prints the test's decision.

## Value

A list with class `"htest"` containing the following components:

- `statistic`:

  the value of the U2 statistic (modified for `dist = "uniform"`)

- `parameter`:

  the number of data

- `p.value`:

  the p-value for the test (`dist = "uniform"` only)

- `p.range`:

  the range of the p-value from the tabulated critical values
  (`dist = "vonmises"` only)

- `estimate`:

  the MLE of the kappa value for the von Mises distribution

- `alternative`:

  a character string describing the alternative hypothesis.

- `method`:

  a character string indicating what type of Watson test was performed.

- `data.name`:

  a character string giving the name(s) of the data.

- `critical`:

  the critical value at the significance level

- `alpha`:

  the significance level

- `reject`:

  logical. whether to reject the null-hypothesis based on the
  significance level.

## Details

### Hypotheses

**Null Hypothesis** (\\H_0\\): The circular sample comes from a
specified theoretical distribution (such as a uniform distribution or a
specific von Mises distribution).

**Alternative Hypothesis** (\\H_1\\): The circular sample does not
follow the specified theoretical distribution.

### Interpretation

To interpret the output of Watson's \\U^2\\ test, compare your
calculated \\U^2\\ test statistic to the critical value from Watson's
goodness-of-fit/homogeneity tables at your chosen significance level
(\\\alpha\\, commonly set to 0.05), or check the resulting p-value:

- If \\U^2\_{\text{calculated}} \> U^2\_{\text{critical}}\\ (or p \<
  \\\alpha\\): Reject the null hypothesis (\\H_0\\). Conclude that the
  data significantly deviates from the theoretical distribution.

- If \\U^2\_{\text{calculated}} \le U^2\_{\text{critical}}\\ (or \\p \ge
  \alpha\\): Fail to reject the null hypothesis (\\H_0\\). There is not
  enough evidence to claim the data deviates from the expected model.

## Note

Watson's test statistic is a rotation-invariant Cramer - von Mises test.

## References

Mardia and Jupp (2000). Directional Statistics. John Wiley and Sons.

## See also

[vonmises](https://tobiste.github.io/tectonicr/reference/vonmises.md)
and [cunif](https://tobiste.github.io/tectonicr/reference/cunif.md)

Other One-Sample-Tests:
[`hodges_ajne_test()`](https://tobiste.github.io/tectonicr/reference/hodges_ajne_test.md),
[`kuiper_test()`](https://tobiste.github.io/tectonicr/reference/kuiper_test.md),
[`norm_chisq()`](https://tobiste.github.io/tectonicr/reference/norm_chisq.md),
[`rao_spacing_test()`](https://tobiste.github.io/tectonicr/reference/rao_spacing_test.md),
[`rayleigh-test`](https://tobiste.github.io/tectonicr/reference/rayleigh-test.md),
[`weighted-rayleigh-test`](https://tobiste.github.io/tectonicr/reference/weighted-rayleigh-test.md)

## Examples

``` r
# Example data from Mardia and Jupp (1999), pp. 93
watson_test(homing, axial = FALSE, alpha = .05)
#> Do Not Reject Null Hypothesis
#> 
#>  Watson's U2 test against uniformity
#> 
#> data:  homing
#> U2 = 0.11536, n = 10, p-value = 0.2049
#> alternative hypothesis: non-uniform distribution
#> 

# San Andreas Fault Data:
data(san_andreas)
data("nuvel1")
PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
sa.por <- PoR_shmax(san_andreas, PoR, "right")
watson_test(sa.por$azi.PoR, alpha = .05)
#> Reject Null Hypothesis
#> 
#>  Watson's U2 test against uniformity
#> 
#> data:  sa.por$azi.PoR
#> U2 = 31.687, n = 1126, p-value < 2.2e-16
#> alternative hypothesis: non-uniform distribution
#> 
watson_test(sa.por$azi.PoR, alpha = .05, dist = "vonmises")
#> Reject Null Hypothesis
#> 
#>  Watson's U2 test against von Mises distribution (P-value < 0.01)
#> 
#> data:  sa.por$azi.PoR
#> U2 = 0.50844, n = 1126
#> alternative hypothesis: Sample does not follow a von Mises distribution
#> sample estimates:
#> Kappa MLE 
#>  2.184977 
#> 
```
