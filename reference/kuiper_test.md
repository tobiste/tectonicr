# Kuiper Test of Circular Uniformity

A statistical test used to determine whether a set of angular or
circular data points (such as times of day, compass directions, or
degrees) are spread out evenly around a circle or if they cluster in
some way.

## Usage

``` r
kuiper_test(x, alpha = 0.05, axial = TRUE, quiet = FALSE)
```

## Arguments

- x:

  numeric vector. Values in degrees

- alpha:

  Significance level of the test. Valid levels are `0.01`, `0.05`, and
  `0.1`. This argument may be omitted (`NULL`, the default), in which
  case, a range for the p-value will be returned.

- axial:

  logical. Whether the data are axial, i.e. \\\pi\\-periodical (`TRUE`,
  the default) or directional, i.e. \\2 \pi\\-periodical (`FALSE`). In
  case of axial data, the angles will be doubled for the test.

- quiet:

  logical. Prints the test's decision.

## Value

A list with class `"htest"` containing the following components:

- `statistic`:

  the value of the D statistic

- `parameter`:

  the number of data

- `p.value`:

  the p-value for the test

- `estimate`:

  the D+ and D- values

- `alternative`:

  a character string describing the alternative hypothesis.

- `method`:

  a character string indicating what type of Kuiper test was performed.

- `data.name`:

  a character string giving the name(s) of the data.

- `alpha`:

  the significance level

- `reject`:

  logical. whether to reject the null-hypothesis based on the
  significance level.

## Details

The **Null Hypothesis** (\\H_0\\): The data are distributed completely
uniformly (randomly and evenly) around the circle.

The **Alternative Hypothesis** (\\H_1\\): The data are not uniform and
show a preference, clustering, or pattern somewhere on the circle.

The Test Statistic (V or \\D^{+} + D^{-}\\): It measures the greatest
positive and negative differences between your data's empirical
cumulative distribution and a theoretical uniform distribution.

### Interpreting the Results

- High Test Statistic / Low p-value (\\p \< \alpha\\, typically 0.05):
  You reject the null hypothesis. This means your data are not uniform;
  they have a significant preferred direction, grouping, or non-random
  pattern on the circle.

- Low Test Statistic / High p-value (\\p \ge 0.05\\): You fail to reject
  the null hypothesis. There is no strong evidence to say the data are
  different from a flat, uniform distribution. The points appear random
  across the circle.

## Note

Kuiper's test statistic is a rotation-invariant Kolmogorov-type test
statistic. The critical values of a modified Kuiper's test statistic are
used according to the tabulation given in Stephens (1970).

## See also

Other Tests:
[`ar_test()`](https://tobiste.github.io/tectonicr/reference/ar_test.md),
[`equal_kappa_test()`](https://tobiste.github.io/tectonicr/reference/equal_kappa_test.md),
[`hodges_ajne_test()`](https://tobiste.github.io/tectonicr/reference/hodges_ajne_test.md),
[`norm_chisq()`](https://tobiste.github.io/tectonicr/reference/norm_chisq.md),
[`rao_spacing_test()`](https://tobiste.github.io/tectonicr/reference/rao_spacing_test.md),
[`rayleigh-test`](https://tobiste.github.io/tectonicr/reference/rayleigh-test.md),
[`uniform_scores_test()`](https://tobiste.github.io/tectonicr/reference/uniform_scores_test.md),
[`watson_test()`](https://tobiste.github.io/tectonicr/reference/watson_test.md),
[`watson_two_sample`](https://tobiste.github.io/tectonicr/reference/watson_two_sample.md),
[`watson_wheeler_test_perm()`](https://tobiste.github.io/tectonicr/reference/watson_wheeler_test_perm.md),
[`weighted-rayleigh-test`](https://tobiste.github.io/tectonicr/reference/weighted-rayleigh-test.md)

## Examples

``` r
# Example data from Mardia and Jupp (1999), pp. 93
kuiper_test(homing, alpha = .05)
#> 
#>  Kuiper's test of uniformity
#> 
#> data:  homing
#> V = 2.2621, n = 10, p-value = 0.001398
#> alternative hypothesis: non-uniform distribution
#> sample estimates:
#>        D+        D- 
#> 0.3611111 0.3055556 
#> 

# San Andreas Fault Data:
data(san_andreas)
data("nuvel1")
PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
sa.por <- PoR_shmax(san_andreas, PoR, "right")
kuiper_test(sa.por$azi.PoR, alpha = .05)
#> 
#>  Kuiper's test of uniformity
#> 
#> data:  sa.por$azi.PoR
#> V = 16.605, n = 1126, p-value < 2.2e-16
#> alternative hypothesis: non-uniform distribution
#> sample estimates:
#>          D+          D- 
#> 0.007774912 0.484679833 
#> 
```
