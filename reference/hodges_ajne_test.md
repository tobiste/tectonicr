# Hodges-Ajne Test

Hodges-Ajne test of uniformity on the circle

## Usage

``` r
hodges_ajne_test(x, alpha = 0.05, axial = TRUE)
```

## Arguments

- x:

  numeric vector. Values in degrees

- alpha:

  Significance level of the test ranging between 0 and 1.

- axial:

  logical. Whether the data are axial, i.e. \\\pi\\-periodical (`TRUE`,
  the default) or directional, i.e. \\2 \pi\\-periodical (`FALSE`). In
  case of axial data, the angles will be doubled for the test.

## Value

A list with class `"htest"` containing the following components:

- `statistic`:

  the value of the m statistic

- `parameter`:

  the number of data

- `p.value`:

  the p-value for the test

- `alternative`:

  a character string describing the alternative hypothesis.

- `method`:

  a character string indicating what type of test was performed.

- `data.name`:

  a character string giving the name(s) of the data.

## Details

The statistic \\m\\ is the smallest number of angles in a half-circle.
Its exact distribution is \$\$P(m \le t) = 2^{(1 - n)} (n - 2t) C(n,
t)\$\$ for \\t \< n / 3\\ (Mardia and Jupp, 2000, sec. 6.3.4), which
gives the p-value (1 above \\n / 3\\). Small \\m\\ indicates a preferred
direction.

## See also

Other One-Sample-Tests:
[`kuiper_test()`](https://tobiste.github.io/tectonicr/reference/kuiper_test.md),
[`norm_chisq()`](https://tobiste.github.io/tectonicr/reference/norm_chisq.md),
[`rao_spacing_test()`](https://tobiste.github.io/tectonicr/reference/rao_spacing_test.md),
[`rayleigh-test`](https://tobiste.github.io/tectonicr/reference/rayleigh-test.md),
[`watson_test()`](https://tobiste.github.io/tectonicr/reference/watson_test.md),
[`weighted-rayleigh-test`](https://tobiste.github.io/tectonicr/reference/weighted-rayleigh-test.md)

## Examples

``` r
hodges_ajne_test(rvm(100, mean = 0, kappa = 10))
#> 
#>  Hodges-Ajne test of uniformity
#> 
#> data:  rvm(100, mean = 0, kappa = 10)
#> m = 1, n = 100, p-value < 2.2e-16
#> alternative hypothesis: unimodal distribution
#> 

# San Andreas Fault Data:
data(san_andreas)
data("nuvel1")
PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
sa.por <- PoR_shmax(san_andreas, PoR, "right")
hodges_ajne_test(sa.por$azi.PoR, alpha = .05)
#> 
#>  Hodges-Ajne test of uniformity
#> 
#> data:  sa.por$azi.PoR
#> m = 76, n = 1126, p-value < 2.2e-16
#> alternative hypothesis: unimodal distribution
#> 
```
