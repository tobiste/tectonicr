# Equal \\\kappa\\ Test

Test whether samples of directions share the concentration

## Usage

``` r
equal_kappa_test(
  samples,
  method = c("f", "bartlett"),
  alpha = 0.05,
  axial = TRUE
)
```

## Arguments

- samples:

  list of samples of orientations

- method:

  character. One of

  `"f"`

  :   (default for two samples): \\F = ((n_1 - R_1) / nu_1) / ((n_2 -
      R_2) / nu_2)\\, two-sided, \\F(\nu_1, \nu_2)\\ (MJ eq. 10.5.3;
      Watson, 1956).

  `"bartlett"`

  :   (default for more samples): Bartlett's test of the dispersions
      \\(n_i - R_i) / \nu_i\\, chi-square with \\q - 1\\ degrees of
      freedom (MJ eq. 10.6.21; McFadden & Lowes, 1981, eq. 64-66).

- alpha:

  Significance level of the test. Valid levels are `0.01`, `0.05`, and
  `0.1`. This argument may be omitted (`NULL`, the default), in which
  case, a range for the p-value will be returned.

- axial:

  logical. Whether the data are axial, i.e. \\\pi\\-periodical (`TRUE`,
  the default) or directional, i.e. \\2 \pi\\-periodical (`FALSE`). In
  case of axial data, the angles will be doubled for the test.

## Value

A list with class `"htest"`

## Details

For concentrated data \\2 \kappa (n_i - R_i)\\ is approximately
chi-square with \\\nu_i = (p - 1)(n_i - 1)\\ degrees of freedom in 2
dimensions.

## Note

Both need concentrated data (\\\bar{R} \ge 0.7\\; MJ eq. 10.6.24-10.6.27
and 7.3.24-7.3.25 give variance-stabilised tests for dispersed data.

## See also

Other Tests:
[`ar_test()`](https://tobiste.github.io/tectonicr/reference/ar_test.md),
[`hodges_ajne_test()`](https://tobiste.github.io/tectonicr/reference/hodges_ajne_test.md),
[`kuiper_test()`](https://tobiste.github.io/tectonicr/reference/kuiper_test.md),
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
set.seed(20250411)
x1 <- c(35, 45, 50, 55, 60, 70, 85, 95, 105, 120)
x2 <- c(75, 80, 90, 100, 110, 130, 135, 140, 150, 160, 165)
equal_kappa_test(list(x1, x2))
#> 
#>  F test of equal concentrations
#> 
#> data:  list(x1, x2)
#> F = 0.9935, num df = 9, denom df = 10, p-value = 0.9986
#> alternative hypothesis: true ratio of concentrations is not equal to 1
#> sample estimates:
#> ratio of concentrations 
#>               0.9934996 
#> 

#' # San Andreas Fault Data:
data(san_andreas)
data("nuvel1")
PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
sa.por <- PoR_shmax(san_andreas, PoR, "right")
equal_kappa_test(list(sa.por$azi.PoR, rvm(100, 135, 10)))
#> 
#>  F test of equal concentrations
#> 
#> data:  list(sa.por$azi.PoR, rvm(100, 135, 10))
#> F = 0.99836, num df = 1125, denom df = 99, p-value = 0.9581
#> alternative hypothesis: true ratio of concentrations is not equal to 1
#> sample estimates:
#> ratio of concentrations 
#>               0.9983557 
#> 
```
