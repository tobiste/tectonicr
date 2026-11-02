# Uniform scores test

q-sample uniform-scores test tests the equality of given distributions,
ie. whether samples of directions share the mean direction. The test is
the distribution-free test of Mardia, Watson & Wheeler on the uniform
scores \\2 \pi r_i / n\\ of the ranks in the pooled sample.

## Usage

``` r
uniform_scores_test(samples, alpha = 0.05, axial = TRUE)
```

## Arguments

- samples:

  list of samples of orientations

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

\$\$W = 2 \sum\_{i=1}^q R_i^2 / n_i\$\$ with \\2(q - 1)\\ degrees of
freedom, \\R^\* = 2 (n - 1) R_1^2 / (n_1 n_2)\\ for two samples

## See also

Other Tests:
[`ar_test()`](https://tobiste.github.io/tectonicr/reference/ar_test.md),
[`equal_kappa_test()`](https://tobiste.github.io/tectonicr/reference/equal_kappa_test.md),
[`hodges_ajne_test()`](https://tobiste.github.io/tectonicr/reference/hodges_ajne_test.md),
[`kuiper_test()`](https://tobiste.github.io/tectonicr/reference/kuiper_test.md),
[`norm_chisq()`](https://tobiste.github.io/tectonicr/reference/norm_chisq.md),
[`rao_spacing_test()`](https://tobiste.github.io/tectonicr/reference/rao_spacing_test.md),
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
uniform_scores_test(list(x1, x2))
#> 
#>  Uniform scores test
#> 
#> data:  list(x1, x2)
#> W = 3.6783, n = 21, p-value = 0.159
#> alternative hypothesis: Samples come from different distributions
#> 

#' # San Andreas Fault Data:
data(san_andreas)
data("nuvel1")
PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
sa.por <- PoR_shmax(san_andreas, PoR, "right")
uniform_scores_test(list(sa.por$azi.PoR, rvm(100, 135, 10)))
#> 
#>  Uniform scores test
#> 
#> data:  list(sa.por$azi.PoR, rvm(100, 135, 10))
#> W = 0.9193, n = 1226, p-value = 0.6315
#> alternative hypothesis: Samples come from different distributions
#> 
```
