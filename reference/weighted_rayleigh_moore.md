# Moore (1980) rank-weighted Rayleigh test

Rank-based modification of the Rayleigh test for weighted vector data
(Moore 1980). Vectors are ranked by weight and the rank replaces the
vector length. Only the order of `w` matters, not its magnitude.

## Usage

``` r
weighted_rayleigh_moore(
  x,
  mu = NULL,
  w = NULL,
  axial = TRUE,
  alpha = 0.05,
  quiet = FALSE,
  ties = "random"
)

weighted_rayleigh_moore_perm(
  x,
  mu = NULL,
  w = NULL,
  axial = TRUE,
  ties = "random",
  n_perm = 10000L
)
```

## Arguments

- x:

  numeric vector. Values in degrees

- mu:

  (optional) The specified or known mean direction (in degrees) in
  alternative hypothesis.

- w:

  numeric vector of weights, same length as `x`. If `NULL`, the standard
  Rayleigh test is returned.

- axial:

  logical. Whether the data are axial, i.e. \\\pi\\-periodical (`TRUE`,
  the default) or directional, i.e. \\2 \pi\\-periodical (`FALSE`). In
  case of axial data, the angles will be doubled for the test.

- alpha:

  Significance level of the test ranging between 0 and 1.

- quiet:

  logical. Prints the test's decision.

- ties:

  tie-breaking passed to [`rank()`](https://rdrr.io/r/base/rank.html).
  `"random"` follows Moore (1980).

- n_perm:

  integer. Number of simulated uniform samples.

## Note

Random tie-breaking makes results depend on the seed when weights tie.
Call [`set.seed()`](https://rdrr.io/r/base/Random.html) or use
`ties = "first"` if you need reproducibility.

## References

Moore, B. R. (1980) A modification of the Rayleigh test for vector data.
Biometrika 67, 175-180.

## See also

Other One-Sample-Tests:
[`hodges_ajne_test()`](https://tobiste.github.io/tectonicr/reference/hodges_ajne_test.md),
[`kuiper_test()`](https://tobiste.github.io/tectonicr/reference/kuiper_test.md),
[`norm_chisq()`](https://tobiste.github.io/tectonicr/reference/norm_chisq.md),
[`rao_spacing_test()`](https://tobiste.github.io/tectonicr/reference/rao_spacing_test.md),
[`rayleigh-test`](https://tobiste.github.io/tectonicr/reference/rayleigh-test.md),
[`watson_test()`](https://tobiste.github.io/tectonicr/reference/watson_test.md),
[`weighted-rayleigh-test`](https://tobiste.github.io/tectonicr/reference/weighted-rayleigh-test.md)

## Examples

``` r
set.seed(20250411)

# Load data
data("cpm_models")
data(san_andreas)
PoR <- equivalent_rotation(cpm_models[["NNR-MORVEL56"]], "na", "pa")
sa.por <- PoR_shmax(san_andreas, PoR, "right")
data("iceland")
PoR.ice <- equivalent_rotation(cpm_models[["NNR-MORVEL56"]], "eu", "na")
ice.por <- PoR_shmax(iceland, PoR.ice, "out")
data("tibet")
PoR.tib <- equivalent_rotation(cpm_models[["NNR-MORVEL56"]], "eu", "in")
tibet.por <- PoR_shmax(tibet, PoR.tib, "in")

# GOF test:
weighted_rayleigh_moore(tibet.por$azi.PoR, mu = 90, w = 1 / tibet$unc)
#> Reject Null Hypothesis
#> 
#>  Moore (1980) rank-weighted V-test of uniformity
#> 
#> data:  tibet.por$azi.PoR
#> S = 22.585, n = 1165, p-value < 2.2e-16
#> alternative hypothesis: unimodal distribution with mean direction 90
#> sample estimates:
#>      C 
#> 366868 
#> 
weighted_rayleigh_moore(ice.por$azi.PoR, mu = 0, w = 1 / iceland$unc)
#> Reject Null Hypothesis
#> 
#>  Moore (1980) rank-weighted V-test of uniformity
#> 
#> data:  ice.por$azi.PoR
#> S = 9.7672, n = 490, p-value < 2.2e-16
#> alternative hypothesis: unimodal distribution with mean direction 0
#> sample estimates:
#>        C 
#> 43316.57 
#> 
weighted_rayleigh_moore(sa.por$azi.PoR, mu = 135, w = 1 / san_andreas$unc)
#> Reject Null Hypothesis
#> 
#>  Moore (1980) rank-weighted V-test of uniformity
#> 
#> data:  sa.por$azi.PoR
#> S = 30.825, n = 1126, p-value < 2.2e-16
#> alternative hypothesis: unimodal distribution with mean direction 135
#> sample estimates:
#>        C 
#> 475800.6 
#> 

weighted_rayleigh_moore_perm(tibet.por$azi.PoR, mu = 90, w = 1 / tibet$unc)
#> 
#>  Moore (1980) rank-weighted V-test (Monte Carlo p-value based on 10000
#>  uniform samples)
#> 
#> data:  tibet.por$azi.PoR
#> C = 360059, n = 1165, n_perm = 10000, p-value = 9.999e-05
#> alternative hypothesis: unimodal distribution with mean direction 90
#> 
weighted_rayleigh_moore_perm(ice.por$azi.PoR, mu = 0, w = 1 / iceland$unc)
#> 
#>  Moore (1980) rank-weighted V-test (Monte Carlo p-value based on 10000
#>  uniform samples)
#> 
#> data:  ice.por$azi.PoR
#> C = 42899, n = 490, n_perm = 10000, p-value = 9.999e-05
#> alternative hypothesis: unimodal distribution with mean direction 0
#> 
weighted_rayleigh_moore_perm(sa.por$azi.PoR, mu = 135, w = 1 / san_andreas$unc)
#> 
#>  Moore (1980) rank-weighted V-test (Monte Carlo p-value based on 10000
#>  uniform samples)
#> 
#> data:  sa.por$azi.PoR
#> C = 475967, n = 1126, n_perm = 10000, p-value = 9.999e-05
#> alternative hypothesis: unimodal distribution with mean direction 135
#> 
```
