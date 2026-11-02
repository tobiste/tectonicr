# Weighted Goodness-of-fit Test for Circular Data

Weighted version of the Rayleigh test (or V0-test) for uniformity
against a distribution with a priori expected von Mises concentration.
`weighted_rayleigh_perm()` estimates p-values by Monte Carlo simulation
of uniform samples with the same weights.

## Usage

``` r
weighted_rayleigh(
  x,
  mu = NULL,
  w = NULL,
  axial = TRUE,
  alpha = 0.05,
  quiet = FALSE
)

weighted_rayleigh_perm(x, mu = NULL, w = NULL, axial = TRUE, n_perm = 1000L)

v0_test(x, mu, w = NULL, axial = TRUE, n_perm = 1000L)
```

## Arguments

- x:

  numeric vector. Values in degrees

- mu:

  (optional) The specified or known mean direction (in degrees) in
  alternative hypothesis.

- w:

  numeric vector weights of length `length(x)`. If `NULL`, the
  non-weighted Rayleigh test is performed.

- axial:

  logical. Whether the data are axial, i.e. \\\pi\\-periodical (`TRUE`,
  the default) or directional, i.e. \\2 \pi\\-periodical (`FALSE`). In
  case of axial data, the angles will be doubled for the test.

- alpha:

  Significance level of the test ranging between 0 and 1.

- quiet:

  logical. Prints the test's decision.

- n_perm:

  integer. Number of simulated uniform samples.

## Value

a list with the components:

- `R` or `C`:

  weighted mean resultant length or, if `mu` is specified, the weighted
  mean cosine of the deviations from `mu`. Large values reject
  uniformity.

- `statistic`:

  Test statistic

- `p.value`:

  significance level of the test statistic

- `n_eff`:

  effective sample size

`weighted_rayleigh_perm()` returns a list of class `"htest"`.

## Details

The Null hypothesis is uniformity (randomness). The alternative is a
unimodal distribution with the specified mean direction `mu` (weighted
V-test), or with an unknown mean direction if `mu` is `NULL` (weighted
Rayleigh test). The null hypothesis is rejected if `p.value < alpha`.

Under the null hypothesis the weighted mean cosine has variance \\\sum
w_i^2 / (2 (\sum w_i)^2)\\, so the asymptotic tests use the effective
sample size \\n\_\text{eff} = (\sum w_i)^2 / \sum w_i^2\\, which equals
\\n\\ for equal weights.

## See also

Other One-Sample-Tests:
[`hodges_ajne_test()`](https://tobiste.github.io/tectonicr/reference/hodges_ajne_test.md),
[`kuiper_test()`](https://tobiste.github.io/tectonicr/reference/kuiper_test.md),
[`norm_chisq()`](https://tobiste.github.io/tectonicr/reference/norm_chisq.md),
[`rao_spacing_test()`](https://tobiste.github.io/tectonicr/reference/rao_spacing_test.md),
[`rayleigh-test`](https://tobiste.github.io/tectonicr/reference/rayleigh-test.md),
[`watson_test()`](https://tobiste.github.io/tectonicr/reference/watson_test.md)

## Examples

``` r
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
weighted_rayleigh(tibet.por$azi.PoR, mu = 90, w = 1 / tibet$unc)
#> Reject Null Hypothesis
#> $C
#> [1] 0.5321474
#> 
#> $statistic
#> [1] 16.75814
#> 
#> $p.value
#> [1] 9.206653e-62
#> 
#> $n_eff
#> [1] 495.8583
#> 
weighted_rayleigh(ice.por$azi.PoR, mu = 0, w = 1 / iceland$unc)
#> Reject Null Hypothesis
#> $C
#> [1] 0.3728874
#> 
#> $statistic
#> [1] 7.950463
#> 
#> $p.value
#> [1] 3.500254e-16
#> 
#> $n_eff
#> [1] 227.3
#> 
weighted_rayleigh(sa.por$azi.PoR, mu = 135, w = 1 / san_andreas$unc)
#> Reject Null Hypothesis
#> $C
#> [1] 0.8042366
#> 
#> $statistic
#> [1] 24.38294
#> 
#> $p.value
#> [1] 1.404208e-128
#> 
#> $n_eff
#> [1] 459.5941
#> 

weighted_rayleigh_perm(tibet.por$azi.PoR, mu = 90, w = 1 / tibet$unc)
#> 
#>  Weighted Rayleigh V-test of uniformity (Monte Carlo p-value based on
#>  1000 uniform samples)
#> 
#> data:  tibet.por$azi.PoR
#> C = 0.53215, n = 1165, n_perm = 1000, p-value = 0.000999
#> alternative hypothesis: unimodal distribution with mean direction 90
#> 
weighted_rayleigh_perm(ice.por$azi.PoR, mu = 0, w = 1 / iceland$unc)
#> 
#>  Weighted Rayleigh V-test of uniformity (Monte Carlo p-value based on
#>  1000 uniform samples)
#> 
#> data:  ice.por$azi.PoR
#> C = 0.37289, n = 490, n_perm = 1000, p-value = 0.000999
#> alternative hypothesis: unimodal distribution with mean direction 0
#> 
weighted_rayleigh_perm(sa.por$azi.PoR, mu = 135, w = 1 / san_andreas$unc)
#> 
#>  Weighted Rayleigh V-test of uniformity (Monte Carlo p-value based on
#>  1000 uniform samples)
#> 
#> data:  sa.por$azi.PoR
#> C = 0.80424, n = 1126, n_perm = 1000, p-value = 0.000999
#> alternative hypothesis: unimodal distribution with mean direction 135
#> 
```
