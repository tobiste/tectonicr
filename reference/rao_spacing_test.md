# Rao Spacing Test

Rao's spacing test of uniformity on the circle

## Usage

``` r
rao_spacing_test(x, alpha = 0.05, n_perm = 9999L, axial = TRUE)
```

## Arguments

- x:

  numeric vector. Values in degrees

- alpha:

  Significance level of the test ranging between 0 and 1.

- n_perm:

  integer. Number of simulated uniform samples.

- axial:

  logical. Whether the data are axial, i.e. \\\pi\\-periodical (`TRUE`,
  the default) or directional, i.e. \\2 \pi\\-periodical (`FALSE`). In
  case of axial data, the angles will be doubled for the test.

## Value

A list with class `"htest"` containing the following components:

- `statistic`:

  the value of the L statistic

- `parameter`:

  the number of data and the number of simulations

- `p.value`:

  the p-value for the test

- `alternative`:

  a character string describing the alternative hypothesis.

- `method`:

  a character string indicating what type of test was performed.

- `data.name`:

  a character string giving the name(s) of the data.

## Details

Test statistic is \$\$L = 1/2 \sum(\|T_i - 360 / n\|)\$\$ in degrees,
with the spacings \\T_i\\ between consecutive angles (Mardia & Jupp
2000, eq. 6.3.45). The p-value is the proportion of `n_perm` simulated
uniform samples with a larger \\L\\

## See also

Other One-Sample-Tests:
[`hodges_ajne_test()`](https://tobiste.github.io/tectonicr/reference/hodges_ajne_test.md),
[`kuiper_test()`](https://tobiste.github.io/tectonicr/reference/kuiper_test.md),
[`norm_chisq()`](https://tobiste.github.io/tectonicr/reference/norm_chisq.md),
[`rayleigh-test`](https://tobiste.github.io/tectonicr/reference/rayleigh-test.md),
[`watson_test()`](https://tobiste.github.io/tectonicr/reference/watson_test.md),
[`weighted-rayleigh-test`](https://tobiste.github.io/tectonicr/reference/weighted-rayleigh-test.md)

## Examples

``` r
set.seed(20250411)
rao_spacing_test(rvm(100, mean = 0, kappa = 10))
#> 
#>  Rao's spacing test of uniformity with simulated p-value (based on 9999
#>  replicates)
#> 
#> data:  rvm(100, mean = 0, kappa = 10)
#> L = 230.53, n = 100, n_perm = 9999, p-value = 1e-04
#> alternative hypothesis: non-uniform distribution
#> 

# San Andreas Fault Data:
data(san_andreas)
data("nuvel1")
PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
sa.por <- PoR_shmax(san_andreas, PoR, "right")
rao_spacing_test(sa.por$azi.PoR, alpha = .05)
#> 
#>  Rao's spacing test of uniformity with simulated p-value (based on 9999
#>  replicates)
#> 
#> data:  sa.por$azi.PoR
#> L = 202.73, n = 1126, n_perm = 9999, p-value = 1e-04
#> alternative hypothesis: non-uniform distribution
#> 
```
