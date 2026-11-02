# Watson-Williams Test of Homogeneity of Means

Performs the Watson-Williams test for homogeneity of means between
several samples of circular data.

## Usage

``` r
watson_williams_test(samples, alpha = 0.05, axial = TRUE)
```

## Arguments

- samples:

  list of samples of orientations

- alpha:

  Significance level of the test ranging between 0 and 1. For
  `dist = "vonmises"` it must be one of `0.1`, `0.05`, `0.01` or `0`
  (only the tabulated p-value range is returned).

- axial:

  logical. Whether the data are axial, i.e. \\\pi\\-periodical (`TRUE`,
  the default) or directional, i.e. \\2 \pi\\-periodical (`FALSE`). In
  case of axial data, the angles will be doubled for the test.

## Value

A list with class `"htest"`

## Details

The Watson-Williams test for the homogeneity of means between two or
more groups is performed and the results are printed. The null
hypothesis is that means are equal across groups.

The assumptions are that: (1) the samples are drawn from populations
with a von Mises distribution; (2) the parameter of concentration has
the same value in all populations; (3) this parameter is sufficiently
large (i.e. \> 2 for two samples, \> 1 for more samples). Assumptions 2
and 3 are checked and a warning is issued if they are not met.

## See also

Other Two-Sample-Tests:
[`ar_test()`](https://tobiste.github.io/tectonicr/reference/ar_test.md),
[`equal_kappa_test()`](https://tobiste.github.io/tectonicr/reference/equal_kappa_test.md),
[`uniform_scores_test()`](https://tobiste.github.io/tectonicr/reference/uniform_scores_test.md),
[`walraff_test()`](https://tobiste.github.io/tectonicr/reference/walraff_test.md),
[`watson_two_sample`](https://tobiste.github.io/tectonicr/reference/watson_two_sample.md),
[`watson_wheeler_test_perm()`](https://tobiste.github.io/tectonicr/reference/watson_wheeler_test_perm.md)

## Examples

``` r
x1 <- c(35, 45, 50, 55, 60, 70, 85, 95, 105, 120)
x2 <- c(75, 80, 90, 100, 110, 130, 135, 140, 150, 160, 165)
watson_williams_test(list(x1, x2), axial = FALSE)
#> 
#>  Watson-Williams test for homogeneity of mean directions
#> 
#> data:  1 and 2
#> F = 13.935, df1 = 1, df2 = 19, p-value = 0.00141
#> sample estimates:
#> mean direction 1 mean direction 2 
#>         71.60242        121.56169 
#> 

data1 <- rvm(n=20, mean = 0, kappa=3)
data2 <- rvm(n=20, mean = 90, kappa=2)
watson_williams_test(list(data1, data2), axial = FALSE)
#> 
#>  Watson-Williams test for homogeneity of mean directions
#> 
#> data:  1 and 2
#> F = 46.167, df1 = 1, df2 = 38, p-value = 4.669e-08
#> sample estimates:
#> mean direction 1 mean direction 2 
#>         6.351635        87.388648 
#> 

# San Andreas Fault Data:
data(san_andreas)
data("nuvel1")
PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
sa.por <- PoR_shmax(san_andreas, PoR, "right")
watson_williams_test(list(sa.por$azi.PoR, rvm(100, 135, 10)))
#> Warning: Concentrations (1.76, 2.77) differ between samples (p = 0.00458); the test might not be applicable
#> 
#>  Watson-Williams test for homogeneity of mean axes
#> 
#> data:  1 and 2
#> F = 13.208, df1 = 1, df2 = 1224, p-value = 0.0002903
#> sample estimates:
#> mean direction 1 mean direction 2 
#>         140.8778         132.5831 
#> 
```
