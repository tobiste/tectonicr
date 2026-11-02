# Equal \\\kappa\\ Test

This function tests for the homogeneity of concentration parameters for
multiple samples of directional data.

## Usage

``` r
equal_kappa_test(samples, method = NULL, alpha = 0.05, axial = TRUE)
```

## Arguments

- samples:

  list of samples of orientations

- method:

  character. One of

  `"f"`

  :   (default for two samples): \\F = ((n_1 - R_1) / \nu_1) / ((n_2 -
      R_2) / \nu_2)\\, two-sided, \\F(\nu_1, \nu_2)\\ (Mardia and Jupp,
      2000 eq. 10.5.3; Watson, 1956).

  `"bartlett"`

  :   (default for more samples): Bartlett's test of the dispersions
      \\(n_i - R_i) / \nu_i\\, chi-square with \\q - 1\\ degrees of
      freedom (Mardia and Jupp, 2000 eq. 10.6.21; McFadden & Lowes,
      1981, eq. 64-66).

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

The samples are assumed to have been drawn from von Mises populations.
The null hypothesis tested is that all populations sampled have the same
concentration parameter, kappa \\\kappa\\.

For concentrated data \\2 \kappa (n_i - R_i)\\ is approximately
chi-square with \\\nu_i = n_i - 1\\ degrees of freedom on the circle,
where \\R_i\\ is the resultant length of sample \\i\\.

## Note

Both need concentrated data (\\\bar{R} \ge 0.7\\); MJ eq.
10.6.24-10.6.27 and 7.3.24-7.3.25 give variance-stabilised tests for
dispersed data.

## References

Jammalamadaka, S. Rao and SenGupta, A. (2001). Topics in Circular
Statistics, Section 5.3, World Scientific Press, Singapore.

Mardia, K. and Jupp, P. (1999). Directional Statistics, Section 7.4,
John Wiley and Sons, England.

## See also

Other Two-Sample-Tests:
[`ar_test()`](https://tobiste.github.io/tectonicr/reference/ar_test.md),
[`uniform_scores_test()`](https://tobiste.github.io/tectonicr/reference/uniform_scores_test.md),
[`walraff_test()`](https://tobiste.github.io/tectonicr/reference/walraff_test.md),
[`watson_two_sample`](https://tobiste.github.io/tectonicr/reference/watson_two_sample.md),
[`watson_wheeler_test_perm()`](https://tobiste.github.io/tectonicr/reference/watson_wheeler_test_perm.md),
[`watson_williams_test()`](https://tobiste.github.io/tectonicr/reference/watson_williams_test.md)

## Examples

``` r
set.seed(20250411)
x1 <- c(35, 45, 50, 55, 60, 70, 85, 95, 105, 120)
x2 <- c(75, 80, 90, 100, 110, 130, 135, 140, 150, 160, 165)
equal_kappa_test(list(x1, x2))
#> 
#>  F test of equal concentrations
#> 
#> data:  1 and 2
#> F = 0.78175, num df = 9, denom df = 10, p-value = 0.7216
#> alternative hypothesis: true ratio of concentrations is not equal to 1
#> sample estimates:
#> ratio of concentrations 
#>                0.781748 
#> 

x1 <- rvm(100, mean = 90, kappa = 20)
x2 <- rvm(100, mean = 90, kappa = 10)
x3 <- rvm(100, mean = 90, kappa = 1)
equal_kappa_test(list(x1, x2, x3))
#> 
#>  Bartlett test of equal concentrations
#> 
#> data:  1 and 2 and 3
#> Bartlett's K-squared = 133.57, df = 2, p-value < 2.2e-16
#> sample estimates:
#>   kappa_1   kappa_2   kappa_3 
#> 5.6313939 2.5644762 0.5629104 
#> 

# San Andreas Fault Data:
data(san_andreas)
data("nuvel1")
PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
sa.por <- PoR_shmax(san_andreas, PoR, "right")
equal_kappa_test(list(sa.por$azi.PoR, rvm(100, 135, 10)))
#> 
#>  F test of equal concentrations
#> 
#> data:  1 and 2
#> F = 1.3227, num df = 1125, denom df = 99, p-value = 0.0764
#> alternative hypothesis: true ratio of concentrations is not equal to 1
#> sample estimates:
#> ratio of concentrations 
#>                 1.32274 
#> 
```
