# Wallraff Test of Angular Distances

Performs the Wallraff test of angular distances or angular dispersion
around the mean.

## Usage

``` r
walraff_test(samples, ref = NULL, axial = TRUE)
```

## Arguments

- samples:

  list of samples of orientations

- ref:

  a vector of angles (in degrees) used as reference to compute the
  angular distances from, in each group. It should contain either one
  element (a common reference for all groups) or as many elements as
  there are groups, *in the same order* as the elements of the list
  `samples`.

  If `ref` is `NULL` (the default), the mean angle of each group will be
  used as reference. In this situation, the Wallraff test becomes a
  comparison of angular dispersion around the mean

- axial:

  logical. Whether the data are axial, i.e. \\\pi\\-periodical (`TRUE`,
  the default) or directional, i.e. \\2 \pi\\-periodical (`FALSE`). In
  case of axial data, the angles will be doubled for the test.

## Value

A list with class `"htest"` containing the following components:

- statistic:

  the chi-squared statistic from
  [`kruskal.test`](https://rdrr.io/r/stats/kruskal.test.html).

- parameter:

  the degrees of freedom for the chi-squared statistic.

- p.value:

  the p-value for the test.

- method:

  a character string containing the name of the test.

- data.name:

  a character string giving the name(s) of the data.

## Details

The Wallraff test of angular distances between two or more groups is
performed and the results are printed. The null hypothesis is that
distances are equal across groups.

The test proceeds by computing the angular distances from a reference
angle, in each group. The angular distance between two angles is the
circular range. Then the distances are compared with a usual rank sum
test (Kruskal-Wallis,
[stats::kruskal.test](https://rdrr.io/r/stats/kruskal.test.html)). When
there are only two groups, the Wilcoxon-Mann-Whitney test could be used
but [stats::wilcox.test](https://rdrr.io/r/stats/wilcox.test.html)
without continuity correction for the p-value is equivalent to
[stats::kruskal.test](https://rdrr.io/r/stats/kruskal.test.html) so only
[stats::kruskal.test](https://rdrr.io/r/stats/kruskal.test.html) is used
here.

The Wallraff test is most frequently used to compare angular dispersion
around the mean, between samples. In this case, the reference angle is
the mean angle of each sample. This is the default here, when no
reference angles are provided.

## References

Batschelet, E (1981). Circular Statistics in Biology. chap. 6.10, p. 124

Zar, J H (2010). Biostatistical analysis. sec. 27.7-8, p. 643

## See also

Other Two-Sample-Tests:
[`ar_test()`](https://tobiste.github.io/tectonicr/reference/ar_test.md),
[`equal_kappa_test()`](https://tobiste.github.io/tectonicr/reference/equal_kappa_test.md),
[`uniform_scores_test()`](https://tobiste.github.io/tectonicr/reference/uniform_scores_test.md),
[`watson_two_sample`](https://tobiste.github.io/tectonicr/reference/watson_two_sample.md),
[`watson_wheeler_test_perm()`](https://tobiste.github.io/tectonicr/reference/watson_wheeler_test_perm.md),
[`watson_williams_test()`](https://tobiste.github.io/tectonicr/reference/watson_williams_test.md)

## Examples

``` r
set.seed(20250411)
x1 <- c(35, 45, 50, 55, 60, 70, 85, 95, 105, 120)
x2 <- c(75, 80, 90, 100, 110, 130, 135, 140, 150, 160, 165)
walraff_test(list(x1, x2))
#> 
#>  Wallraff rank sum test of angular distance to the sample mean
#>  directions
#> 
#> data:  1 and 2
#> Kruskal-Wallis chi-squared = 0.40165, df = 1, p-value = 0.5262
#> 

x1 <- rvm(100, mean = 90, kappa = 20)
x2 <- rvm(100, mean = 90, kappa = 10)
x3 <- rvm(100, mean = 90, kappa = 1)
walraff_test(list(x1, x2, x3))
#> 
#>  Wallraff rank sum test of angular distance to the sample mean
#>  directions
#> 
#> data:  1 and 2 and 3
#> Kruskal-Wallis chi-squared = 116.51, df = 2, p-value < 2.2e-16
#> 

# San Andreas Fault Data:
data(san_andreas)
data("nuvel1")
PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
sa.por <- PoR_shmax(san_andreas, PoR, "right")
walraff_test(list(sa.por$azi.PoR, rvm(100, 135, 10)))
#> 
#>  Wallraff rank sum test of angular distance to the sample mean
#>  directions
#> 
#> data:  1 and 2
#> Kruskal-Wallis chi-squared = 0.082398, df = 1, p-value = 0.7741
#> 
```
