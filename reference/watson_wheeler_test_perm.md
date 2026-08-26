# Watson-Wheeler Test of Homogeneity of Means

Performs the Watson-Wheeler test for homogeneity on two or more samples
of circular data.

## Usage

``` r
watson_wheeler_test_perm(x, y, axial = TRUE, n_perm = 1000L)
```

## Arguments

- x, y:

  numeric vectors. Angles in degrees

- axial:

  logical. Whether the data are axial, i.e. \\\pi\\-periodical (`TRUE`,
  the default) or directional, i.e. \\2 \pi\\-periodical (`FALSE`).

- n_perm:

  integer. Number of permutations

## Value

list

## Details

The Watson-Wheeler (or Mardia-Watson-Wheeler, or uniform score) test is
a non-parametric test to compare two or several samples. The difference
between the samples can be in either the mean or the variance.

The p-value is estimated by assuming that the test statistic follows a
chi-squared distribution. For this approximation to be valid, all groups
must have at least 10 elements.

## See also

Other Tests:
[`ar_test()`](https://tobiste.github.io/tectonicr/reference/ar_test.md),
[`kuiper_test()`](https://tobiste.github.io/tectonicr/reference/kuiper_test.md),
[`norm_chisq()`](https://tobiste.github.io/tectonicr/reference/norm_chisq.md),
[`rayleigh-test`](https://tobiste.github.io/tectonicr/reference/rayleigh-test.md),
[`watson_test()`](https://tobiste.github.io/tectonicr/reference/watson_test.md),
[`watson_two_sample`](https://tobiste.github.io/tectonicr/reference/watson_two_sample.md),
[`weighted-rayleigh-test`](https://tobiste.github.io/tectonicr/reference/weighted-rayleigh-test.md)

## Examples

``` r
set.seed(20250411)
x1 <- c(35, 45, 50, 55, 60, 70, 85, 95, 105, 120)
x2 <- c(75, 80, 90, 100, 110, 130, 135, 140, 150, 160, 165)
watson_wheeler_test_perm(x1, x2, axial = FALSE)
#> $statistic
#> [1] 3.67827
#> 
#> $p.value
#> [1] 0.1688312
#> 

data1 <- rvm(n=20, mean = 0, kappa=3)
data2 <- rvm(n=20, mean = 90, kappa=2)
watson_wheeler_test_perm(data1, data2, axial = FALSE)
#> $statistic
#> [1] 15.96367
#> 
#> $p.value
#> [1] 0.000999001
#> 

# San Andreas Fault Data:
data(san_andreas)
data("nuvel1")
PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
sa.por <- PoR_shmax(san_andreas, PoR, "right")
watson_wheeler_test_perm(sa.por$azi.PoR, rvm(100, 135, 10))
#> $statistic
#> [1] 2.934737
#> 
#> $p.value
#> [1] 0.2547453
#> 
```
