# Kuiper Test of Circular Uniformity

Kuiper's test statistic is a rotation-invariant Kolmogorov-type test
statistic. The critical values of a modified Kuiper's test statistic are
used according to the tabulation given in Stephens (1970).

## Usage

``` r
kuiper_test(x, alpha = 0, axial = TRUE, quiet = FALSE)
```

## Arguments

- x:

  numeric vector containing the circular data which are expressed in
  degrees

- alpha:

  Significance level of the test. Valid levels are `0.01`, `0.05`, and
  `0.1`. This argument may be omitted (`NULL`, the default), in which
  case, a range for the p-value will be returned.

- axial:

  logical. Whether the data are axial, i.e. \\\pi\\-periodical (`TRUE`,
  the default) or circular, i.e. \\2 \pi\\-periodical (`FALSE`).

- quiet:

  logical. Prints the test's decision.

## Value

list containing the test statistic `statistic` and the significance
level `p.value`.

## Details

If `statistic > p.value`, the null hypothesis is rejected. If not,
randomness (uniform distribution) cannot be excluded.

## Examples

``` r
# Example data from Mardia and Jupp (1999), pp. 93
pidgeon_homing <- c(55, 60, 65, 95, 100, 110, 260, 275, 285, 295)
kuiper_test(pidgeon_homing, alpha = .05)
#> Reject Null Hypothesis
#> $statistic
#> [1] 2.262115
#> 
#> $p.value
#> [1] 1.747
#> 

# San Andreas Fault Data:
data(san_andreas)
data("nuvel1")
PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
sa.por <- PoR_shmax(san_andreas, PoR, "right")
kuiper_test(sa.por$azi.PoR, alpha = .05)
#> Reject Null Hypothesis
#> $statistic
#> [1] 16.60463
#> 
#> $p.value
#> [1] 1.747
#> 
```
