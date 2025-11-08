# Numerical values to World Stress Map Quality Ranking

Assigns numeric values of the precision (sd.) of each measurement to the
categorical quality ranking of the World Stress Map (A, B, C, D, E, X).

## Usage

``` r
parse_wsm_quality(x)

quantise_wsm_quality(x)
```

## Arguments

- x:

  Either a string or a character/factor vector of WSM quality ranking

## Value

`"numeric"`. the standard deviation of stress azimuth

## References

Heidbach, O., Barth, A., M\<U+00FC\>ller, B., Reinecker, J.,
Stephansson, O., Tingay, M., Zang, A. (2016). WSM quality ranking
scheme, database description and analysis guidelines for stress
indicator. *World Stress Map Technical Report* **16-01**, GFZ German
Research Centre for Geosciences.
[doi:10.2312/wsm.2016.001](https://doi.org/10.2312/wsm.2016.001)

## Examples

``` r
parse_wsm_quality(c("A", "B", "C", "D", NA, "E", "X"))
#>    A    B    C    D <NA>    E    X 
#>   15   20   25   40   NA   90  180 
data("san_andreas")
head(parse_wsm_quality(san_andreas$quality))
#>  C  C  C  C  C  C 
#> 25 25 25 25 25 25 
```
