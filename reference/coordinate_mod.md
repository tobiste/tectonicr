# Coordinate Correction

Corrects the longitudes or latitudes to value between -180&deg and
180&deg or -90&deg and 90°

## Usage

``` r
longitude_modulo(x)

latitude_modulo(x)
```

## Arguments

- x:

  Longitude(s) or latitude(s) in degrees

## Value

numeric

## Examples

``` r
longitude_modulo(-361 + 5 * 360) # -1
#> [1] -1
latitude_modulo(-91 + 5 * 180) # 89
#> [1] 89
```
