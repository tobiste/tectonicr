# PoR coordinate reference system

Create the reference system transformed in Euler pole coordinate

## Usage

``` r
PoR_crs(x)
```

## Arguments

- x:

  `"data.frame"` or `"euler.pole"` object containing the geographical
  coordinates of the Euler pole

## Value

Object of class `crs`

## Details

The PoR coordinate reference system is oblique transformation of the
geographical coordinate system with the Euler pole coordinates being the
the translation factors.

## See also

[`sf::st_crs()`](https://r-spatial.github.io/sf/reference/st_crs.html)

## Examples

``` r
data("nuvel1")
por <- subset(nuvel1, nuvel1$plate.rot == "na") # North America relative to Pacific plate
PoR_crs(por)
#> Coordinate Reference System:
#>   User input: +proj=ob_tran +o_proj=longlat +datum=WGS84 +o_lat_p=48.709 +o_lon_p=-78.167 
#>   wkt:
#> GEOGCRS["unnamed",
#>     BASEGEOGCRS["unknown",
#>         DATUM["World Geodetic System 1984",
#>             ELLIPSOID["WGS 84",6378137,298.257223563,
#>                 LENGTHUNIT["metre",1]],
#>             ID["EPSG",6326]],
#>         PRIMEM["Greenwich",0,
#>             ANGLEUNIT["degree",0.0174532925199433],
#>             ID["EPSG",8901]]],
#>     DERIVINGCONVERSION["unknown",
#>         METHOD["PROJ ob_tran o_proj=longlat"],
#>         PARAMETER["o_lat_p",48.709,
#>             ANGLEUNIT["degree",0.0174532925199433,
#>                 ID["EPSG",9122]]],
#>         PARAMETER["o_lon_p",-78.167,
#>             ANGLEUNIT["degree",0.0174532925199433,
#>                 ID["EPSG",9122]]]],
#>     CS[ellipsoidal,2],
#>         AXIS["longitude",east,
#>             ORDER[1],
#>             ANGLEUNIT["degree",0.0174532925199433,
#>                 ID["EPSG",9122]]],
#>         AXIS["latitude",north,
#>             ORDER[2],
#>             ANGLEUNIT["degree",0.0174532925199433,
#>                 ID["EPSG",9122]]]]
```
