# Extract azimuths of line segments

Extract azimuths of line segments

## Usage

``` r
line_azimuth(x, warn = TRUE)

lines_azimuths(x)
```

## Arguments

- x:

  sf object of type `"LINESTRING"` or `"MULTILINESTRING"`

- warn:

  logical; if `TRUE`, warn if `"MULTILINESTRING"` (default).

## Value

sf object of type `"POINT"` with the columns and entries of the first
row of `x`

## Details

It is recommended to perform `line_azimuth()` on single line objects,
i.e. type `"LINESTRING"`, instead of `"MULTILINESTRING"`. This is
because the azimuth of the last point of a line will be calculated to
the first point of the next line otherwise. This will cause a warning
message (if `warn = TRUE`). For `"MULTILINESTRING"` objects, use
`lines_azimuths()`.

## Examples

``` r
data("plates")

# one line:
subset(plates, pair == "af-eu") |>
  smoothr::densify() |>
  line_azimuth() |>
  head()
#> Warning: MULTILINESTRING object is not recommended
#> Simple feature collection with 6 features and 9 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: -5.80797 ymin: 34.00205 xmax: -5.602275 ymax: 34.0191
#> Geodetic CRS:  WGS 84
#>        azi  pair plateA plateB       type displacement             name  nameA
#> 1 95.69929 af-eu     af     eu convergent           in af-eu_convergent Africa
#> 2 95.69906 af-eu     af     eu convergent           in af-eu_convergent Africa
#> 3 95.69883 af-eu     af     eu convergent           in af-eu_convergent Africa
#> 4 95.69861 af-eu     af     eu convergent           in af-eu_convergent Africa
#> 5 95.69838 af-eu     af     eu convergent           in af-eu_convergent Africa
#> 6 95.69815 af-eu     af     eu convergent           in af-eu_convergent Africa
#>     nameB                   geometry
#> 1 Eurasia   POINT (-5.80797 34.0191)
#> 2 Eurasia POINT (-5.766831 34.01569)
#> 3 Eurasia POINT (-5.725692 34.01228)
#> 4 Eurasia POINT (-5.684553 34.00887)
#> 5 Eurasia POINT (-5.643414 34.00546)
#> 6 Eurasia POINT (-5.602275 34.00205)

# multiple lines:
lines_azimuths(plates[1:5, ]) |> head()
#> Simple feature collection with 6 features and 9 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: -0.4379 ymin: -54.8518 xmax: 2.35975 ymax: -54.0374
#> Geodetic CRS:  WGS 84
#>         azi  pair plateA plateB      type displacement            name  nameA
#> 1  52.98880 af-an     af     an divergent          out af-an_divergent Africa
#> 2  51.23433 af-an     af     an divergent          out af-an_divergent Africa
#> 3 141.83793 af-an     af     an divergent          out af-an_divergent Africa
#> 4  44.60702 af-an     af     an divergent          out af-an_divergent Africa
#> 5  47.34529 af-an     af     an divergent          out af-an_divergent Africa
#> 6  45.71276 af-an     af     an divergent          out af-an_divergent Africa
#>        nameB                    geometry
#> 1 Antarctica    POINT (-0.4379 -54.8518)
#> 2 Antarctica POINT (-0.0388257 -54.6772)
#> 3 Antarctica   POINT (0.443182 -54.4512)
#> 4 Antarctica   POINT (0.964534 -54.8322)
#> 5 Antarctica     POINT (1.69481 -54.399)
#> 6 Antarctica    POINT (2.35975 -54.0374)
```
