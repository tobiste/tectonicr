# Coordinates of the Pole of Rotation Reference System

Retrieve the PoR equivalent coordinates of an object

## Usage

``` r
PoR_coordinates(x, PoR)
```

## Arguments

- x:

  Can be either a `"data.frame"` containing `lat` and `lon` coordinates
  of a point in the geographical CRS or the `lat.PoR`, `lon.PoR`) of the
  point in the PoR CRS, a two-column matrix containing the lat and lon
  coordinates, a `sf` object, or a `raster` object.

- PoR:

  Pole of Rotation. `"data.frame"` containing the geographical
  coordinates of the Euler pole

## Value

`PoR_coordinates()` returns `data.frame` with the PoR coordinates
(`lat.PoR`, `lon.PoR`).

## Examples

``` r
data("nuvel1")
por <- subset(nuvel1, nuvel1$plate.rot == "na") # North America relative to Pacific plate
data("san_andreas")

# coordinates from sf object
san_andreas.por_sf <- PoR_coordinates(san_andreas, por)
head(san_andreas.por_sf)
#>     lon.PoR  lat.PoR
#> 1 -85.46364 59.05548
#> 2 -78.16683 60.50445
#> 3 -84.55317 59.38025
#> 4 -85.74431 59.73632
#> 5 -84.31053 61.68018
#> 6 -80.60932 63.39714

# coordinates from data.frame
san_andreas.por_df <- PoR_coordinates(sf::st_drop_geometry(san_andreas), por)
head(san_andreas.por_df)
#>    lat.PoR   lon.PoR
#> 1 59.05548 -85.46364
#> 2 60.50445 -78.16683
#> 3 59.38025 -84.55317
#> 4 59.73632 -85.74431
#> 5 61.68018 -84.31053
#> 6 63.39714 -80.60932
```
