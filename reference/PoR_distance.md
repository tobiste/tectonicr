# Distance to Pole of Rotation

Retrieve the (angular) distance to the PoR (Euler pole).

## Usage

``` r
PoR_distance(x, PoR, FUN = orthodrome)
```

## Arguments

- x:

  `sf` or `data.frame` containing lat and lon coordinates (`lat`, `lon`)

- PoR:

  Pole of Rotation. `"data.frame"` or object of class `"euler.pole"`
  containing the geographical coordinates of the Euler pole

- FUN:

  function to calculate the great-circle distance.
  [`orthodrome()`](https://tobiste.github.io/tectonicr/reference/spherical_angle.md),
  [`haversine()`](https://tobiste.github.io/tectonicr/reference/spherical_angle.md)
  (the default), or
  [`vincenty()`](https://tobiste.github.io/tectonicr/reference/spherical_angle.md).

## Value

numeric vector

## Examples

``` r
data("nuvel1")
por <- subset(nuvel1, nuvel1$plate.rot == "na") # North America relative to Pacific plate
data("san_andreas")

# distance form sf object
PoR_distance(san_andreas, por) |> head()
#> [1] 30.94452 29.49555 30.61975 30.26368 28.31982 26.60286

# distance form data.frame
PoR_distance(sf::st_drop_geometry(san_andreas), por) |> head()
#> [1] 30.94452 29.49555 30.61975 30.26368 28.31982 26.60286
PoR_distance(sf::st_drop_geometry(san_andreas), por, FUN = orthodrome) |> head()
#> [1] 30.94452 29.49555 30.61975 30.26368 28.31982 26.60286
PoR_distance(sf::st_drop_geometry(san_andreas), por, FUN = vincenty) |> head()
#> [1] 30.94452 29.49555 30.61975 30.26368 28.31982 26.60286
```
