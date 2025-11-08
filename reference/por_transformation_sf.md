# Conversion between PoR to geographical coordinates of sf data

Transform spatial objects from PoR to geographical coordinate reference
system and vice versa.

## Usage

``` r
PoR_to_geographical_sf(x, PoR)

geographical_to_PoR_sf(x, PoR)
```

## Arguments

- x:

  `sf`, `SpatRast`, or `Raster*` object of the data points in
  geographical or PoR coordinate system

- PoR:

  Pole of Rotation. `"data.frame"` or object of class `"euler.pole"`
  containing the geographical coordinates of the Euler pole

## Value

`sf` or `SpatRast` object of the data points in the transformed
geographical or PoR coordinate system

## Details

The PoR coordinate reference system is oblique transformation of the
geographical coordinate system with the Euler pole coordinates being the
translation factors.
