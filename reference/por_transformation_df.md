# Conversion between spherical PoR to geographical coordinate system of data.frames

Transformation from spherical PoR to geographical coordinate system and
vice versa

## Usage

``` r
geographical_to_PoR_df(x, PoR)

PoR_to_geographical_df(x, PoR)
```

## Arguments

- x:

  `"data.frame"` containing `lat` and `lon` coordinates of a point in
  the geographical CRS or the `lat.PoR`, `lon.PoR`) of the point in the
  PoR CRS.

- PoR:

  Pole of Rotation. `"data.frame"` containing the geographical
  coordinates of the Euler pole

## Value

`"data.frame"` with the transformed coordinates (`lat.PoR` and `lon.PoR`
for PoR CRS, or `lat` and `lon` for geographical CRS).
