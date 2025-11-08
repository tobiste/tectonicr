# Conversion between PoR to geographical coordinate reference system of raster data

Helper function to transform raster data set from PoR to geographical
coordinates

## Usage

``` r
geographical_to_PoR_raster(x, PoR)

PoR_to_geographical_raster(x, PoR)
```

## Arguments

- x:

  `"SpatRaster"` or `"RasterLayer"`

- PoR:

  Pole of Rotation. `"data.frame"` or object of class `"euler.pole"`
  containing the geographical coordinates of the Euler pole

## Value

terra "SpatRaster" object
