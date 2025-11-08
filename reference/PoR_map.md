# Map of data in Pole of Rotation reference frame

Transforms the spatial data and azimuths into the PoR reference frame
and shows them in a map

## Usage

``` r
PoR_map(
  x,
  PoR,
  pb = NULL,
  type = c("none", "in", "out", "right", "left"),
  show.deviation = FALSE,
  ...
)
```

## Arguments

- x, pb:

  `sf` objects of the data points and the plate boundary geometries in
  the geographical coordinate system

- PoR:

  Pole of Rotation. `"data.frame"` or object of class `"euler.pole"`
  containing the geographical coordinates of the Pole of Rotation

- type:

  Character. Type of plate boundary (optional). Can be `"out"`, `"in"`,
  `"right"`, or `"left"` for outward, inward, right-lateral, or
  left-lateral moving plate boundaries, respectively. If `"none"` (the
  default), only the PoR-equivalent azimuth is returned.

- show.deviation:

  logical. Whether the data should be color-coded according to the
  deviation from the prediction, or according to the stress regime? Is
  ignored if `type=='none'`.

- ...:

  optional arguments passed to
  [`tectonicr.colors()`](https://tobiste.github.io/tectonicr/reference/tectonicr.colors.md)

## Value

plot

## See also

[`PoR_shmax()`](https://tobiste.github.io/tectonicr/reference/PoR_azi.md),
[`axes()`](https://tobiste.github.io/tectonicr/reference/axes.md),
[`tectonicr.colors()`](https://tobiste.github.io/tectonicr/reference/tectonicr.colors.md)

## Examples

``` r
data("nuvel1")
na_pa <- subset(nuvel1, nuvel1$plate.rot == "na")

data("plates")
plate_boundary <- subset(plates, plates$pair == "na-pa")

data("san_andreas")
PoR_map(san_andreas, PoR = na_pa, pb = plate_boundary, type = "right", show.deviation = TRUE)
```
