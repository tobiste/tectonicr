# Equivalent rotation

Transforms a sequence of rotations into a new reference system

## Usage

``` r
equivalent_rotation(x, fixed, rot)
```

## Arguments

- x:

  Object of class `"data.frame"` containing the Euler poles of plate
  rotations:

  `plate.rot`

  :   Moving plate

  `lat`, `lon`

  :   coordinates of Euler pole

  `angle`

  :   Angle of rotation

  `plate.fix`

  :   Fixed plate

- fixed:

  plate that will be regarded as fixed. Has to be one out of
  `x$plate.fix`

- rot:

  (optional) plate that will be regarded as rotating. Has to be one out
  of `x$plate.rot`.

## Value

sequence of plate rotations in new reference system. Same object class
as `x`

## See also

[`relative_rotation()`](https://tobiste.github.io/tectonicr/reference/relative_rotation.md)

## Examples

``` r
data(nuvel1) # load the NUVEL1 rotation parameters

# all nuvel1 rotation equivalent to fixed Africa:
equivalent_rotation(nuvel1, fixed = "af")
#>    plate.rot        lat        lon     angle plate.fix
#> af        af  90.000000    0.00000 0.0000000        af
#> an        an  -5.729623  141.14853 0.1330674        af
#> ar        ar  24.358850   24.41733 0.4163188        af
#> au        au  12.575099   50.25300 0.6604014        af
#> ca        ca -64.951860 -164.64172 0.1551761        af
#> co        co  17.709052 -121.08039 1.3658002        af
#> eu        eu -21.224307  159.74610 0.1283940        af
#> in        in  23.801567   28.95037 0.4299189        af
#> nz        nz  43.352648 -113.62807 0.4948982        af
#> na        na -79.039887 -140.84695 0.2489005        af
#> sa        sa -62.620310  140.61897 0.3237372        af
#> jf        jf -36.503136   69.97118 0.8851710        af
#> ph        ph -54.222891  -24.76214 1.0236334        af
#> pa        pa -59.160000  106.82600 0.9695000        af
# relative plate motion between Eurasia and India:
equivalent_rotation(nuvel1, "eu", "in") # lat = 24.58, lon = 18.07, angle = 0.528
#>    plate.rot      lat      lon     angle plate.fix
#> in        in 24.57682 18.07207 0.5281698        eu
```
