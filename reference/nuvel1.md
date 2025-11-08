# NUVEL-1 Global model of current plate motions

NNR-NUVEL-1 global model of current plate motions by DeMets et al. 1990

## Usage

``` r
data('nuvel1')
```

## Format

An object of class `data.frame`

- plate.name:

  The rotating plate

- plate.rot:

  The abbreviation of the plate's name

- lat,lon:

  Coordinates of the Pole of Rotation

- angle:

  The amount of rotation (angle in 1 Myr)

- plate.fix:

  The anchored plate, i.e. `plate.rot` moves relative to `plate.fix`

- source:

  Reference to underlying study

## References

DeMets, C., Gordon, R. G., Argus, D. F., & Stein, S. (1990). Current
plate motions. *Geophysical Journal International*, **101**(2), 425-478.
[doi:10.1111/j.1365-246X.1990.tb06579.x](https://doi.org/10.1111/j.1365-246X.1990.tb06579.x)
.

## Examples

``` r
data("nuvel1")
head("nuvel1")
#> [1] "nuvel1"
```
