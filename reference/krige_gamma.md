# Variogram model of the trace variogram

Fit the variogram gamma \\\gamma\\ values for given distances and
variogram models.

## Usage

``` r
krige_gamma(h, model, nugget, range_km, sill = 1)
```

## Arguments

- h:

  numeric. Distance values (in km) for which the gamma values should be
  fitted.

- model:

  Variogram model: `"spherical"` (default), `"exponential"` or
  `"gaussian"`.

- nugget:

  numeric in \\\[0, 1)\\. Nugget of the trace variogram.

- range_km:

  numeric. Range of the variogram model in km.

- sill:

  numeric in \\(\text{nugget}, 1\]\\. Sill of the trace variogram.

## Value

numeric vector of the fitted gamma \\\gamma\\ values in the variogram.
Same length as `h`

## See also

[`circular_krige()`](https://tobiste.github.io/tectonicr/reference/krige.md)
and
[`variogram_axial()`](https://tobiste.github.io/tectonicr/reference/variogram_axial.md)
