# Normalized Chi-Squared Test for Circular Data

A quantitative comparison between the predicted and observed directions
of \\\sigma\_{Hmax}\\ is obtained by the calculation of the average
azimuth and by a normalized \\\chi^2\\ test.

## Usage

``` r
norm_chisq(obs, prd, unc)
```

## Arguments

- obs:

  Numeric vector containing the observed azimuth of \\\sigma\_{Hmax}\\,
  same length as `prd`

- prd:

  Numeric vector containing the modeled azimuths of \\\sigma\_{Hmax}\\,
  i.e. the return object from
  [`model_shmax()`](https://tobiste.github.io/tectonicr/reference/model_shmax.md)

- unc:

  Uncertainty of observed \\\sigma\_{Hmax}\\, either a numeric vector or
  a number

## Value

Numeric vector

## Details

The normalized \\\chi^2\\ test is \$\$ {Norm} \chi^2_i = = \frac{
\sum^M\_{i = 1} \left( \frac{\alpha_i - \alpha\_{{predict}}}{\sigma_i}
\right) ^2} {\sum^M\_{i = 1} \left( \frac{90}{\sigma_i} \right) ^2 }\$\$
The value of the chi-squared test statistic is a number between 0 and 1
indicating the quality of the predicted \\\sigma\_{Hmax}\\ directions.
Low values (\\\le 0.15\\) indicate good agreement, high values (\\\>
0.7\\) indicate a systematic misfit between predicted and observed
\\\sigma\_{Hmax}\\ directions.

## References

Wdowinski, S., 1998, A theory of intraplate tectonics. *Journal of
Geophysical Research: Solid Earth*, **103**, 5037-5059, doi:
10.1029/97JB03390.

## Examples

``` r
data("nuvel1")
PoR <- subset(nuvel1, nuvel1$plate.rot == "na") # North America relative to
# Pacific plate
data(san_andreas)
point <- data.frame(lat = 45, lon = 20)
prd <- model_shmax(point, PoR)
norm_chisq(obs = c(50, 40, 42), prd = prd$sc, unc = c(10, NA, 5))
#> [1] 0.001426232

data(san_andreas)
prd2 <- PoR_shmax(san_andreas, PoR, type = "right")
norm_chisq(obs = prd2$azi.PoR, 135, unc = san_andreas$unc)
#> [1] 0.0529064
```
