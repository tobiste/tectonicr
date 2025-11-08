# Apply Rolling Functions using Circular Statistics

A generic function for applying a function to rolling margins of an
array.

## Usage

``` r
roll_circstats(
  x,
  w = NULL,
  FUN,
  axial = TRUE,
  na.rm = TRUE,
  width = NULL,
  by.column = FALSE,
  partial = TRUE,
  fill = NA,
  ...
)
```

## Arguments

- x:

  numeric vector. Values in degrees.

- w:

  (optional) Weights. A vector of positive numbers and of the same
  length as `x`.

- FUN:

  the function to be applied

- axial:

  logical. Whether the data are axial, i.e. pi-periodical (`TRUE`, the
  default) or directional, i.e. \\2 \pi\\-periodical (`FALSE`).

- na.rm:

  logical value indicating whether `NA` values in `x` should be stripped
  before the computation proceeds.

- width:

  integer specifying the window width (in numbers of observations) which
  is aligned to the original sample according to the `align` argument.
  If `NULL`, an optimal width is calculated.

- by.column:

  logical. If `TRUE`, `FUN` is applied to each column separately.

- partial:

  logical or numeric. If `FALSE` then `FUN` is only applied when all
  indexes of the rolling window are within the observed time range. If
  `TRUE` (default), then the subset of indexes that are in range are
  passed to `FUN`. A numeric argument to partial can be used to
  determine the minimal window size for partial computations. See below
  for more details.

- fill:

  a three-component vector or list (recycled otherwise) providing
  filling values at the left/within/to the right of the data range. See
  the fill argument of
  [`zoo::na.fill()`](https://rdrr.io/pkg/zoo/man/na.fill.html) for
  details

- ...:

  optional arguments passed to
  [`zoo::rollapply()`](https://rdrr.io/pkg/zoo/man/rollapply.html)

## Value

numeric vector with the results of the rolling function.

## Note

If the rolling statistics are applied to values that are a function of
distance it is recommended to sort the values first.

## Examples

``` r
data("plates")
plate_boundary <- subset(plates, plates$pair == "na-pa")
data("san_andreas")
PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
distance <- distance_from_pb(
  x = san_andreas,
  PoR = PoR,
  pb = plate_boundary,
  tangential = TRUE
)
dat <- san_andreas[order(distance), ]
roll_circstats(dat$azi, w = 1 / dat$unc, circular_mean, width = 51) |> head()
#> [1] 172.460149 175.401782 178.547220 176.670176 179.709559   1.298709
```
