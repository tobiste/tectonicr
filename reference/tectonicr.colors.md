# Colors for input variables

assigns colors to continuous or categorical values for plotting

## Usage

``` r
tectonicr.colors(
  x,
  n = 10,
  pal = NULL,
  categorical = FALSE,
  na.value = "grey",
  ...
)
```

## Arguments

- x:

  values for color assignment

- n:

  integer. number of colors for continuous colors (i.e. \`categorical =
  FALSE“).

- pal:

  either a named vector specifying the colors for categorical values, or
  a color function. If `NULL`, default colors are
  [`RColorBrewer::brewer.pal()`](https://rdrr.io/pkg/RColorBrewer/man/ColorBrewer.html)
  (`categorical = TRUE`) and
  [`viridis::viridis()`](https://sjmgarnier.github.io/viridisLite/reference/viridis.html)
  (`categorical = FALSE`).

- categorical:

  logical.

- na.value:

  color for `NA` values (categorical).

- ...:

  optional arguments passed to palette function

## Value

named color vector

## Examples

``` r
val1 <- c("N", "S", "T", "T", NA)
tectonicr.colors(val1, categorical = TRUE)
#>         N         S         T         T      <NA> 
#> "#66C2A5" "#FC8D62" "#8DA0CB" "#8DA0CB"    "grey" 
tectonicr.colors(val1, pal = stress_colors(), categorical = TRUE)
#>         N         S         T         T      <NA> 
#> "#D55E00" "#009E73" "#0072B2" "#0072B2"    "grey" 

val2 <- runif(10)
tectonicr.colors(val2, n = 5)
#>   (0.6,0.7]   (0.3,0.4]   (0.3,0.4]   [0.1,0.2]   (0.6,0.7]   [0.1,0.2] 
#> "#27AD81FF" "#3B528BFF" "#3B528BFF" "#440154FF" "#27AD81FF" "#440154FF" 
#>     (0.9,1]   (0.8,0.9]     (0.9,1]   (0.7,0.8] 
#> "#FDE725FF" "#AADC32FF" "#FDE725FF" "#5DC863FF" 
```
