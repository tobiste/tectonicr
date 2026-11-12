# Theoretical Plate Tectonic Stress Paths

Construct \\\sigma\_{Hmax}\\ lines that are following small circles,
great circles, or loxodromes of an Euler pole for the relative plate
motion.

## Usage

``` r
eulerpole_paths(x, type = c("sc", "gc", "ld"), n = 10L, angle = 45, cw)

eulerpole_smallcircles(x, n = 10L)

eulerpole_greatcircles(x, n = 10L)

eulerpole_loxodromes(x, n = 10L, angle = 45, cw)
```

## Arguments

- x:

  Either an object of class `"euler.pole"` or `"data.frame"` containing
  coordinates of Euler pole in lat, lon, and rotation angle (optional).

- type:

  character. String specifying the type of curves to export. Either
  `"sm"` for small circles (default), `"gc"` for great circles, or
  `"ld"` for loxodromes.

- n:

  integer. Number of equally spaced curves; `n = 10` by default (angular
  distance between curves: `180 / n`)

- angle:

  numeric. Direction of loxodromes; `angle = 45` by default.

- cw:

  logical. Sense of loxodromes: `TRUE` for clockwise loxodromes
  (left-lateral displaced plate boundaries). `FALSE` for
  counterclockwise loxodromes (right-lateral displaced plate
  boundaries).

## Value

`sf` object

## Details

Maximum horizontal stress can be aligned to three types of curves
related to relative plate motion:

- Small circles:

  Lines that have a constant distance to the Euler pole. If `x` contains
  `angle`, output additionally gives absolute velocity on small circle
  (°/Myr -\> km/Myr).

- Great circles:

  Paths of the shortest distance between the Euler pole and its
  antipodal position.

- Loxodromes:

  Lines of constant bearing, i.e. curves cutting small circles at a
  constant angle.

## Author

Tobias Stephan

## Examples

``` r
data("nuvel1")
por <- subset(nuvel1, nuvel1$plate.rot == "na") # North America relative to
# Pacific plate

eulerpole_smallcircles(por)
#> Simple feature collection with 11 features and 2 fields
#> Geometry type: MULTILINESTRING
#> Dimension:     XY
#> Bounding box:  xmin: -180 ymin: -84.709 xmax: 180 ymax: 84.709
#> Geodetic CRS:  WGS 84
#> # A tibble: 11 × 3
#>                                                           geometry     d abs_vel
#>  *                                           <MULTILINESTRING [°]> <dbl>   <dbl>
#>  1 ((101.833 -48.709, 101.833 -48.709, 101.833 -48.709, 101.833 -…     0     0  
#>  2 ((101.833 -30.709, 102.4799 -30.71571, 103.1265 -30.73582, 103…   -18    26.9
#>  3 ((101.833 -12.709, 102.9175 -12.72024, 104.0016 -12.75396, 105…   -36    51.2
#>  4 ((101.833 5.291, 103.2953 5.275842, 104.757 5.230387, 106.2172…   -54    70.4
#>  5 ((101.833 23.291, 103.6965 23.27168, 105.5586 23.21377, 107.41…   -72    82.8
#>  6 ((101.833 41.291, 104.228 41.26618, 106.6194 41.19178, 109.003…    90    87.1
#>  7 ((101.833 59.291, 105.1832 59.25627, 108.521 59.15234, 111.834…    72    82.8
#>  8 ((101.833 77.291, 108.4309 77.22258, 114.9037 77.0195, 121.141…    54    70.4
#>  9 ((-78.167 84.709, -89.46284 84.59139, -99.79473 84.25303, -108…    36    51.2
#> 10 ((-78.167 66.709, -79.57283 66.69442, -80.97315 66.65076, -82.…    18    26.9
#> 11 ((-78.167 48.709, -78.167 48.709, -78.167 48.709, -78.167 48.7…     0     0  
eulerpole_greatcircles(por)
#> Simple feature collection with 11 features and 1 field
#> Geometry type: MULTILINESTRING
#> Dimension:     XY
#> Bounding box:  xmin: -180 ymin: -89.991 xmax: 180 ymax: 89.991
#> Geodetic CRS:  WGS 84
#> # A tibble: 11 × 2
#>        d                                                                geometry
#>  * <dbl>                                                   <MULTILINESTRING [°]>
#>  1     0 ((101.833 -48.709, 101.833 -48.609, 101.833 -48.509, 101.833 -48.409, …
#>  2    36 ((101.833 -48.709, 101.7441 -48.62806, 101.6554 -48.54706, 101.5671 -4…
#>  3    72 ((101.833 -48.709, 101.689 -48.67801, 101.5451 -48.64684, 101.4014 -48…
#>  4   108 ((101.833 -48.709, 101.6888 -48.73981, 101.5444 -48.77044, 101.3998 -4…
#>  5   144 ((101.833 -48.709, 101.7438 -48.78987, 101.6543 -48.87067, 101.5645 -4…
#>  6     0 ((101.833 -48.709, 101.833 -48.809, 101.833 -48.909, 101.833 -49.009, …
#>  7    36 ((101.833 -48.709, 101.9222 -48.78987, 102.0117 -48.87067, 102.1015 -4…
#>  8    72 ((101.833 -48.709, 101.9772 -48.73981, 102.1216 -48.77044, 102.2662 -4…
#>  9   108 ((101.833 -48.709, 101.977 -48.67801, 102.1209 -48.64684, 102.2646 -48…
#> 10   144 ((101.833 -48.709, 101.9219 -48.62806, 102.0106 -48.54706, 102.0989 -4…
#> 11     0 ((101.833 -48.709, 101.833 -48.609, 101.833 -48.509, 101.833 -48.409, …
eulerpole_loxodromes(x = por, angle = 45, n = 10, cw = FALSE)
#> Simple feature collection with 15 features and 1 field
#> Geometry type: MULTILINESTRING
#> Dimension:     XY
#> Bounding box:  xmin: -180 ymin: -82.82221 xmax: 180 ymax: 82.82221
#> Geodetic CRS:  WGS 84
#> # A tibble: 15 × 2
#>        d                                                                geometry
#>  * <dbl>                                                   <MULTILINESTRING [°]>
#>  1   108 ((-78.07054 66.58404, -78.01693 66.51323, -77.96404 66.44238, -77.9118…
#>  2   144 ((101.6024 77.35336, 101.3391 77.42378, 101.0733 77.49398, 100.805 77.…
#>  3     0 ((101.833 41.291, 101.7388 41.36167, 101.6444 41.43227, 101.5498 41.50…
#>  4    36 ((101.3044 5.933699, 101.2458 6.003953, 101.1871 6.074159, 101.1283 6.…
#>  5    72 ((101.6185 -30.1275, 101.5917 -30.05698, 101.5648 -29.98648, 101.5378 …
#>  6   108 ((101.7399 -48.5214, 101.7063 -48.4543, 101.6725 -48.38724, 101.6385 -…
#>  7   144 ((101.5909 -48.59325, 101.5042 -48.55193, 101.4174 -48.5107, 101.3307 …
#>  8     0 ((101.5338 -48.70929, 101.4266 -48.70954, 101.3195 -48.70986, 101.2123…
#>  9    36 ((101.591 -48.82535, 101.5043 -48.86716, 101.4177 -48.90905, 101.3311 …
#> 10    72 ((101.7412 -48.89695, 101.7086 -48.96434, 101.6762 -49.03178, 101.644 …
#> 11   108 ((101.9268 -48.89653, 101.961 -48.96356, 101.9956 -49.03054, 102.0305 …
#> 12   144 ((102.0762 -48.82424, 102.1639 -48.86512, 102.2519 -48.9058, 102.3401 …
#> 13     0 ((102.1322 -48.70793, 102.2393 -48.70703, 102.3465 -48.70586, 102.4536…
#> 14    36 ((102.0739 -48.59215, 102.1596 -48.54991, 102.245 -48.50747, 102.3301 …
#> 15    72 ((101.9241 -48.52098, 101.9561 -48.45352, 101.9878 -48.38601, 102.0191…
eulerpole_loxodromes(x = por, angle = 30, cw = TRUE)
#> Simple feature collection with 13 features and 1 field
#> Geometry type: MULTILINESTRING
#> Dimension:     XY
#> Bounding box:  xmin: -180 ymin: -85.21953 xmax: 180 ymax: 85.21953
#> Geodetic CRS:  WGS 84
#> # A tibble: 13 × 2
#>        d                                                                geometry
#>  * <dbl>                                                   <MULTILINESTRING [°]>
#>  1   144 ((101.5141 -47.93824, 101.4813 -47.85449, 101.4487 -47.77068, 101.4165…
#>  2     0 ((100.8943 -48.20753, 100.7944 -48.15215, 100.6949 -48.09657, 100.5958…
#>  3    36 ((100.6237 -48.6678, 100.4929 -48.66181, 100.3621 -48.65551, 100.2314 …
#>  4    72 ((100.8139 -49.14575, 100.7018 -49.19187, 100.5894 -49.23775, 100.4766…
#>  5   108 ((101.4023 -49.45697, 101.3538 -49.53763, 101.3048 -49.61821, 101.2555…
#>  6   144 ((102.1618 -49.47885, 102.1969 -49.5624, 102.2319 -49.64599, 102.2667 …
#>  7     0 ((102.7905 -49.20269, 102.8946 -49.25631, 102.9989 -49.30995, 103.1032…
#>  8    36 ((103.044 -48.73753, 103.1752 -48.74063, 103.3065 -48.74373, 103.4377 …
#>  9    72 ((102.8344 -48.26342, 102.9424 -48.21527, 103.0504 -48.16714, 103.1584…
#> 10   108 ((102.251 -47.95948, 102.2967 -47.87845, 102.3424 -47.79746, 102.3882 …
#> 11   144 ((101.833 -21.06283, 101.8579 -20.97623, 101.883 -20.88965, 101.9081 -…
#> 12     0 ((101.833 41.291, 101.8996 41.37758, 101.9664 41.46413, 102.0334 41.55…
#> 13    36 ((-78.167 76.35517, -78.26446 76.26854, -78.36015 76.18185, -78.45412 …
eulerpole_smallcircles(data.frame(lat = 30, lon = 10))
#> Simple feature collection with 11 features and 1 field
#> Geometry type: MULTILINESTRING
#> Dimension:     XY
#> Bounding box:  xmin: -180 ymin: -84 xmax: 180 ymax: 84
#> Geodetic CRS:  WGS 84
#> # A tibble: 11 × 2
#>                                                                   geometry     d
#>  *                                                   <MULTILINESTRING [°]> <dbl>
#>  1 ((-170 -30, -170 -30, -170 -30, -170 -30, -170 -30, -170 -30, -170 -30…     0
#>  2 ((-170 -12, -169.4314 -12.00774, -168.8632 -12.03093, -168.2959 -12.06…   -18
#>  3 ((-170 6, -168.9363 5.985529, -167.8735 5.942134, -166.8123 5.869864, …   -36
#>  4 ((-170 24, -168.4063 23.97832, -166.8145 23.91332, -165.2266 23.80513,…   -54
#>  5 ((-170 42, -167.6973 41.96867, -165.3999 41.87482, -163.1132 41.7188, …   -72
#>  6 ((-170 60, -166.4035 59.95107, -162.8282 59.80475, -159.2943 59.56246,…    90
#>  7 ((-170 78, -161.8144 77.88851, -153.9057 77.56006, -146.4959 77.03116,…    72
#>  8 ((10 84, -3.638991 83.81339, -15.74781 83.28465, -25.60532 82.4858, -3…    54
#>  9 ((10 66, 7.401912 65.96464, 4.822514 65.85889, 2.279844 65.68372, -0.2…    36
#> 10 ((10 48, 9.169016 47.98869, 8.339768 47.9548, 7.51398 47.8984, 6.69335…    18
#> 11 ((10 30, 10 30, 10 30, 10 30, 10 30, 10 30, 10 30, 10 30, 10 30, 10 30…     0
```
