# Circular Distance Matrix Computation

This function computes and returns the distance matrix computed by using
the circular distance measure to compute the distances between the rows
of a data matrix.

## Usage

``` r
circular_dist(x, axial = TRUE, diag = FALSE, upper = FALSE)
```

## Arguments

- x:

  numeric. Angles in degrees.

- axial:

  logical. Whether the data are axial, i.e. pi-periodical (`TRUE`, the
  default) or directional, i.e. \\2 \pi\\-periodical (`FALSE`).

- diag:

  logical value indicating whether the diagonal of the distance matrix
  should be printed by `print.dist`.

- upper:

  logical value indicating whether the upper triangle of the distance
  matrix should be printed by `print.dist`.

## Value

returns an object of class `"dist"`.

## See also

[`circular_distance()`](https://tobiste.github.io/tectonicr/reference/dispersion.md)

## Examples

``` r
# Calculate the distance matrix of the homing data set
data("homing")
homing_dist <- circular_dist(homing, axial = FALSE)
print(homing_dist)
#>              1           2           3           4           5           6
#> 2  0.001902651                                                            
#> 3  0.007596123 0.001902651                                                
#> 4  0.116977778 0.090423978 0.066987298                                    
#> 5  0.146446609 0.116977778 0.090423978 0.001902651                        
#> 6  0.213211782 0.178606195 0.146446609 0.017037087 0.007596123            
#> 7  0.953153894 0.969846310 0.982962913 0.982962913 0.969846310 0.933012702
#> 8  0.883022222 0.909576022 0.933012702 1.000000000 0.998097349 0.982962913
#> 9  0.821393805 0.853553391 0.883022222 0.992403877 0.998097349 0.998097349
#> 10 0.750000000 0.786788218 0.821393805 0.969846310 0.982962913 0.998097349
#>              7           8           9
#> 2                                     
#> 3                                     
#> 4                                     
#> 5                                     
#> 6                                     
#> 7                                     
#> 8  0.017037087                        
#> 9  0.046846106 0.007596123            
#> 10 0.090423978 0.030153690 0.007596123

# Cluster the homing data set based on the distance matrix
homing_cluster <- stats::kmeans(homing_dist, centers = 3)
rose_line(homing, axial = FALSE, col = homing_cluster$cluster, lwd = 2, add = FALSE)
```
