# Plate Stress Dummy Grid

Helper functions to create a dummy grid for small circles, great
circles, and loxodromes of an Euler pole

## Usage

``` r
smallcircle_dummy(n)

greatcircle_dummy(n)

loxodrome_dummy(n, angle, cw)
```

## Arguments

- n:

  Number of curves

- angle:

  Direction of loxodromes (in degree)

- cw:

  logical. Sense of loxodromes: `TRUE` for clockwise loxodromes
  (right-lateral displaced plate boundaries). `FALSE` for
  counterclockwise loxodromes (left-lateral displaced plate boundaries).

## Value

`data.frame`
