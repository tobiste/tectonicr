# Product of quaternions

Helper function for multiplication of two quaternions. Concatenation of
two rotations R1 followed by R2

## Usage

``` r
product_Q4(q1, q2, normalize = FALSE)
```

## Arguments

- q1, q2:

  two objects of class `"quaternion"`. first rotation R1 expressed by q1
  followed by second rotation R2 expressed by q2

- normalize:

  logical. Whether a quaternion normalization should be applied (TRUE)
  or not (FALSE, the default).

## Value

object of class `"quaternion"`

## Note

Multiplication is not commutative.
