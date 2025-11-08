# Weighting Factors

Helper function to transform uncertainty angles into weighting factors

## Usage

``` r
weighting(
  x,
  method = c("linear-inverse", "inverse", "cosine", "none"),
  max.err = 90
)
```

## Arguments

- x:

  numeric. Uncertainty angle in degrees.

- method:

  character. One of `"linear-inverse"` (the default), `"inverse"`,
  `"cosine"`, or `"none"` (no transformation).

- max.err:

  numeric. The maximum expected error for x (90 by default).

## Value

numeric

## Details

Linear inverse: \\w = 1 - x/\sigma\\, where \\\sigma\\ is the maximum
error expected for \\x\\ (e.g. \\90^\circ\\).

Inverse: \\w = 1/x\\

Cosine: \\w = \cos{x}\\

## Examples

``` r
x <- seq(0, 90, 1)

plot(x, weighting(x, "inverse"), col = 1, type = "l",
xlab = "Uncertainty angle in degrees", ylab = "weight")
lines(x, weighting(x, "cosine"), col = 2)
lines(x, weighting(x, "linear-inverse"), col = 3)
legend("topright", col = 1:3, lty = 1,
legend = c("inverse", "cosine", "linear-inverse"))
```
