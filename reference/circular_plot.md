# Circular plot

Circular plot

## Usage

``` r
circular_plot(
  main = NULL,
  labels = TRUE,
  at = seq(0, 360 - 45, 45),
  cborder = TRUE,
  ...
)
```

## Arguments

- main:

  Character string specifying the title of the plot.

- labels:

  Either a logical value indicating whether to plot labels next to the
  tick marks, or a vector of labels for the tick marks.

- at:

  Optional vector of angles at which tick marks should be plotted. Set
  `at=numeric(0)` to suppress tick marks.

- cborder:

  logical. Border of rose plot.

- ...:

  optional arguments passed to
  [`plot.default()`](https://rdrr.io/r/graphics/plot.default.html)

## Value

none

## Note

Polar diagram where angles increase clockwise.
