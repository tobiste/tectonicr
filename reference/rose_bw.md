# Selecting optimal number of bins and width for rose diagrams

\#' @rdname rose_bw rose_binwidth \<- function(n, axial = TRUE, ...) if
(isTRUE(axial)) r \<- 180 else r \<- 360r / rose_bins(n)

## Usage

``` r
optimal_rose_bins(
  n,
  rule = c("sturges", "rice"),
  full_circle = TRUE,
  min_bins = 4L
)
```

## Arguments

- n:

  Integer. number of data

- rule:

  character. Rule to select optimal bin width, one of "sturges" (the
  default) or "rice". The optimal number according to Sturges' rule is
  \\1 + \log2(n)\\, while for Rice's rule it is \\2 n^{1/3}\\ (less
  conservative for large `n`).

- full_circle:

  logical. Whether bin width should be calculated for the full (`TRUE`,
  the default) or the half circle (`FALSE`).

- min_bins:

  integer. The minimum number of bins. `4` by default.
