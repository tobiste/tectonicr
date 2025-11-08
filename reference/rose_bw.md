# Selecting optimal number of bins and width for rose diagrams

Selecting optimal number of bins and width for rose diagrams

## Usage

``` r
rose_bins(n, round = FALSE)

rose_binwidth(n, axial = TRUE, ...)
```

## Arguments

- n:

  Integer. number of data

- round:

  Logical. Whether bin width is round to zero digits (`round=TRUE`, the
  default) or as is (`FALSE`).

- axial:

  Logical. Whether data are uniaxial (`axial=FALSE`) or biaxial (`TRUE`,
  the default).

- ...:

  Additional arguments passed to `rose_bw()`.
