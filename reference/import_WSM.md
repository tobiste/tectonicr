# World Stress Map Database (WSM)

Download WSM2025 or WSM2016 database from the GFZ sever and applies
optional filters. If `destdir` is specified, the data can be reloaded in
a later R session using `load_WSM()` using the same arguments.

## Usage

``` r
download_WSM(
  destdir = tempdir(),
  load = TRUE,
  version = c("2025", "2016"),
  ...
)

load_WSM(
  file,
  quality = c("A", "B", "C", "D", "E", "X"),
  lat_range = c(-90, 90),
  lon_range = c(-180, 180),
  depth_range = c(-Inf, Inf),
  type = c("BO", "BOC", "BOT", "BS", "DIF", "FMA", "FMF", "FMS", "GFI", "GFM", "GFS",
    "GVA", "HF", "HFG", "HFM", "HFH", "HFP", "HFS", "OC", "PC", "SWB", "SWL", "SWS"),
  regime = c("N", "NS", "T", "TS", "S", NA)
)

download_WSM2016(destdir = tempdir(), load = TRUE, ...)

load_WSM2016(
  file,
  quality = c("A", "B", "C", "D", "E"),
  lat_range = c(-90, 90),
  lon_range = c(-180, 180),
  depth_range = c(-Inf, Inf),
  type = c("BO", "BOC", "BOT", "BS", "DIF", "FMA", "FMF", "FMS", "GFI", "GFM", "GFS",
    "GVA", "HF", "HFG", "HFM", "HFP", "OC", "PC", "SWB", "SWL", "SWS"),
  regime = c("N", "NS", "T", "TS", "S", NA)
)
```

## Source

<https://datapub.gfz.de/download/10.5880.WSM.2025.001-Scbwez/WSM_Database_2025.csv>

<https://datapub.gfz-potsdam.de/download/10.5880.WSM.2016.001/wsm2016.csv>

## Arguments

- destdir:

  where to save files, defaults to
  [`base::tempdir()`](https://rdrr.io/r/base/tempfile.html),
  [`base::getwd()`](https://rdrr.io/r/base/getwd.html) is also possible.

- load:

  `TRUE` load the dataset into R, `FALSE` return the file name of the
  downloaded object.

- version:

  character. Version of the World stress map database. Either `"2025"`
  (default) or `"2016"`

- ...:

  (optional) arguments passed to `load_WSM()`

- file:

  the name of the file which the data are to be read from.

- quality:

  a character vectors containing the quality levels to be included.
  Includes all quality ranks (A-X) by default.

- lat_range, lon_range:

  two-element numeric vectors giving the range of latitudes and
  longitudes (in degrees).

- depth_range:

  two-element numeric vectors giving the depth interval (in km)

- type:

  a character vectors containing the methods of stress inversion to be
  included. Includes all methods by default. See WSM2016 manual for used
  acronyms.

- regime:

  a character vectors containing the stress regimes to be included.
  Acronyms: `"N"` - normal, `"T"` - thrust, `"S"` - strike-slip,
  `"NS"` - oblique normal, `"TS"` - oblique thrust, and `NA` - unknown
  faulting

## Value

`sf` object of and the parsed numeric uncertainty (`unc`) based on the
reported standard deviation and the quality rank. If `load=FALSE`, the
path to the downloaded file is returned.

## Note

Because of R-compatibility and easy readability reasons, the downloaded
dataset is a modified version of the original, WSM server version: All
column names have been changed from uppercase (in the original dataset)
to lowercase characters. Unknown azimuth values are represented by `NA`
values instead of `999` in the original. Unknown regimes are represented
by `NA` instead of "U" in the original.

## References

Heidbach, O., M. Rajabi, X. Cui, K. Fuchs, B. M\<U+00FC\>ller, J.
Reinecker, K. Reiter, M. Tingay, F. Wenzel, F. Xie, M. O. Ziegler, M.-L.
Zoback, and M. D. Zoback (2018): The World Stress Map database release
2016: Crustal stress pattern across scales. *Tectonophysics*, **744**,
484-498,
[doi:10.1016/j.tecto.2018.07.007](https://doi.org/10.1016/j.tecto.2018.07.007)
.

Heidbach, Oliver; Rajabi, Mojtaba; Di Giacomo, Domenico; Harris, James;
Lammers, Steffi; Morawietz, Sophia; Pierdominici, Simona; Reiter,
Karsten; von Specht, Sebastian; Storchak, Dmitry; Ziegler, Moritz O.
(2025): World Stress Map Database Release 2025. GFZ Data Services.
[doi:10.5880/WSM.2025.001](https://doi.org/10.5880/WSM.2025.001)

## Examples

``` r
if (FALSE) { # \dontrun{
download_WSM(
  quality = c("A", "B", "C"), lat_range = c(51, 72),
  lon_range = c(-180, -130), depth_range = c(0, 10), type = "FMS"
)
} # }
```
