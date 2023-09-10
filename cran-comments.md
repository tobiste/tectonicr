## Test environments
- R-hub windows-x86_64-devel (r-devel)
- R-hub ubuntu-gcc-release (r-release)
- R-hub fedora-clang-devel (r-devel)

## R CMD check results
> On windows-x86_64-devel (r-devel)
  checking CRAN incoming feasibility ... [13s] NOTE
  Maintainer: 'Tobias Stephan <tobias.stephan1@yahoo.com>'
  
  New submission

> On windows-x86_64-devel (r-devel)
  checking examples ... [62s] NOTE
  Examples with CPU (user + system) or elapsed time > 5s
                   user system elapsed
  distance_from_pb 7.14   0.05    7.19
  rolling_test     5.99   0.05    6.03

> On windows-x86_64-devel (r-devel)
  checking HTML version of manual ... [12s] NOTE
  Skipping checking math rendering: package 'V8' unavailable

> On windows-x86_64-devel (r-devel)
  checking for non-standard things in the check directory ... NOTE
  Found the following files/directories:
    ''NULL''

> On windows-x86_64-devel (r-devel)
  checking for detritus in the temp directory ... NOTE
  Found the following files/directories:
    'lastMiKTeXException'

> On ubuntu-gcc-release (r-release)
  checking CRAN incoming feasibility ... [6s/21s] NOTE
  Maintainer: ‘Tobias Stephan <tobias.stephan1@yahoo.com>’
  
  New submission

> On ubuntu-gcc-release (r-release)
  checking examples ... [48s/48s] NOTE
  Examples with CPU (user + system) or elapsed time > 5s
                    user system elapsed
  distance_from_pb 5.639   0.00   5.640
  rolling_test     5.065   0.02   5.086

> On fedora-clang-devel (r-devel)
  checking CRAN incoming feasibility ... [8s/27s] NOTE
  Maintainer: ‘Tobias Stephan <tobias.stephan1@yahoo.com>’
  
  New submission

> On fedora-clang-devel (r-devel)
  checking examples ... [80s/89s] NOTE
  Examples with CPU (user + system) or elapsed time > 5s
                          user system elapsed
  distance_from_pb      11.249  0.024  12.203
  rolling_test           9.161  0.051  10.385
  quick_plot             6.750  0.016   7.572
  projected_pb_strike    5.735  0.004   6.974
  roll_circstats         5.660  0.008   6.528
  por_transformation_df  5.371  0.000   6.060
  PoR_coordinates        5.337  0.032   5.663

> On fedora-clang-devel (r-devel)
  checking HTML version of manual ... NOTE
  Skipping checking HTML validation: no command 'tidy' found
  Skipping checking math rendering: package 'V8' unavailable

0 errors ✔ | 0 warnings ✔ | 10 notes ✖
