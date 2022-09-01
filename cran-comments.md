## Test environments
- R-hub windows-x86_64-devel (r-devel)
- R-hub ubuntu-gcc-release (r-release)
- R-hub fedora-clang-devel (r-devel)

## R CMD check results
> On windows-x86_64-devel (r-devel)
  checking CRAN incoming feasibility ... [48s] NOTE
  Maintainer: 'Tobias Stephan <tobias.stephan1@yahoo.com>'
  
  New submission
  
  Found the following (possibly) invalid URLs:
    URL: https://doi.org/10.1007/978-3-030-26050-7_435-1
      From: man/relative_rotation.Rd
      Status: 404
      Message: Not Found
    URL: https://doi.org/10.1029/2001GC000252
      From: man/pb2002.Rd
            man/plates.Rd
      Status: 503
      Message: Service Unavailable
    URL: https://doi.org/10.1029/97JB03390
      From: man/norm_chisq.Rd
      Status: 503
      Message: Service Unavailable

> On windows-x86_64-devel (r-devel)
  checking for detritus in the temp directory ... NOTE
  Found the following files/directories:
    'lastMiKTeXException'

> On ubuntu-gcc-release (r-release)
  checking CRAN incoming feasibility ... NOTE
  Maintainer: ‘Tobias Stephan <tobias.stephan1@yahoo.com>’
  
  New submission
  
  Found the following (possibly) invalid URLs:
    URL: https://doi.org/10.1007/978-3-030-26050-7_435-1
      From: man/relative_rotation.Rd
      Status: 404
      Message: Not Found
    URL: https://doi.org/10.1029/2001GC000252
      From: man/pb2002.Rd
            man/plates.Rd
      Status: 503
      Message: Service Unavailable
    URL: https://doi.org/10.1029/97JB03390
      From: man/norm_chisq.Rd
      Status: 503
      Message: Service Unavailable

> On fedora-clang-devel (r-devel)
  checking CRAN incoming feasibility ... [7s/49s] NOTE
  Maintainer: ‘Tobias Stephan <tobias.stephan1@yahoo.com>’
  
  New submission
  
  Found the following (possibly) invalid URLs:
    URL: https://doi.org/10.1007/978-3-030-26050-7_435-1
      From: man/relative_rotation.Rd
      Status: 404
      Message: Not Found
    URL: https://doi.org/10.1029/2001GC000252
      From: man/pb2002.Rd
            man/plates.Rd
      Status: 503
      Message: Service Unavailable
    URL: https://doi.org/10.1029/97JB03390
      From: man/norm_chisq.Rd
      Status: 503
      Message: Service Unavailable

> On fedora-clang-devel (r-devel)
  checking examples ... [29s/58s] NOTE
  Examples with CPU (user + system) or elapsed time > 5s
               user system elapsed
  stress_paths  2.4  0.028   5.219

> On fedora-clang-devel (r-devel)
  checking HTML version of manual ... NOTE
  Skipping checking HTML validation: no command 'tidy' found
  Skipping checking math rendering: package 'V8' unavailable
