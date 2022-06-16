## Test environments
- R-hub windows-x86_64-devel (r-devel)
- R-hub ubuntu-gcc-release (r-release)
- R-hub fedora-clang-devel (r-devel)

## R CMD check results
There were no ERRORs or WARNINGs. 

There were 3 NOTEs:
  
* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Tobias Stephan <tobias.stephan1@yahoo.com>'
  
  New submission

* checking examples ... NOTE
  Examples with CPU (user + system) or elapsed time > 5s
               user system elapsed
  stress_paths 8.12   0.26    8.39

* checking for detritus in the temp directory ... NOTE
  Found the following files/directories:
    'lastMiKTeXException'
