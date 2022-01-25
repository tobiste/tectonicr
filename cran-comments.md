## R CMD check results
There were 1 ERRORs, 1 WARNING, and 3 NOTEs

### ERRORS:

* checking PDF version of manual without hyperrefs or index ... ERROR

### WARNINGS:

* checking PDF version of manual ... WARNING
LaTeX errors when creating PDF version.
This typically indicates Rd problems.
LaTeX errors found:
! LaTeX Error: Environment bmatrix undefined.

See the LaTeX manual or LaTeX Companion for explanation.
Type  H <return>  for immediate help.
 ...                                              
! Misplaced alignment tab character &.
<argument> ..._{x}^{2}\left (1-\cos \psi \right )&
                                                  u_{x}u_{y}\left (1-\cos \p...
l.752 ...eft(1-\cos \psi \right)\end{bmatrix}} }{}
                                                  
! Misplaced alignment tab character &.
<argument> ...-\cos \psi \right )-u_{z}\sin \psi &
                                                  u_{x}u_{z}\left (1-\cos \p...
l.752 ...eft(1-\cos \psi \right)\end{bmatrix}} }{}
                                                  
! Misplaced alignment tab character &.
<argument> ...-\cos \psi \right )+u_{z}\sin \psi &
                                                  \cos \psi +u_{y}^{2}\left ...
l.752 ...eft(1-\cos \psi \right)\end{bmatrix}} }{}
                                                  
! Misplaced alignment tab character &.
<argument> ..._{y}^{2}\left (1-\cos \psi \right )&
                                                  u_{y}u_{z}\left (1-\cos \p...
l.752 ...eft(1-\cos \psi \right)\end{bmatrix}} }{}
                                                  
! Misplaced alignment tab character &.
<argument> ...-\cos \psi \right )-u_{y}\sin \psi &
                                                  u_{z}u_{y}\left (1-\cos \p...
l.752 ...eft(1-\cos \psi \right)\end{bmatrix}} }{}
                                                  
! Misplaced alignment tab character &.
<argument> ...-\cos \psi \right )+u_{x}\sin \psi &
                                                  \cos \psi +u_{z}^{2}\left ...
l.752 ...eft(1-\cos \psi \right)\end{bmatrix}} }{}
                                                  
! LaTeX Error: \begin{list} on input line 748 ended by \end{bmatrix}.

See the LaTeX manual or LaTeX Companion for explanation.
Type  H <return>  for immediate help.
 ...                                              

### NOTES:

* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Tobias Stephan <tobias.stephan1@yahoo.com>’

New submission

Version contains large components (0.0.0.9000)

Possibly mis-spelled words in DESCRIPTION:
  Wdowinski (10:43)

The Title field should be in title case. Current version is:
‘Modeling the direction of the maximum horizontal stress using relative plate motion’
In title case that is:
‘Modeling the Direction of the Maximum Horizontal Stress using Relative Plate Motion’

The Description field should not start with the package name,
  'This package' or similar.

* checking R code for possible problems ... NOTE
eulerpole_smallcircles: no visible binding for global variable
  ‘small_circle’
Undefined global functions or variables:
  small_circle

* checking for non-standard things in the check directory ... NOTE
Found the following files/directories:
  ‘PlateTectonicStressR-manual.tex’
