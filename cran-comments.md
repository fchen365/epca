## Test environments
* local R installation, R 4.0.0
* ubuntu 16.04 (on travis-ci), R 4.0.0
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 0 note

* This is a new release.

## Changes since submission on May 20, 2020 

* Rephrase the description field in the DESCRIPTION file.
* Add small executable examples
* Replace all T and F by TRUE and FALSE.
* Check \value field in .Rd files for all exported methods.
* Avoid writing to anywhere apart from the R session’s temporary directory.
