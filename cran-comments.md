## Resubmission
This is a resubmission. In this version I have:
* Updated http URL to avoid https redirect.

Previous version of the submission added rmarkdown as dependency.
It also addressed concerns by Prof Brian Ripley regarding rgl. rgl package was moved
to Suggest and the package now uses it conditionally.

## Test environments
* local Manjaro, R 4.0.5
* Ubuntu 16.04 (on travis-ci), R 3.4, R-oldrel, R-release and R-devel.
* Win-builder checks: devel, release and oldrelease.

## R CMD check results
There was 1 NOTE (only in Ubuntu systems):
checking installed package size ... NOTE
installed size is  6.6Mb
sub-directories of 1Mb or more:
  libs   6.0Mb

## revdepcheck results
We checked 4 reverse dependencies (3 from CRAN + 1 from Bioconductor), comparing
R CMD check results across CRAN and dev versions of this package.

* We saw 0 new problems
* We failed to check 0 packages

