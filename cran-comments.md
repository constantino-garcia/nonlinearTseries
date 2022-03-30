* This submission fixes the following warning on gcc-12:
"void operator delete(void*)’ called on pointer returned from a mismatched allocation function [-Wmismatched-new-delete]" 

## Test environments
* local Manjaro, R 4.1.2
* Ubuntu 16.04 (on travis-ci), R-oldrel, R-release and R-devel.
* Win-builder checks: devel, release and oldrelease.
* Docker image rocker/r-edge:gcc-12

## R CMD check results
There was 1 NOTE (only in Ubuntu systems):
checking installed package size ... NOTE
installed size is  6.6Mb
sub-directories of 1Mb or more:
  libs   6.0Mb


