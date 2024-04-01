## Test environments
* local Manjaro, R 4.3.3
* macos-latest (release) on Github-Actions
* Ubuntu 22.04.4 (on Github-Actions), R-oldrel, R-release and R-devel.
* Win-builder checks: devel, release and oldrelease.

## R CMD check results
There was 1 NOTE (only in Linux systems):
checking installed package size ... NOTE
  installed size is  7.1Mb
  sub-directories of 1Mb or more:
    libs   6.3Mb


## revdepcheck results
We checked 5 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages
