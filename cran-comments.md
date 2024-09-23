## Test environments
* local Endeavour OS, R 4.4.1
* macos-latest (release) on Github-Actions
* Ubuntu-latest (on Github-Actions), R-oldrel, R-release and R-devel.
* Windows-latest (on Github-Actions).
* Win-builder checks: r-devel.

## R CMD check results
There was 1 NOTE (only in Linux systems):
checking installed package size ... NOTE
  installed size is  8.0Mb
  sub-directories of 1Mb or more:
    libs   7.2Mb


## revdepcheck results
We checked 4 strong reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages
