## R CMD check results

0 errors | 0 warnings | 1 note

* There was 1 NOTE about the package size which is due to compiled code (7.2 MB); no large data files are included.
* This is a new release.

## Test environments
* local Windows 11 install, R 4.4.1
* Ubuntu 22.04 (on GitHub Actions), R-devel, R-release, R-oldrel
* win-builder (devel and release)

## Downstream dependencies
There are currently no downstream dependencies.

## Comments
* Replaced deprecated `arma::is_finite(val)` calls with `std::isfinite(val)` to comply with Armadillo 15.
* Bumped version from 0.3.5 to 0.3.6.
