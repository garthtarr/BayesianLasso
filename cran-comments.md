## R CMD check results

0 errors | 0 warnings | 1 note

## Test environments
* local Windows 11 install, R 4.4.2
* Ubuntu 22.04 (on GitHub Actions), R-devel, R-release, R-oldrel
* win-builder (release)

## Downstream dependencies
There are currently no downstream dependencies.

## Comments
This is a resubmission of BayesianLasso version 0.3.9.

The previous CRAN version was archived because it depended on the archived package `RcppClock`.
That dependency has now been completely removed.

I rebuilt the package from a clean source tree and verified that the submitted source tarball contains no reference to `RcppClock`.

Additional updates in this resubmission include:
* corrected the `Date` field in `DESCRIPTION`
* removed non-standard top-level files
* added a compressed logo (`man/figures/logo.png`, <100 KB)
* replaced deprecated `arma::is_finite(val)` calls with `std::isfinite(val)`
* removed the non-CRAN dependency `bayeslm` from `DESCRIPTION`

Checks performed:
* `devtools::check(args = "--as-cran")`
* win-builder (release)

All checks pass with only minor NOTES.