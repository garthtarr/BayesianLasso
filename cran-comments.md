## R CMD check results

0 errors | 0 warnings | 0 notes

* All local and win-builder checks passed cleanly.


## Test environments
* local Windows 11 install, R 4.4.1
* Ubuntu 22.04 (on GitHub Actions), R-devel, R-release, R-oldrel
* win-builder (devel and release)

## Downstream dependencies
There are currently no downstream dependencies.

## Comments
* Resubmission of version 0.3.6 (BayesianLasso)
* Corrected the `Date` field in the DESCRIPTION file and removed non-standard top-level files (`README.html`, `pic.png`)
* Added compressed logo (`man/figures/logo.png`, <100 KB)
* Replaced deprecated `arma::is_finite(val)` calls with `std::isfinite(val)` to comply with Armadillo ≥ 15  
  (this fix was introduced in 0.3.6, resolving the 0.3.5 NOTE)
* Package passes all checks: 0 errors, 0 warnings, 0 notes on local and win-builder
