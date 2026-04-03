## R CMD check results

0 errors | 0 warnings | 1 note

* checking for future file timestamps: unable to verify current time
  This is a local environment issue and not a package problem.

* This is a resubmission. Changes since last submission:
  - Fixed the vignette build issue: prebuilt vignette outputs are now 
    included in inst/doc/ as required.
  - Expanded the vignette with a complete simulated workflow example, 
    convergence diagnostics (trace plots, ACF plots, R-hat, mixing 
    statistics), and interpretation guidance, addressing reviewer comments.
  - Removed real dataset dependencies (lars, Ecdat) from the vignette 
    to ensure CRAN build reproducibility.