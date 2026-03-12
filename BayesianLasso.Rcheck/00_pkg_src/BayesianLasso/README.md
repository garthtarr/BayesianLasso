
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BayesianLasso <img src="man/figures/logo.png" align="right" width="195"/>

<!-- badges: start -->

<!-- badges: end -->

BayesianLasso is an R package for efficient Bayesian inference in sparse
linear regression models using the Bayesian Lasso. It includes optimized
Gibbs sampling algorithms and utilities for working with the Lasso
distribution.

## Installation

You can install the development version of BayesianLasso from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("garthtarr/BayesianLasso")
```

## Features

- Efficient Gibbs samplers for Bayesian Lasso (e.g.,
  Modified_Hans_Gibbs, Modified_PC_Gibbs)

- Support for drawing from the Lasso distribution

- Utilities for computing moments and densities

## Example Usage

These are basic examples which show you how to solve a common problem:

``` r
library(BayesianLasso)
## basic example code

# Simulated data
set.seed(123)
X <- matrix(rnorm(100), 20, 5)
y <- rnorm(20)
beta_init <- rep(1, 5)

# Run modified Hans Gibbs sampler
result <- Modified_Hans_Gibbs(
  X = X,
  y = y,
  beta_init = beta_init,
  a1 = 0.01,
  b1 = 0.01,
  u1 = 0.01,
  v1 = 0.01,
  nsamples = 100,
  lambda_init = 0.1,
  sigma2_init = 1,
  verbose = 0, tune_lambda2 = TRUE, rao_blackwellization = FALSE
)

str(result)
#> List of 6
#>  $ mBeta   : num [1:100, 1:5] 0.33112 0.06944 0.00723 -0.00512 0.31311 ...
#>  $ vsigma2 : num [1:100, 1] 0.787 0.457 0.597 1.219 0.735 ...
#>  $ vlambda2: num [1:100, 1] 64.45 76.32 40.19 41.73 5.21 ...
#>  $ mA      : num[0 , 0 ] 
#>  $ mB      : num[0 , 0 ] 
#>  $ mC      : num[0 , 0 ]
```

The `Modified_Hans_Gibbs()` function returns a list with the following
components:

- `mBeta`: MCMC samples of the regression coefficients
  $\boldsymbol{\beta}$, stored as a matrix with `nsamples` rows and `p`
  columns.
- `vsigma2`: MCMC samples of the error variance $\sigma^2$.
- `vlambda2`: MCMC samples of the global shrinkage parameter
  $\lambda^2$.
- `mA`, `mB`, `mC`: Matrices containing the MCMC samples of the Lasso
  distribution parameters $A_j$, $B_j$, and $C_j$ for each coefficient
  $\beta_j$, where each row corresponds to one MCMC iteration and each
  column to a regression coefficient.

## Lasso Distribution Functions

The package provides functions for working with the Lasso distribution:

- zlasso(): Normalizing constant

- dlasso(): Density function

- plasso(): CDF

- qlasso(): Quantile function

- rlasso(): Random generation

- elasso(): Expected value

- vlasso(): Variance

- mlasso(): Mode

- MillsRatio(): Mills ratio

## Citation

If you use this package in your work, please cite it appropriately.
Citation information can be found using:

``` r
citation("BayesianLasso")
#> To cite package 'BayesianLasso' in publications use:
#> 
#>   Ormerod J, Davoudabadi M, Tarr G, Mueller S, Tidswell J (2025).
#>   _Bayesian Lasso Regression and Tools for the Lasso Distribution_. R
#>   package version 0.3.0, <https://garthtarr.github.io/BayesianLasso/>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {Bayesian Lasso Regression and Tools for the Lasso Distribution},
#>     author = {John Ormerod and Mohammad Javad Davoudabadi and Garth Tarr and Samuel Mueller and Jonathon Tidswell},
#>     year = {2025},
#>     note = {R package version 0.3.0},
#>     url = {https://garthtarr.github.io/BayesianLasso/},
#>   }
```
