# Normalize Response and Covariates

This function centers and (optionally) scales the response vector and
each column of the design matrix using the population variance. It is
used to prepare data for Bayesian Lasso regression.

## Usage

``` r
normalize(y, X, scale = TRUE)
```

## Arguments

- y:

  A numeric response vector.

- X:

  A numeric matrix or data frame of covariates (design matrix).

- scale:

  Logical; if `TRUE`, variables are scaled to have unit population
  variance (default is `TRUE`).

## Value

A list with the following elements:

- `vy`: Normalized response vector.

- `mX`: Normalized design matrix.

- `mu.y`: Mean of the response vector.

- `sigma2.y`: Population variance of the response vector.

- `mu.x`: Vector of column means of `X`.

- `sigma2.x`: Vector of population variances for columns of `X`.

## Examples

``` r
set.seed(1)
X <- matrix(rnorm(100 * 10), 100, 10)
beta <- c(2, -3, rep(0, 8))
y <- as.vector(X %*% beta + rnorm(100))
norm_result <- normalize(y, X)
```
