# A better Bayesian Lasso

This package will grow over time. It currently includes functions for the 
Lasso Distribution.

The corresponding paper is [NOT YET] available on [arXiv]().

## Installation

The `BayesianLasso` package can be installed from GitHub using

```
remotes::install_github("garthtarr/BayesianLasso")
```

## Usage
zlasso(a, b, c, logarithm)
dlasso(x, a, b, c, logarithm)
plasso(q, a, b, c)
qlasso(p, a, b, c)
rlasso(n, a, b, c)
elasso(a, b, c)
vlasso(a, b, c)
mlasso(a, b, c)
MillsRatio(d)

### Lasso distribution
Normalizing constant Z, probability density function, cumulative distribution function,
quantile function, and random generation for the Lasso distribution with parameters \code{a}, \code{b}, and \code{c}.
In addition, Mills ratio, mean, and variance of the Lasso distribution are provided.

If \eqn{X \sim \text{Lasso}(a, b, c)} then its density function is:
\deqn{
p(x;a,b,c) = Z^{-1} \exp\left(-\frac{1}{2} a x^2 + bx - c|x| \right)
}
where \eqn{x \in \mathbb{R}}, \eqn{a > 0}, \eqn{b \in \mathbb{R}}, \eqn{c > 0}, and \eqn{Z} is the normalizing constant.

More details are included for the CDF, quantile function, and normalizing constant in the original documentation.

### Bayesian Lasso

To see how this can be used to fit Bayesian Lasso models see the vignette.
