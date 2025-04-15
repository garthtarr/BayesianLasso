#' @title The Lasso Distribution
#' @name LassoDistribution
#' @aliases zlasso dlasso plasso qlasso rlasso elasso vlasso mlasso
#' @description
#' Normalizing constant Z, probability density function, cumulative distribution function,
#' quantile function, and random generation for the Lasso distribution with parameters \code{a}, \code{b}, and \code{c}.
#' In addition, Mills ratio, mean, and variance of the Lasso distribution are provided.
#' 
#' @usage
#' zlasso(a, b, c, logarithm)
#' dlasso(x, a, b, c, logarithm)
#' plasso(q, a, b, c)
#' qlasso(p, a, b, c)
#' rlasso(n, a, b, c)
#' elasso(a, b, c)
#' vlasso(a, b, c)
#' mlasso(a, b, c)
#'
#' @details
#' If \eqn{X \sim \text{Lasso}(a, b, c)} then its density function is:
#' \deqn{
#' p(x;a,b,c) = Z^{-1} \exp\left(-\frac{1}{2} a x^2 + bx - c|x| \right)
#' }
#' where \eqn{x \in \mathbb{R}}, \eqn{a > 0}, \eqn{b \in \mathbb{R}}, \eqn{c > 0}, and \eqn{Z} is the normalizing constant.
#'
#' More details are included for the CDF, quantile function, and normalizing constant in the original documentation.
#'
#' @param x,q Vector of quantiles (vectorized).
#' @param p Vector of probabilities.
#' @param a  Vector of precision parameter which must be non-negative.
#' @param b  Vector of off set parameter.
#' @param c  Vector of tuning parameter which must be non-negative values.
#' @param n Number of observations.
#' @param logarithm Logical. If \code{TRUE}, probabilities are returned on the log scale.
#'
#' @return
#' \itemize{
#'   \item \code{zlasso}: the normalizing constant Z.
#'   \item \code{dlasso}: density function.
#'   \item \code{plasso}: cumulative distribution function.
#'   \item \code{qlasso}: quantile function.
#'   \item \code{rlasso}: random generation.
#'   \item \code{elasso}: expected value.
#'   \item \code{vlasso}: variance.
#'   \item \code{mlasso}: mode.
#' }
#'
#' @examples
#' a <- 2; b <- 1; c <- 3
#' x <- seq(-3, 3, length.out = 1000)
#' plot(x, dlasso(x, a, b, c, logarithm = FALSE), type = 'l')
#'
#' r <- rlasso(1000, a, b, c)
#' hist(r, breaks = 50, probability = TRUE, col = "grey", border = "white")
#' lines(x, dlasso(x, a, b, c, logarithm = FALSE), col = "blue")
#'
#' plasso(0, a, b, c)
#' qlasso(0.25, a, b, c)
#' elasso(a, b, c)
#' vlasso(a, b, c)
#' mlasso(a, b, c)
#' 
#' 
#' 
#' # The Modified_Hans_Gibbs() function uses the Lasso distribution to draw 
#' # samples from the full conditional distribution of the regression coefficients.
#' 
#' 
#' 
#' y <- 1:20
#' X <- matrix(c(1:20,12:31,7:26),20,3,byrow = T)
#' 
#' a1 <- b1 <- u1 <- v1 <- 0.01
#' sigma2_init <- 1
#' lambda_init <- 0.1
#' beta_init <- rep(1, ncol(X))
#' nsamples <- 1000
#' verbose <- 100
#' 
#' Output <- Modified_Hans_Gibbs(
#'   X, y, a1, b1, u1, v1,
#'   nsamples, beta_init, lambda_init, sigma2_init, verbose
#' )
#' 
#' colMeans(Output$mBeta)
#' mean(Output$vlambda2)
#' 
#' @export zlasso
#' @export dlasso
#' @export plasso
#' @export qlasso
#' @export rlasso
#' @export elasso
#' @export vlasso
#' @export mlasso
NULL