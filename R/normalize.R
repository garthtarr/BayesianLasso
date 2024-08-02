
normalize <- function(y, X) {
  n <- length(y)
  p <- ncol(X)
  
  mu.y <- mean(y)
  sigma2.y <- (n - 1) * var(y) / n
  vy <- (y - mu.y) / sqrt(sigma2.y)
  
  # Normalise covariates
  mX <- matrix(0, n, p)
  mu.x <- c()
  sigma2.x <- c()
  for (j in 1:p)
  {
    mu.x[j] <- mean(X[, j])
    sigma2.x[j] <- (n - 1) * var(X[, j]) / n
    mX[, j] <- (X[, j] - mu.x[j]) / sqrt(sigma2.x[j])
  }
  
  return(list(vy = vy, mX = mX, mu.y = mu.y, sigma2.y = sigma2.y, mu.x = mu.x, sigma2.x = sigma2.x))
}