test_that("Modified_Hans_Gibbs returns expected components and dimensions", {
  set.seed(1)
  n <- 30; p <- 10
  X <- matrix(rnorm(n*p), n, p)
  y <- rnorm(n)
  
  fit <- Modified_Hans_Gibbs(
    y = y, X = X,
    beta_init = rep(0, p),
    a1 = 0.01, b1 = 0.01, u1 = 0.01, v1 = 0.01,
    nsamples = 50,
    lambda_init = 1, sigma2_init = 1,
    thin = 1, 
    verbose = 0,
    tune_lambda2 = TRUE,
    rao_blackwellization = FALSE
  )
  
  expect_true(is.list(fit))
  expect_true(all(c("mBeta","vsigma2","vlambda2") %in% names(fit)))
  expect_equal(ncol(fit$mBeta), p)
  expect_equal(nrow(fit$mBeta), 50)
  expect_equal(length(fit$vsigma2), 50)
  expect_equal(length(fit$vlambda2), 50)
})