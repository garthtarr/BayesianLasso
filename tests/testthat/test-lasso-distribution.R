test_that("qlasso is monotone in p and returns finite values", {
  p <- c(0.01, 0.1, 0.5, 0.9, 0.99)
  x <- qlasso(p, a = 1.2, b = 0.3, c = 0.7)
  expect_length(x, length(p))
  expect_true(all(is.finite(x)))
  expect_true(is.unsorted(x, strictly = FALSE) == FALSE) # nondecreasing
})

test_that("plasso and qlasso approximately invert each other on a grid", {
  set.seed(1)
  x <- seq(-3, 3, length.out = 50)
  p <- plasso(x, a = 1.2, b = 0.3, c = 0.7)
  x2 <- qlasso(p, a = 1.2, b = 0.3, c = 0.7)
  expect_true(max(abs(x - x2)) < 1e-6)
})