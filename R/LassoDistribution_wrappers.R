#' @export
dlasso <- function(x, a, b, c, logarithm = FALSE) {
  val <- .Call(`_BayesianLasso_dlasso`, as.numeric(x), a, b, c, logarithm)
  if (length(x) == 1L) val <- as.vector(val)[1] else val <- as.vector(val)
  val
}

#' @export
plasso <- function(q, a, b, c) {
  val <- .Call(`_BayesianLasso_plasso`, as.numeric(q), a, b, c)
  if (length(q) == 1L) val <- as.vector(val)[1] else val <- as.vector(val)
  val
}

#' @export
qlasso <- function(p, a, b, c) {
  val <- .Call(`_BayesianLasso_qlasso`, as.numeric(p), a, b, c)
  if (length(p) == 1L) val <- as.vector(val)[1] else val <- as.vector(val)
  val
}

#' @export
elasso <- function(a, b, c) {
  as.vector(.Call(`_BayesianLasso_elasso`, a, b, c))[1]
}

#' @export
vlasso <- function(a, b, c) {
  as.vector(.Call(`_BayesianLasso_vlasso`, a, b, c))[1]
}

#' @export
mlasso <- function(a, b, c) {
  as.vector(.Call(`_BayesianLasso_mlasso`, a, b, c))[1]
}

#' @export
rlasso <- function(n, a, b, c) {
  val <- .Call(`_BayesianLasso_rlasso`, as.integer(n), a, b, c)
  as.vector(val)
}
