# Local quantile regression (check-loss with kernel weights)

#' Local quantile estimator at a point
#' @param Y numeric vector, response
#' @param X numeric matrix/data.frame, predictors (n x d)
#' @param x numeric vector of length d, target point
#' @param tau quantile level in (0,1)
#' @param h bandwidth (scalar or length-d vector; here treated as scalar)
#' @return numeric scalar, estimated conditional tau-quantile at x
#' @export
#' @importFrom matrixStats rowProds
#' @importFrom stats optim
local_quantile <- function(Y, X, x, tau, h) {
  if (is.vector(X)) X <- matrix(X, ncol = 1)
  n <- length(Y)
  d <- ncol(X)
  x <- as.vector(x)
  U <- sweep(X, 2, x, "-")  # n x d
  weights <- matrixStats::rowProds(kernel_gau(U / h))

  loss <- function(b) {
    residual <- Y - b[1] - U %*% b[-1]
    sum((tau - (residual < 0)) * residual * weights)
  }

  init_b <- rep(0, d + 1)
  opt <- stats::optim(init_b, loss)
  opt$par[1]
}

#' Vectorized local quantile for multiple target points
#' @export
local_quantile_vec <- function(Y, X, X0, tau, h) {
  apply(X0, 1, function(x0) local_quantile(Y, X, x0, tau, h))
}
