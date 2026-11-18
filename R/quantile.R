# Local quantile regression (check-loss with kernel weights)

#' Local linear quantile estimator at a point
#'
#' Minimizes the Gaussian-kernel-weighted check loss of a local linear fit
#' around `x` with [stats::optim()] and returns the intercept.
#'
#' @param Y Numeric response vector of length n.
#' @param X Numeric matrix or data frame of covariates (n x d).
#' @param x Numeric vector of length d, the target point.
#' @param tau Quantile level in (0, 1).
#' @param h Bandwidth (scalar).
#' @return Numeric scalar, the estimated conditional `tau`-quantile at `x`.
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

#' Local linear quantile estimates at several points
#'
#' Applies [local_quantile()] to each row of `X0`.
#'
#' @inheritParams local_quantile
#' @param X0 Numeric matrix of target points, one per row.
#' @return Numeric vector of estimated conditional quantiles.
#' @export
local_quantile_vec <- function(Y, X, X0, tau, h) {
  apply(X0, 1, function(x0) local_quantile(Y, X, x0, tau, h))
}
