# Kernel density and conditional density utilities

#' Kernel density estimate of the covariate density
#'
#' Product Gaussian kernel estimate of \eqn{f_X(x)}.
#'
#' @param X Numeric matrix or data frame of covariates (n x d). A vector is
#'   treated as a single covariate.
#' @param x Numeric vector of length d, the evaluation point.
#' @param h Bandwidth (scalar).
#' @return Numeric scalar, the estimated density at `x`.
#' @export
fx_den <- function(X, x, h) {
  if (is.vector(X)) X <- matrix(X, ncol = 1)
  d <- ncol(X)
  x <- as.vector(x)
  kernel_products <- matrixStats::rowProds(kernel_gau(sweep(X, 2, x, "-") / h))
  mean(kernel_products) / (h^d)
}

#' Joint kernel density of (Y, X) at (Q, x)
#'
#' @inheritParams fqx
#' @return Numeric scalar.
#' @keywords internal
fq <- function(Y, X, x, Q, h) {
  if (is.vector(X)) X <- matrix(X, ncol = 1)
  d <- ncol(X)
  x <- as.vector(x)
  kernel_products <- matrixStats::rowProds(kernel_gau(sweep(X, 2, x, "-") / h))
  kernel_q <- kernel_ep((Y - Q) / h)
  joint_prod <- kernel_products * kernel_q
  mean(joint_prod) / (h^(d + 1))
}

#' Conditional density of Y given X via a kernel ratio
#'
#' Estimates \eqn{f_{Y|X}(Q \mid x)} as the ratio of a joint kernel density
#' (Epanechnikov kernel in Y, product Gaussian kernel in X) to the covariate
#' kernel density.
#'
#' @param Y Numeric response vector of length n.
#' @param X Numeric matrix or data frame of covariates (n x d). A vector is
#'   treated as a single covariate.
#' @param x Numeric vector of length d, the conditioning point.
#' @param Q Numeric scalar, the value of Y at which the density is evaluated.
#' @param h Bandwidth (scalar).
#' @return Numeric scalar, the estimated conditional density.
#' @seealso [fqx2()] for the indicator-based version used by the estimator.
#' @export
fqx <- function(Y, X, x, Q, h) {
  if (is.vector(X)) X <- matrix(X, ncol = 1)
  x <- as.vector(x)
  kernel_products <- matrixStats::rowProds(kernel_gau(sweep(X, 2, x, "-") / h))
  f_hat_x <- sum(kernel_products)
  kernel_q <- kernel_ep((Y - Q) / h)
  joint_prod <- kernel_products * kernel_q
  f_hat_q <- sum(joint_prod)
  f_hat_q / (h * f_hat_x)
}

#' Conditional density of Y given X via an indicator approximation
#'
#' Estimates \eqn{f_{Y|X}(Q \mid x)} by the kernel-weighted share of
#' observations with \eqn{|Y - Q| < h}, divided by \eqn{2h}. This is the
#' density estimate used inside [local_eq_estimate()].
#'
#' @inheritParams fqx
#' @return Numeric scalar, the estimated conditional density.
#' @export
fqx2 <- function(Y, X, x, Q, h) {
  if (is.vector(X)) X <- matrix(X, ncol = 1)
  x <- as.vector(x)
  kernel_products <- matrixStats::rowProds(kernel_gau(sweep(X, 2, x, "-") / h))
  indicator <- abs(Y - Q) < h
  num <- sum(kernel_products * indicator)
  den <- sum(kernel_products)
  num / (den * 2 * h)
}

#' Vectorized conditional density estimates
#'
#' Applies [fqx2()] row by row, pairing the i-th row of `x_mat` with the i-th
#' element of `Q_vec`.
#'
#' @inheritParams fqx
#' @param x_mat Numeric matrix of conditioning points, one per row.
#' @param Q_vec Numeric vector with one Y value per row of `x_mat`.
#' @return Numeric vector of estimated conditional densities.
#' @export
fqx2_vec <- function(Y, X, x_mat, Q_vec, h) {
  mapply(function(i, Q) {
    fqx2(Y, X, x_mat[i, ], Q, h)
  }, seq_len(nrow(x_mat)), Q_vec)
}
