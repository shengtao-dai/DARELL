# Kernel density and conditional density utilities

#' Kernel density estimate of f_X(x)

fx_den <- function(X, x, h) {
  if (is.vector(X)) X <- matrix(X, ncol = 1)
  d <- ncol(X)
  x <- as.vector(x)
  kernel_products <- matrixStats::rowProds(kernel_gau(sweep(X, 2, x, "-") / h))
  mean(kernel_products) / (h^d)
}

#' Joint kernel for (Y,X) evaluated at (Q,x)

fq <- function(Y, X, x, Q, h) {
  if (is.vector(X)) X <- matrix(X, ncol = 1)
  d <- ncol(X)
  x <- as.vector(x)
  kernel_products <- matrixStats::rowProds(kernel_gau(sweep(X, 2, x, "-") / h))
  kernel_q <- kernel_ep((Y - Q) / h)
  joint_prod <- kernel_products * kernel_q
  mean(joint_prod) / (h^(d + 1))
}

#' Conditional density f_{Y|X}(Q|x) via kernel ratio

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

#' Conditional density via indicator approximation

fqx2 <- function(Y, X, x, Q, h) {
  if (is.vector(X)) X <- matrix(X, ncol = 1)
  x <- as.vector(x)
  kernel_products <- matrixStats::rowProds(kernel_gau(sweep(X, 2, x, "-") / h))
  indicator <- abs(Y - Q) < h
  num <- sum(kernel_products * indicator)
  den <- sum(kernel_products)
  num / (den * 2 * h)
}

#' Vectorized variants

fqx2_vec <- function(Y, X, x_mat, Q_vec, h) {
  mapply(function(i, Q) {
    fqx2(Y, X, x_mat[i, ], Q, h)
  }, seq_len(nrow(x_mat)), Q_vec)
}
