# Cross validation to select K

#' Cross-validated selection of the number of quantiles K
#'
#' For each candidate K, computes the out-of-fold mean squared error of the
#' local linear DARE prediction and returns the K with the smallest average.
#' Folds are drawn after `set.seed(1)`.
#'
#' @param Y Numeric response vector of length n.
#' @param X Numeric matrix or data frame of covariates (n x d).
#' @param K_grid Integer vector of candidate K values.
#' @param L Integer, number of folds.
#' @param h Bandwidth (scalar).
#' @return The selected K, with attributes `avg_MSE` (average MSE per
#'   candidate) and `MSE_mat` (L x length(K_grid) matrix of fold MSEs).
#' @export
cv_select_K_local <- function(Y, X, K_grid, L, h) {
  n <- length(Y)
  if (is.vector(X)) X <- matrix(X, ncol = 1)
  set.seed(1L)
  folds <- sample(rep(seq_len(L), length.out = n))
  MSE_mat <- matrix(NA_real_, nrow = L, ncol = length(K_grid))

  for (l in seq_len(L)) {
    S_l <- which(folds == l)
    S_minus_l <- which(folds != l)
    X_train <- X[S_minus_l, , drop = FALSE]
    Y_train <- Y[S_minus_l]
    X_test  <- X[S_l, , drop = FALSE]
    Y_test  <- Y[S_l]

    for (k_idx in seq_along(K_grid)) {
      K <- K_grid[k_idx]
      tau_grid <- seq_len(K) / (K + 1)
      m_hat <- numeric(length(S_l))

      for (i in seq_along(S_l)) {
        x0 <- X_test[i, , drop = TRUE]
        Q_hat <- sapply(tau_grid, function(t) local_quantile(Y_train, X_train, x0, t, h))
        fqx_hat <- sapply(Q_hat, function(Q) fqx2(Y_train, X_train, x0, Q, h))
        fqx_hat_adj <- base::c(0, fqx_hat, 0)
        num_w <- fqx_hat * (-diff(diff(fqx_hat_adj)))
        deno_w <- sum(diff(fqx_hat_adj)^2)
        w_eq <- 0.5 * (num_w + rev(num_w)) / (deno_w + 1e-320)
        m_hat[i] <- sum(Q_hat * w_eq)
      }

      MSE_mat[l, k_idx] <- mean((Y_test - m_hat)^2)
    }
  }

  avg_MSE <- colMeans(MSE_mat)
  K_opt <- K_grid[which.min(avg_MSE)]
  attr(K_opt, "avg_MSE") <- avg_MSE
  attr(K_opt, "MSE_mat") <- MSE_mat
  K_opt
}
