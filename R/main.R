# Main estimator

#' Local equivalent-weight estimator (precision-style CI)
#'
#' This version matches the user's current script: it returns
#' cov_hat_local_eq = 1 / SE and uses CI as m_hat ± 1.96 * (1/SE).
#' See also \code{local_eq_estimate_se()} for the conventional SE-based CI.
#'
#' @param Y numeric vector (response)
#' @param X numeric matrix/data.frame (n x d)
#' @param x0 numeric vector of length d (target point)
#' @param c_band numeric scalar, bandwidth constant
#' @param K_grid integer vector of candidate K
#' @param L integer, number of folds for CV
#' @return list with m_hat_local_eq, cov_hat_local_eq, m_hat_local_eq_lower, m_hat_local_eq_upper
#' @export
local_eq_estimate <- function(Y, X, x0, c_band = 1, K_grid, L = 5) {
  if (is.vector(X)) X <- matrix(X, ncol = length(x0))
  n <- nrow(X)
  d <- length(x0)

  h   <- c_band * n^(-1 / (d + 1))
  con <- 0.2821^d
  ridge <- 1e-320

  K_opt <- cv_select_K_local(Y, X, K_grid, L, h)
  tau_grid <- seq_len(K_opt) / (K_opt + 1)

  Q_hat_vec <- sapply(tau_grid, function(t) local_quantile(Y, X, x0, t, h))
  fqx_hat   <- sapply(Q_hat_vec, function(Q) fqx2(Y, X, x0, Q, h))

  fqx_hat_adj <- base::c(0, fqx_hat, 0)
  num_w <- fqx_hat * (-diff(diff(fqx_hat_adj)))
  deno_w <- sum(diff(fqx_hat_adj)^2) + ridge
  w_eq <- 0.5 * (num_w + rev(num_w)) / deno_w

  m_hat_local_eq <- sum(Q_hat_vec * w_eq)

  fx0_hat <- fx_den(X, x0, h)
  se_hat <- sqrt(n * h^d * (K_opt + 1) * deno_w * fx0_hat / con)

  cov_hat_local_eq <- 1 / se_hat
  m_hat_local_eq_upper <- m_hat_local_eq + 1.96 * cov_hat_local_eq
  m_hat_local_eq_lower <- m_hat_local_eq - 1.96 * cov_hat_local_eq

  list(
    m_hat_local_eq       = m_hat_local_eq,
    cov_hat_local_eq     = cov_hat_local_eq,
    m_hat_local_eq_lower = m_hat_local_eq_lower,
    m_hat_local_eq_upper = m_hat_local_eq_upper
  )
}

