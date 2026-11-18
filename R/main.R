# Main estimator

#' Distributionally adaptive regression estimator (local linear)
#'
#' Estimates the conditional mean \eqn{m(x_0) = E[Y \mid X = x_0]} by a
#' weighted combination of local linear quantile estimates at the levels
#' \eqn{\tau_k = k / (K + 1)}, \eqn{k = 1, \dots, K}. The weights are built
#' from the estimated conditional densities at those quantiles, and K is chosen
#' from `K_grid` by `L`-fold cross-validation with [cv_select_K_local()].
#'
#' The bandwidth is \eqn{h = c_{band} \, n^{-1/(d+1)}}. Cross-validation folds
#' are drawn after `set.seed(1)`, so results are reproducible, but the global
#' random number generator is reset.
#'
#' @param Y Numeric response vector of length n.
#' @param X Numeric matrix or data frame of covariates (n x d).
#' @param x0 Numeric vector of length d, the target point.
#' @param c_band Numeric scalar, the bandwidth constant.
#' @param K_grid Integer vector of candidate numbers of quantiles K.
#' @param L Integer, number of cross-validation folds.
#' @return A list with
#'   * `m_hat_local_eq`: the estimate of \eqn{m(x_0)};
#'   * `cov_hat_local_eq`: its estimated standard error;
#'   * `m_hat_local_eq_lower`, `m_hat_local_eq_upper`: the 95% confidence
#'     interval, `m_hat_local_eq` -/+ 1.96 * `cov_hat_local_eq`.
#' @references Dai, S. (2025). Distributionally Adaptive Regression
#'   Estimators. Working paper. \url{https://ssrn.com/abstract=5683602}
#' @examples
#' set.seed(123)
#' n <- 100
#' X <- cbind(rnorm(n), rnorm(n))
#' Y <- 0.5 + X %*% c(1, -1) + rt(n, df = 3)  # heavy-tailed errors
#' res <- local_eq_estimate(Y, X, x0 = c(0, 0), K_grid = c(1, 3, 5), L = 2)
#' res$m_hat_local_eq  # true value is 0.5
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

