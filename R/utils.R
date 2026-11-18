#' Variance factor of an equally weighted quantile average
#'
#' For quantile levels \eqn{\tau_j = j / (K + 1)}, \eqn{j = 1, \dots, K},
#' computes \eqn{K^{-2} \sum_{j,k} (\min(\tau_j, \tau_k) - \tau_j \tau_k)
#' / (f_j f_k)}, the asymptotic variance factor of the equally weighted average
#' of K quantile estimators.
#'
#' @param f Numeric vector of length K, the (conditional) densities evaluated
#'   at the K quantiles.
#' @return Numeric scalar.
#' @export
inf_cq_local <- function(f) {
  K <- length(f)
  tau <- (1:K) / (K + 1)
  tau_j <- matrix(rep(tau, each = K), nrow = K)
  tau_k <- matrix(rep(tau, K), nrow = K)
  f_j <- matrix(rep(f, each = K), nrow = K)
  f_k <- matrix(rep(f, K), nrow = K)
  numerator <- pmin(tau_j, tau_k) - tau_j * tau_k
  denominator <- f_j * f_k
  sum(numerator / denominator) / K^2
}
