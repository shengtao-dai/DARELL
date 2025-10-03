
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
