# Kernels

#' Epanechnikov kernel
#' @keywords internal
kernel_ep <- function(u) 0.75 * (1 - u^2) * (abs(u) <= 1)

#' Gaussian kernel (standard normal pdf)
#' @keywords internal
kernel_gau <- function(u) exp(-0.5 * u^2) / sqrt(2 * pi)
