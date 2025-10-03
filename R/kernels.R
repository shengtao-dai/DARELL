# Kernels

#' Epanechnikov kernel

kernel_ep <- function(u) 0.75 * (1 - u^2) * (abs(u) <= 1)

#' Gaussian kernel

kernel_gau <- function(u) exp(-0.5 * u^2) / sqrt(2 * pi)
