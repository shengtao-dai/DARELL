test_that("basic run on small data", {
  set.seed(123)
  n <- 50
  X <- cbind(rnorm(n), rnorm(n))
  beta <- c(1, -1)
  Y <- 0.5 + X %*% beta + rnorm(n)
  x0 <- c(0, 0)
  K_grid <- c(1, 4, 9)

  expect_error(local_eq_estimate(Y, X, x0, c_band = 1.0, K_grid = K_grid, L = 2), NA)
})
