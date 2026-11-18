# Distributionally adaptive regression estimator (local linear)

Estimates the conditional mean \\m(x_0) = E\[Y \mid X = x_0\]\\ by a
weighted combination of local linear quantile estimates at the levels
\\\tau_k = k / (K + 1)\\, \\k = 1, \dots, K\\. The weights are built
from the estimated conditional densities at those quantiles, and K is
chosen from `K_grid` by `L`-fold cross-validation with
[`cv_select_K_local()`](https://shengtao-dai.github.io/DARELL/reference/cv_select_K_local.md).

## Usage

``` r
local_eq_estimate(Y, X, x0, c_band = 1, K_grid, L = 5)
```

## Arguments

- Y:

  Numeric response vector of length n.

- X:

  Numeric matrix or data frame of covariates (n x d).

- x0:

  Numeric vector of length d, the target point.

- c_band:

  Numeric scalar, the bandwidth constant.

- K_grid:

  Integer vector of candidate numbers of quantiles K.

- L:

  Integer, number of cross-validation folds.

## Value

A list with

- `m_hat_local_eq`: the estimate of \\m(x_0)\\;

- `cov_hat_local_eq`: its estimated standard error;

- `m_hat_local_eq_lower`, `m_hat_local_eq_upper`: the 95% confidence
  interval, `m_hat_local_eq` -/+ 1.96 \* `cov_hat_local_eq`.

## Details

The bandwidth is \\h = c\_{band} \\ n^{-1/(d+1)}\\. Cross-validation
folds are drawn after `set.seed(1)`, so results are reproducible, but
the global random number generator is reset.

## References

Dai, S. (2025). Distributionally Adaptive Regression Estimators. Working
paper. <https://ssrn.com/abstract=5683602>

## Examples

``` r
set.seed(123)
n <- 100
X <- cbind(rnorm(n), rnorm(n))
Y <- 0.5 + X %*% c(1, -1) + rt(n, df = 3)  # heavy-tailed errors
res <- local_eq_estimate(Y, X, x0 = c(0, 0), K_grid = c(1, 3, 5), L = 2)
res$m_hat_local_eq  # true value is 0.5
#> [1] 0.6438242
```
