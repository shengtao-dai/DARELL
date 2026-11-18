# DARELL

Distributionally Adaptive Regression Estimator under Local Linear
Regression.

## Overview

Least squares regression can be very inefficient when the error
distribution is far from Gaussian, for example with heavy tails or
outliers. The **distributionally adaptive regression estimator (DARE)**
of [Dai (2025)](https://ssrn.com/abstract=5683602) estimates the
conditional mean $`E[Y \mid X = x_0]`$ as a weighted combination of
conditional quantiles at the levels $`k/(K+1)`$, $`k = 1, \dots, K`$.
The weights adapt to the unknown error distribution, and the number of
quantiles $`K`$ is chosen by cross-validation. Under regularity
conditions DARE attains the asymptotic variance of the infeasible
likelihood-based estimator while keeping the convergence rate of
conventional estimators.

`DARELL` implements the **local linear** version of DARE. Conditional
quantiles at the target point are estimated by kernel-weighted local
linear quantile regression, and the weights come from kernel estimates
of the conditional density at each quantile. The estimate comes with a
standard error and a 95% confidence interval.

Documentation: <https://shengtao-dai.github.io/DARELL/>

## Installation

``` r

# From r-universe (no compiler needed)
install.packages("DARELL", repos = c("https://shengtao-dai.r-universe.dev", "https://cloud.r-project.org"))

# Or the development version from GitHub
# install.packages("pak")
pak::pak("shengtao-dai/DARELL")
```

## Example

``` r

library(DARELL)

# Simulated data with heavy-tailed (t, 3 df) errors; the true m(0, 0) is 0.5
set.seed(123)
n <- 200
X <- cbind(rnorm(n), rnorm(n))
Y <- as.numeric(0.5 + X %*% c(1, -1) + rt(n, df = 3))

res <- local_eq_estimate(Y, X, x0 = c(0, 0), K_grid = c(1, 3, 5), L = 2)
unlist(res)
#>       m_hat_local_eq     cov_hat_local_eq m_hat_local_eq_lower
#>           0.49267987           0.29562481          -0.08674475
#> m_hat_local_eq_upper
#>           1.07210450
```

## Main functions

| Function | Purpose |
|----|----|
| [`local_eq_estimate()`](https://shengtao-dai.github.io/DARELL/reference/local_eq_estimate.md) | DARE estimate of E\[Y \| X = x0\], with standard error and 95% CI |
| [`cv_select_K_local()`](https://shengtao-dai.github.io/DARELL/reference/cv_select_K_local.md) | Cross-validated choice of the number of quantiles K |
| [`local_quantile()`](https://shengtao-dai.github.io/DARELL/reference/local_quantile.md), [`local_quantile_vec()`](https://shengtao-dai.github.io/DARELL/reference/local_quantile_vec.md) | Local linear quantile regression at one or several points |
| [`fqx2()`](https://shengtao-dai.github.io/DARELL/reference/fqx2.md), [`fqx()`](https://shengtao-dai.github.io/DARELL/reference/fqx.md), [`fqx2_vec()`](https://shengtao-dai.github.io/DARELL/reference/fqx2_vec.md), [`fx_den()`](https://shengtao-dai.github.io/DARELL/reference/fx_den.md) | Kernel estimates of conditional and covariate densities |
| [`inf_cq_local()`](https://shengtao-dai.github.io/DARELL/reference/inf_cq_local.md) | Variance factor of an equally weighted quantile average |

## Related packages

- [DARELL](https://github.com/shengtao-dai/DARELL): DARE with local
  linear regression (this package)
- [DARERF](https://github.com/shengtao-dai/DARERF): DARE with random
  forests
- [DARESR](https://github.com/shengtao-dai/DARESR): DARE with series
  (B-spline) expansion

## Citation

If you use `DARELL` in your work, please cite:

> Dai, S. (2025). *Distributionally Adaptive Regression Estimators*.
> Working paper. <https://ssrn.com/abstract=5683602>

``` bibtex
@unpublished{dai2025dare,
  title  = {Distributionally Adaptive Regression Estimators},
  author = {Dai, Shengtao},
  year   = {2025},
  note   = {Working paper, SSRN 5683602},
  url    = {https://ssrn.com/abstract=5683602}
}
```

You can also run `citation("DARELL")` in R.

## License

MIT © Shengtao Dai
