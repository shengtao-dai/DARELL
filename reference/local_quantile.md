# Local linear quantile estimator at a point

Minimizes the Gaussian-kernel-weighted check loss of a local linear fit
around `x` with [`stats::optim()`](https://rdrr.io/r/stats/optim.html)
and returns the intercept.

## Usage

``` r
local_quantile(Y, X, x, tau, h)
```

## Arguments

- Y:

  Numeric response vector of length n.

- X:

  Numeric matrix or data frame of covariates (n x d).

- x:

  Numeric vector of length d, the target point.

- tau:

  Quantile level in (0, 1).

- h:

  Bandwidth (scalar).

## Value

Numeric scalar, the estimated conditional `tau`-quantile at `x`.
