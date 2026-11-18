# Conditional density of Y given X via an indicator approximation

Estimates \\f\_{Y\|X}(Q \mid x)\\ by the kernel-weighted share of
observations with \\\|Y - Q\| \< h\\, divided by \\2h\\. This is the
density estimate used inside
[`local_eq_estimate()`](https://shengtao-dai.github.io/DARELL/reference/local_eq_estimate.md).

## Usage

``` r
fqx2(Y, X, x, Q, h)
```

## Arguments

- Y:

  Numeric response vector of length n.

- X:

  Numeric matrix or data frame of covariates (n x d). A vector is
  treated as a single covariate.

- x:

  Numeric vector of length d, the conditioning point.

- Q:

  Numeric scalar, the value of Y at which the density is evaluated.

- h:

  Bandwidth (scalar).

## Value

Numeric scalar, the estimated conditional density.
