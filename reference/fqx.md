# Conditional density of Y given X via a kernel ratio

Estimates \\f\_{Y\|X}(Q \mid x)\\ as the ratio of a joint kernel density
(Epanechnikov kernel in Y, product Gaussian kernel in X) to the
covariate kernel density.

## Usage

``` r
fqx(Y, X, x, Q, h)
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

## See also

[`fqx2()`](https://shengtao-dai.github.io/DARELL/reference/fqx2.md) for
the indicator-based version used by the estimator.
