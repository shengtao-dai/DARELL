# Kernel density estimate of the covariate density

Product Gaussian kernel estimate of \\f_X(x)\\.

## Usage

``` r
fx_den(X, x, h)
```

## Arguments

- X:

  Numeric matrix or data frame of covariates (n x d). A vector is
  treated as a single covariate.

- x:

  Numeric vector of length d, the evaluation point.

- h:

  Bandwidth (scalar).

## Value

Numeric scalar, the estimated density at `x`.
