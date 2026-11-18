# Joint kernel density of (Y, X) at (Q, x)

Joint kernel density of (Y, X) at (Q, x)

## Usage

``` r
fq(Y, X, x, Q, h)
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

Numeric scalar.
