# Vectorized conditional density estimates

Applies
[`fqx2()`](https://shengtao-dai.github.io/DARELL/reference/fqx2.md) row
by row, pairing the i-th row of `x_mat` with the i-th element of
`Q_vec`.

## Usage

``` r
fqx2_vec(Y, X, x_mat, Q_vec, h)
```

## Arguments

- Y:

  Numeric response vector of length n.

- X:

  Numeric matrix or data frame of covariates (n x d). A vector is
  treated as a single covariate.

- x_mat:

  Numeric matrix of conditioning points, one per row.

- Q_vec:

  Numeric vector with one Y value per row of `x_mat`.

- h:

  Bandwidth (scalar).

## Value

Numeric vector of estimated conditional densities.
