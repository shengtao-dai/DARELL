# Local linear quantile estimates at several points

Applies
[`local_quantile()`](https://shengtao-dai.github.io/DARELL/reference/local_quantile.md)
to each row of `X0`.

## Usage

``` r
local_quantile_vec(Y, X, X0, tau, h)
```

## Arguments

- Y:

  Numeric response vector of length n.

- X:

  Numeric matrix or data frame of covariates (n x d).

- X0:

  Numeric matrix of target points, one per row.

- tau:

  Quantile level in (0, 1).

- h:

  Bandwidth (scalar).

## Value

Numeric vector of estimated conditional quantiles.
