# Cross-validated selection of the number of quantiles K

For each candidate K, computes the out-of-fold mean squared error of the
local linear DARE prediction and returns the K with the smallest
average. Folds are drawn after `set.seed(1)`.

## Usage

``` r
cv_select_K_local(Y, X, K_grid, L, h)
```

## Arguments

- Y:

  Numeric response vector of length n.

- X:

  Numeric matrix or data frame of covariates (n x d).

- K_grid:

  Integer vector of candidate K values.

- L:

  Integer, number of folds.

- h:

  Bandwidth (scalar).

## Value

The selected K, with attributes `avg_MSE` (average MSE per candidate)
and `MSE_mat` (L x length(K_grid) matrix of fold MSEs).
