# Variance factor of an equally weighted quantile average

For quantile levels \\\tau_j = j / (K + 1)\\, \\j = 1, \dots, K\\,
computes \\K^{-2} \sum\_{j,k} (\min(\tau_j, \tau_k) - \tau_j \tau_k) /
(f_j f_k)\\, the asymptotic variance factor of the equally weighted
average of K quantile estimators.

## Usage

``` r
inf_cq_local(f)
```

## Arguments

- f:

  Numeric vector of length K, the (conditional) densities evaluated at
  the K quantiles.

## Value

Numeric scalar.
