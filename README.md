# DARELL

An R package for a distributionally adaptive regression estimator under Local Linear Regression.

## Installation

```r
# Install devtools if needed
# install.packages("devtools")

devtools::install_local("DARELL", upgrade = "never")
# or build first:
# devtools::build("DARELL"); devtools::install("DARELL_0.1.0.tar.gz")
```

## Usage

```r
library(DARELL)

# Suppose Y (numeric), X (matrix/data.frame), and x0 (numeric vector)
res <- local_eq_estimate(Y, X, x0, c_band = 1.0, K_grid = c(1,4,9), L = 5)
str(res)
```



