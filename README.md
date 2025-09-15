# DARELL

An R package for a distributionally adaptive regression estimator under Local Linear Regression.

## Installation

```r
# Option 1: pak (recommended)
install.packages("pak")
pak::pak("shengtao-dai/DARELL")

# Option 2: remotes
install.packages("remotes")
remotes::install_github("shengtao-dai/DARELL")

# Option 3: devtools
install.packages("devtools")
devtools::install_github("shengtao-dai/DARELL")

```

## Usage

```r
library(DARELL)

# Suppose Y (numeric), X (matrix/data.frame), and x0 (numeric vector)
res <- local_eq_estimate(Y, X, x0, c_band = 1.0, K_grid = c(1,4,9), L = 5)
str(res)
```



