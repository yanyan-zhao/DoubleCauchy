
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DoubleCauchy

<!-- badges: start -->

<!-- badges: end -->

The goal of DoubleCauchy is to provide a stable and adaptive test for
simultaneously testing regression coefficients of generalized linear
models of high dimensional data. The flexible framework proposed can be
directly applied to leverage side information for polygenic signal
detection.

## Installation

You can install the released version of DoubleCauchy with:

``` r
devtools::install_github("yanyan-zhao/DoubleCauchy")
```

## Example

``` r
library(DoubleCauchy)
## basic example code
```

Example 1: “AdapSide” function to leverage side information.

``` r
 J = 100
 n = 50
 beta=c(rep(2,3),rep(0,J-3))
 data=example_data_generate(n, J, beta)
 Y=data$Y
 G=data$G
 lam=0.1
 weights=c(J:1)
 pow.param=c(1:5)
 AdapSide(Y, G, cor.est="pdsoft", lam, weights, pow.param)
#>          p.1          p.2          p.3          p.4          p.5     p.cauchy 
#> 3.154025e-03 3.133633e-04 5.999152e-05 2.411928e-05 1.863413e-05 4.336407e-05
```

Example 2: A stable and adaptive test based on repeated sample splitting

``` r
 J = 100
 n = 50
 beta=c(rep(2,3),rep(0,J-3))
 data=example_data_generate(n, J, beta)
 Y=data$Y
 G=data$G
 lam=0.1
 weights=c(J:1)
 pow.param=c(1:5)
 n1=20
 m=10
 DoubleCauchy(n1, m, Y, G, varselec.method="DCSIS", J2=10, cor.est="pdsoft", lam=lam, pow.param=pow.param)
#>          p.1          p.2          p.3          p.4          p.5     p.cauchy 
#> 1.764015e-04 5.792850e-04 1.181325e-04 3.664634e-04 9.273655e-05 1.702288e-04
```

## Reference

Zhao and Sun (2020). A stable and adaptive polygenic signal detection
method based on repeated sample splitting.
