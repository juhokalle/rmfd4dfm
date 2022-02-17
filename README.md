
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Estimation of Impulse-Response Functions with Dynamic Factor Models: A New Parametrization

<!-- badges: start -->
<!-- badges: end -->

This repository contains the associated R-package and data to the paper
[Estimation of Impulse-Response Functions with Dynamic Factor Models: A
New Parametrization](https://arxiv.org/pdf/2202.00310) authored by Juho
Koistinen and [Bernd
Funovits](https://sites.google.com/site/berndfunovits/) (2022).
Specifically, the functions within the package allow the user to
estimate the structural dynamic factor model, where with the common
component parametrized as the right matrix fraction description in
echelon form (RMFD-E), which is introduced in [Section
2.2](https://arxiv.org/pdf/2202.00310.pdf#subsection.2.2) of the
associated paper. This document also shows how to replicate the
empirical exercise of [Section
5](https://arxiv.org/pdf/2202.00310.pdf#section.5).

The package builds on the R-packages **rationalmatrices** and **RLDM**,
authored by Bernd Funovits and Wolfgang Scherrer. Since these packages
might change, the parts which are necessary for the analysis in the
associated article are extracted to R files **\~/rmfd4dfm/R/zz_ratmat**
and **\~/rmfd4dfm/R/zz_rldm**. Importantly, please reference the
**rationalmatrices** and the **RLDM** packages should you use their
functionalities and not this package.

## Installation

The package can be installed using the following command as is the
standard with packages hosted by GitHub

``` r
# install.packages("devtools")
devtools::install_github("juhokalle/rmfd4dfm")
```

## The argument in the paper

Here I will provide a brief and informal introduction to the argument
made in the paper as to why the methodology introduced in the paper

## Data

## Replication of the monetary policy example

-   The idea of the empirical example
-   The functions used in the exercise
-   Obtaining the results

``` r
library(rmfd4dfm)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/master/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />
