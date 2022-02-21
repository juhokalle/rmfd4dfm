
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
standard with packages hosted by GitHub.

``` r
# install.packages("devtools")
devtools::install_github("juhokalle/rmfd4dfm")
```

## The Model

Here I will provide a brief and hopefully accessible introduction to the
modelling approach presented in the paper. The aim is to identify and
estimate impulse-response functions (IRFs) using dynamic factor models
(DFMs). The model is given as

## Data

## Replication of the monetary policy example

-   The idea of the empirical example
-   The functions used in the exercise
-   Obtaining the results
