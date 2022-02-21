Estimation of Impulse-Response Functions with Dynamic Factor Models: A
New Parametrization
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

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

## The Model and Argument

Here I will provide a brief and hopefully accessible introduction to the
modelling approach presented in the paper, while the details are
available from the paper. The aim is to identify and estimate
impulse-response functions (IRFs) using dynamic factor models (DFMs),
with the emphasis placed on macroeconomic applications and and analysis.
The model of interest is given as

![
x\_{t} =d\_{0}z\_{t}^{\*}+d\_{1}z\_{t-1}^{\*}+\\cdots+d\_{s}z\_{t-s}^{\*}+\\xi\_{t}\\\\
z\_{t}^{\*} =c\_{1}z\_{t-1}^{\*}+\\cdots+c\_{p}z\_{t-p}^{\*}+\\varepsilon\_{t}, \\\\
\\varepsilon\_{t}=Hu\_{t},
](https://latex.codecogs.com/png.latex?%0Ax_%7Bt%7D%20%3Dd_%7B0%7Dz_%7Bt%7D%5E%7B%2A%7D%2Bd_%7B1%7Dz_%7Bt-1%7D%5E%7B%2A%7D%2B%5Ccdots%2Bd_%7Bs%7Dz_%7Bt-s%7D%5E%7B%2A%7D%2B%5Cxi_%7Bt%7D%5C%5C%0Az_%7Bt%7D%5E%7B%2A%7D%20%3Dc_%7B1%7Dz_%7Bt-1%7D%5E%7B%2A%7D%2B%5Ccdots%2Bc_%7Bp%7Dz_%7Bt-p%7D%5E%7B%2A%7D%2B%5Cvarepsilon_%7Bt%7D%2C%20%5C%5C%0A%5Cvarepsilon_%7Bt%7D%3DHu_%7Bt%7D%2C%0A "
x_{t} =d_{0}z_{t}^{*}+d_{1}z_{t-1}^{*}+\cdots+d_{s}z_{t-s}^{*}+\xi_{t}\\
z_{t}^{*} =c_{1}z_{t-1}^{*}+\cdots+c_{p}z_{t-p}^{*}+\varepsilon_{t}, \\
\varepsilon_{t}=Hu_{t},
")

where

-   ![x_t](https://latex.codecogs.com/png.latex?x_t "x_t") is an
    ![n](https://latex.codecogs.com/png.latex?n "n")-dimensional
    observed time series, with usually
    ![n>100](https://latex.codecogs.com/png.latex?n%3E100 "n>100").

-   ![z_t^\*](https://latex.codecogs.com/png.latex?z_t%5E%2A "z_t^*") is
    a ![q](https://latex.codecogs.com/png.latex?q "q")-dimensional
    dynamic factor process, with usually
    ![q\<10](https://latex.codecogs.com/png.latex?q%3C10 "q<10").

-   ![xi_t](https://latex.codecogs.com/png.latex?xi_t "xi_t") is an
    ![n](https://latex.codecogs.com/png.latex?n "n")-dimensional
    idiosyncratic term, we assume
    ![\\mathbb E(\\xi_t \\xi_t')=\\sigma\_\\xi^2I_n](https://latex.codecogs.com/png.latex?%5Cmathbb%20E%28%5Cxi_t%20%5Cxi_t%27%29%3D%5Csigma_%5Cxi%5E2I_n "\mathbb E(\xi_t \xi_t')=\sigma_\xi^2I_n").

-   ![\\varepsilon_t](https://latex.codecogs.com/png.latex?%5Cvarepsilon_t "\varepsilon_t")
    is the innovation to
    ![z_t^\*](https://latex.codecogs.com/png.latex?z_t%5E%2A "z_t^*"),
    with
    ![\\mathbb E(\\varepsilon_t \\varepsilon_t')=\\Sigma\_\\varepsilon](https://latex.codecogs.com/png.latex?%5Cmathbb%20E%28%5Cvarepsilon_t%20%5Cvarepsilon_t%27%29%3D%5CSigma_%5Cvarepsilon "\mathbb E(\varepsilon_t \varepsilon_t')=\Sigma_\varepsilon").

-   ![d_i](https://latex.codecogs.com/png.latex?d_i "d_i") and
    ![c_j](https://latex.codecogs.com/png.latex?c_j "c_j"),
    ![i=0,\\ldots,s](https://latex.codecogs.com/png.latex?i%3D0%2C%5Cldots%2Cs "i=0,\ldots,s")
    ![j=1,\\ldots,p](https://latex.codecogs.com/png.latex?j%3D1%2C%5Cldots%2Cp "j=1,\ldots,p"),
    are
    ![n \\times q](https://latex.codecogs.com/png.latex?n%20%5Ctimes%20q "n \times q")
    and
    ![q \\times q](https://latex.codecogs.com/png.latex?q%20%5Ctimes%20q "q \times q")
    parameter matrices, values of which are at the center of interest as
    the IRFs are constructed from these matrices.

-   ![u_t](https://latex.codecogs.com/png.latex?u_t "u_t") is the
    ![q](https://latex.codecogs.com/png.latex?q "q")-dimensional
    structural shock process, with
    ![\\mathbb E (u_t u_t') = I_q](https://latex.codecogs.com/png.latex?%5Cmathbb%20E%20%28u_t%20u_t%27%29%20%3D%20I_q "\mathbb E (u_t u_t') = I_q").

-   ![H](https://latex.codecogs.com/png.latex?H "H") is
    ![(q\\times q)](https://latex.codecogs.com/png.latex?%28q%5Ctimes%20q%29 "(q\times q)")-dimensional
    structural shock impact matrix, which we identify as
    ![H=chol(\\Sigma\_\\varepsilon)](https://latex.codecogs.com/png.latex?H%3Dchol%28%5CSigma_%5Cvarepsilon%29 "H=chol(\Sigma_\varepsilon)"),
    where
    ![chol(\\Sigma\_\\varepsilon)](https://latex.codecogs.com/png.latex?chol%28%5CSigma_%5Cvarepsilon%29 "chol(\Sigma_\varepsilon)")
    is the lower triangular Cholesky factor of
    ![\\Sigma\_\\varepsilon](https://latex.codecogs.com/png.latex?%5CSigma_%5Cvarepsilon "\Sigma_\varepsilon").
    Note that any other identification method for uncovering
    ![H](https://latex.codecogs.com/png.latex?H "H") is valid as in the
    structural vector autoregression (SVAR) analysis.

We can write the dynamic factor process compactly as
![z_t^\*=c(L)^{-1}\\varepsilon_t](https://latex.codecogs.com/png.latex?z_t%5E%2A%3Dc%28L%29%5E%7B-1%7D%5Cvarepsilon_t "z_t^*=c(L)^{-1}\varepsilon_t")
and substitute this and the equation for the innovations
![\\varepsilon_t](https://latex.codecogs.com/png.latex?%5Cvarepsilon_t "\varepsilon_t")
into the equation for the observations to get
![x_t=d(L)c(L)^{-1}H u_t + \\xi_t](https://latex.codecogs.com/png.latex?x_t%3Dd%28L%29c%28L%29%5E%7B-1%7DH%20u_t%20%2B%20%5Cxi_t "x_t=d(L)c(L)^{-1}H u_t + \xi_t").
If we assume that the idiosyncratic component accounts for measurements
errors or sectoral dynamics that pertain to a small number of variables,
we can use the structural IRF
![k(L)H](https://latex.codecogs.com/png.latex?k%28L%29H "k(L)H") to
study the shock propagation from
![q](https://latex.codecogs.com/png.latex?q "q") structural shocks to
![n](https://latex.codecogs.com/png.latex?n "n") observed macroeconomic
variables, with
![k(L)=d(L)c(L)^{-1}](https://latex.codecogs.com/png.latex?k%28L%29%3Dd%28L%29c%28L%29%5E%7B-1%7D "k(L)=d(L)c(L)^{-1}").
The IRF ![k(L)](https://latex.codecogs.com/png.latex?k%28L%29 "k(L)") is
not identified without further restrictions since one can always
post-multiply
![d(L)](https://latex.codecogs.com/png.latex?d%28L%29 "d(L)") and
![c(L)](https://latex.codecogs.com/png.latex?c%28L%29 "c(L)") by some
![q\\times q](https://latex.codecogs.com/png.latex?q%5Ctimes%20q "q\times q")
polynomial matrix
![m(L)](https://latex.codecogs.com/png.latex?m%28L%29 "m(L)") such that
it “cancels out”
![k(L)=d(L)c(L)^{-1}=\[d(L)m(L)\]\[c(L)m(L)\]^{-1}=\\bar d(L) \\bar c(L)^{-1}](https://latex.codecogs.com/png.latex?k%28L%29%3Dd%28L%29c%28L%29%5E%7B-1%7D%3D%5Bd%28L%29m%28L%29%5D%5Bc%28L%29m%28L%29%5D%5E%7B-1%7D%3D%5Cbar%20d%28L%29%20%5Cbar%20c%28L%29%5E%7B-1%7D "k(L)=d(L)c(L)^{-1}=[d(L)m(L)][c(L)m(L)]^{-1}=\bar d(L) \bar c(L)^{-1}").
The problem is that the researcher cannot distinguish between
![d(L)c(L)^{-1}](https://latex.codecogs.com/png.latex?d%28L%29c%28L%29%5E%7B-1%7D "d(L)c(L)^{-1}")
and
![\\bar d(L) \\bar c(L)^{-1}](https://latex.codecogs.com/png.latex?%5Cbar%20d%28L%29%20%5Cbar%20c%28L%29%5E%7B-1%7D "\bar d(L) \bar c(L)^{-1}")
from the first and second moments of the data, and therefore conclusions
drawn from the structural IRF
![k(z)H](https://latex.codecogs.com/png.latex?k%28z%29H "k(z)H") are
meaningless.

Our insight is that
![k(L)](https://latex.codecogs.com/png.latex?k%28L%29 "k(L)") can be
identified, that is, the set of matrices
![m(L)](https://latex.codecogs.com/png.latex?m%28L%29 "m(L)") can be
narrowed down to ![I_q](https://latex.codecogs.com/png.latex?I_q "I_q"),
using the identification restrictions that are standard in the
literature dealing with the identification of the vector autoregressive
moving average (VARMA) models. In this model class, the IRF is given as
![k(L)=a(L)^{-1}b(L)](https://latex.codecogs.com/png.latex?k%28L%29%3Da%28L%29%5E%7B-1%7Db%28L%29 "k(L)=a(L)^{-1}b(L)")
for some AR and MA lag polynomials
![a(L)](https://latex.codecogs.com/png.latex?a%28L%29 "a(L)") and
![b(L)](https://latex.codecogs.com/png.latex?b%28L%29 "b(L)") of
suitable dimensions, and here
![a(L)](https://latex.codecogs.com/png.latex?a%28L%29 "a(L)") and
![b(L)](https://latex.codecogs.com/png.latex?b%28L%29 "b(L)") can be
pre-multiplied by some lag polynomial
![m(L)](https://latex.codecogs.com/png.latex?m%28L%29 "m(L)") to obtain
an observationally equivalent IRF:
![k(L)=a(L)^{-1}b(L)=\[m(L)a(L)\]^{-1}\[m(L)b(L)\]=\\bar a(L)^{-1}\\bar b(L)](https://latex.codecogs.com/png.latex?k%28L%29%3Da%28L%29%5E%7B-1%7Db%28L%29%3D%5Bm%28L%29a%28L%29%5D%5E%7B-1%7D%5Bm%28L%29b%28L%29%5D%3D%5Cbar%20a%28L%29%5E%7B-1%7D%5Cbar%20b%28L%29 "k(L)=a(L)^{-1}b(L)=[m(L)a(L)]^{-1}[m(L)b(L)]=\bar a(L)^{-1}\bar b(L)").
One popular identification approach for the VARMA model is to use the
canonical echelon form parametrization, which involves zero and unity
restrictions in
![a(L)](https://latex.codecogs.com/png.latex?a%28L%29 "a(L)") and
![b(L)](https://latex.codecogs.com/png.latex?b%28L%29 "b(L)"). By noting
that the identification restrictions for the IRF in right matrix
fraction description (RMFD, for short) model,
i.e. ![k(L)=d(L)c(L)^{-1}](https://latex.codecogs.com/png.latex?k%28L%29%3Dd%28L%29c%28L%29%5E%7B-1%7D "k(L)=d(L)c(L)^{-1}"),
are equivalent to those for the VARMA model after transposing, the
derivation is quite straightforward using the existing results for the
VARMA model (for details, see [Appendix
A](https://arxiv.org/pdf/2202.00310.pdf#section*.2)).

The structure of the RMFD model in echelon form is quite complicated,
and while the functions of this package `rmfd4dfm` code these
restrictions without the need for the user to trouble herself with the
implementation, one might still wonder what is the point of going
through through the painstaking process of deriving and coding the
restrictions in the first place, as there are alternatives in the
literature? To make the point, let us introduce a popular alternative
used in the literature and compare it with our approach. One way to deal
with the identification and estimation of the DFM presented above is to
stack the dynamic factors into a
![r=q(s+1)](https://latex.codecogs.com/png.latex?r%3Dq%28s%2B1%29 "r=q(s+1)")-dimensional
vector
![z_t = \\left(z_t^{\*'}, z\_{t-1}^{\*'}, \\cdots, z\_{t-s}^{\*'} \\right)'](https://latex.codecogs.com/png.latex?z_t%20%3D%20%5Cleft%28z_t%5E%7B%2A%27%7D%2C%20z_%7Bt-1%7D%5E%7B%2A%27%7D%2C%20%5Ccdots%2C%20z_%7Bt-s%7D%5E%7B%2A%27%7D%20%5Cright%29%27 "z_t = \left(z_t^{*'}, z_{t-1}^{*'}, \cdots, z_{t-s}^{*'} \right)'"),
and write the equation for
![x_t](https://latex.codecogs.com/png.latex?x_t "x_t") as
![x_t=Dz_t + \\xi_t](https://latex.codecogs.com/png.latex?x_t%3DDz_t%20%2B%20%5Cxi_t "x_t=Dz_t + \xi_t"),
where
![D=\\left(d_0,\\cdots,d_s\\right)](https://latex.codecogs.com/png.latex?D%3D%5Cleft%28d_0%2C%5Ccdots%2Cd_s%5Cright%29 "D=\left(d_0,\cdots,d_s\right)"),
which makes the model amenable to the principal components (PC) methods.
In the second step, the static factor process is modelled as a VAR
process, A(L)z_t = B\_t, where
![A(L)](https://latex.codecogs.com/png.latex?A%28L%29 "A(L)") is an
![r \\times r](https://latex.codecogs.com/png.latex?r%20%5Ctimes%20r "r \times r")
lag polynomial, and ![B](https://latex.codecogs.com/png.latex?B "B") is
an
![r\\times q](https://latex.codecogs.com/png.latex?r%5Ctimes%20q "r\times q")
constant matrix. After the estimation of D and A(L), the IRF
![k(L)=DA(L)^{-1}B](https://latex.codecogs.com/png.latex?k%28L%29%3DDA%28L%29%5E%7B-1%7DB "k(L)=DA(L)^{-1}B")
can be constructed straightforwardly.

This alternative using PCs is straightforward from an estimation point
of view, but let us point two associated restrictive features. First,
note that whenever
![s>0](https://latex.codecogs.com/png.latex?s%3E0 "s>0"), the minimal
number of parameter restrictions needed is
![r=q(s+1)](https://latex.codecogs.com/png.latex?r%3Dq%28s%2B1%29 "r=q(s+1)"),
which is higher than that of needed in our parametrization,
i.e. ![q](https://latex.codecogs.com/png.latex?q "q"). Second, we note
that the reliable estimation of VAR on
![z_t](https://latex.codecogs.com/png.latex?z_t "z_t") can be difficult
if the lag order is misspecified, which is due to the singularity of the
covariance matrix of the innovation to
![z_t](https://latex.codecogs.com/png.latex?z_t "z_t"),
i.e. ![\\mathbb E (B\\varepsilon_t \\varepsilon_t'B')=B\\Sigma\_\\varepsilon B'](https://latex.codecogs.com/png.latex?%5Cmathbb%20E%20%28B%5Cvarepsilon_t%20%5Cvarepsilon_t%27B%27%29%3DB%5CSigma_%5Cvarepsilon%20B%27 "\mathbb E (B\varepsilon_t \varepsilon_t'B')=B\Sigma_\varepsilon B'")
has rank ![q](https://latex.codecogs.com/png.latex?q "q"), which is
smaller than
![r=(s+1)q](https://latex.codecogs.com/png.latex?r%3D%28s%2B1%29q "r=(s+1)q")
whenever ![s>0](https://latex.codecogs.com/png.latex?s%3E0 "s>0").

## Data

## Replication of the monetary policy example

-   The idea of the empirical example
-   The functions used in the exercise
-   Obtaining the results
