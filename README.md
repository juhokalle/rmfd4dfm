Estimation of Impulse-Response Functions with Dynamic Factor Models: A
New Parametrization
================

<!-- badges: start -->

[![R-CMD-check](https://github.com/juhokalle/rmfd4dfm/workflows/R-CMD-check/badge.svg)](https://github.com/juhokalle/rmfd4dfm/actions)
<!-- badges: end -->

This repository contains the associated R-package and data to the paper
[Estimation of Impulse-Response Functions with Dynamic Factor Models: A
New Parametrization](https://arxiv.org/pdf/2202.00310) authored by [Juho
Koistinen](https://juhokois.netlify.app/) and [Bernd
Funovits](https://sites.google.com/site/berndfunovits/) (2022).
Specifically, the functions within the package allow the user to
estimate a structural dynamic factor model (DFM), with the common
component parametrized as a right matrix fraction description in echelon
form (RMFD-E), which is introduced in [Section
2.2](https://arxiv.org/pdf/2202.00310.pdf#subsection.2.2) of the
associated paper. This document also shows how to replicate the
empirical exercise of [Section
4](https://arxiv.org/pdf/2202.00310.pdf#section.4).

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
with the emphasis placed on macroeconomic applications and analysis. The
model of interest is given as

![
 (1): \\quad x\_{t} =d\_{0}z\_{t}^{\*}+d\_{1}z\_{t-1}^{\*}+\\cdots+d\_{s}z\_{t-s}^{\*}+\\xi\_{t}
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%20%281%29%3A%20%5Cquad%20x_%7Bt%7D%20%3Dd_%7B0%7Dz_%7Bt%7D%5E%7B%2A%7D%2Bd_%7B1%7Dz_%7Bt-1%7D%5E%7B%2A%7D%2B%5Ccdots%2Bd_%7Bs%7Dz_%7Bt-s%7D%5E%7B%2A%7D%2B%5Cxi_%7Bt%7D%0A "
 (1): \quad x_{t} =d_{0}z_{t}^{*}+d_{1}z_{t-1}^{*}+\cdots+d_{s}z_{t-s}^{*}+\xi_{t}
")

![
(2): \\quad z\_{t}^{\*} =c\_{1}z\_{t-1}^{\*}+\\cdots+c\_{p}z\_{t-p}^{\*}+\\varepsilon\_{t}
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%282%29%3A%20%5Cquad%20z_%7Bt%7D%5E%7B%2A%7D%20%3Dc_%7B1%7Dz_%7Bt-1%7D%5E%7B%2A%7D%2B%5Ccdots%2Bc_%7Bp%7Dz_%7Bt-p%7D%5E%7B%2A%7D%2B%5Cvarepsilon_%7Bt%7D%0A "
(2): \quad z_{t}^{*} =c_{1}z_{t-1}^{*}+\cdots+c_{p}z_{t-p}^{*}+\varepsilon_{t}
")

![
(3): \\quad \\varepsilon\_{t}=Hu\_{t}
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%283%29%3A%20%5Cquad%20%5Cvarepsilon_%7Bt%7D%3DHu_%7Bt%7D%0A "
(3): \quad \varepsilon_{t}=Hu_{t}
")

where

-   ![x_t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_t "x_t")
    is an
    ![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n")-dimensional
    observed time series, with usually
    ![n\>100](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n%3E100 "n>100").
-   ![z_t^\*](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;z_t%5E%2A "z_t^*")
    is a
    ![q](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;q "q")-dimensional
    dynamic factor process, with usually
    ![q\<10](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;q%3C10 "q<10").
-   ![x_t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_t "x_t")
    is an
    ![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n")-dimensional
    idiosyncratic term, we assume
    ![\\mathbb E(\\xi_t \\xi_t')=\\sigma\_\\xi^2I_n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbb%20E%28%5Cxi_t%20%5Cxi_t%27%29%3D%5Csigma_%5Cxi%5E2I_n "\mathbb E(\xi_t \xi_t')=\sigma_\xi^2I_n").
-   ![\\varepsilon_t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cvarepsilon_t "\varepsilon_t")
    is the innovation to
    ![z_t^\*](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;z_t%5E%2A "z_t^*"),
    with
    ![\\mathbb E(\\varepsilon_t \\varepsilon_t')=\\Sigma\_\\varepsilon](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbb%20E%28%5Cvarepsilon_t%20%5Cvarepsilon_t%27%29%3D%5CSigma_%5Cvarepsilon "\mathbb E(\varepsilon_t \varepsilon_t')=\Sigma_\varepsilon").
-   ![\\{d_i\\}\_{i=0}^s](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5C%7Bd_i%5C%7D_%7Bi%3D0%7D%5Es "\{d_i\}_{i=0}^s")
    and
    ![\\{c_j\\}\_{i=1}^p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5C%7Bc_j%5C%7D_%7Bi%3D1%7D%5Ep "\{c_j\}_{i=1}^p")
    are
    ![n \\times q](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n%20%5Ctimes%20q "n \times q")
    and
    ![q \\times q](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;q%20%5Ctimes%20q "q \times q")
    parameter matrices, values of which are at the center of interest as
    the IRFs are constructed from these matrices.
-   ![u_t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;u_t "u_t")
    is the
    ![q](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;q "q")-dimensional
    structural shock process, with
    ![\\mathbb E (u_t u_t') = I_q](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbb%20E%20%28u_t%20u_t%27%29%20%3D%20I_q "\mathbb E (u_t u_t') = I_q").
-   ![H](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;H "H")
    is
    ![(q\\times q)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28q%5Ctimes%20q%29 "(q\times q)")-dimensional
    structural impact multiplier matrix, which we identify as
    ![H=chol(\\Sigma\_\\varepsilon)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;H%3Dchol%28%5CSigma_%5Cvarepsilon%29 "H=chol(\Sigma_\varepsilon)"),
    where
    ![chol(\\Sigma\_\\varepsilon)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;chol%28%5CSigma_%5Cvarepsilon%29 "chol(\Sigma_\varepsilon)")
    is the lower triangular Cholesky factor of
    ![\\Sigma\_\\varepsilon](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CSigma_%5Cvarepsilon "\Sigma_\varepsilon").
    Note that any other identification method for uncovering
    ![H](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;H "H")
    is valid as in the structural vector autoregression (SVAR) analysis.

We can write Eq. (2) compactly as
![z_t^\*=c(L)^{-1}\\varepsilon_t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;z_t%5E%2A%3Dc%28L%29%5E%7B-1%7D%5Cvarepsilon_t "z_t^*=c(L)^{-1}\varepsilon_t")
and substitute this and Eq. (3) into Eq. (1) to get

![
x_t=d(L)c(L)^{-1}H u_t + \\xi_t,
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0Ax_t%3Dd%28L%29c%28L%29%5E%7B-1%7DH%20u_t%20%2B%20%5Cxi_t%2C%0A "
x_t=d(L)c(L)^{-1}H u_t + \xi_t,
")

where

![
c(L)=I_q-c_1 L - \\cdots - c_p L^p
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0Ac%28L%29%3DI_q-c_1%20L%20-%20%5Ccdots%20-%20c_p%20L%5Ep%0A "
c(L)=I_q-c_1 L - \cdots - c_p L^p
")

and

![
d(L)=d_0 + d_1 L + \\cdots + d_sL^s.
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0Ad%28L%29%3Dd_0%20%2B%20d_1%20L%20%2B%20%5Ccdots%20%2B%20d_sL%5Es.%0A "
d(L)=d_0 + d_1 L + \cdots + d_sL^s.
")

If we assume that the idiosyncratic component accounts for measurements
errors or sectoral dynamics that pertain to a small number of variables,
we can use the structural IRF
![k(L)H](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;k%28L%29H "k(L)H")
to study the shock propagation from
![q](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;q "q")
structural shocks to
![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n")
observed macroeconomic variables, with
![k(L)=d(L)c(L)^{-1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;k%28L%29%3Dd%28L%29c%28L%29%5E%7B-1%7D "k(L)=d(L)c(L)^{-1}").
The IRF
![k(L)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;k%28L%29 "k(L)")
is not identified without further restrictions since one can always
post-multiply
![d(L)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;d%28L%29 "d(L)")
and
![c(L)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;c%28L%29 "c(L)")
by some
![q\\times q](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;q%5Ctimes%20q "q\times q")
polynomial matrix
![m(L)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;m%28L%29 "m(L)")
such that it “cancels out”:

![
k(L)=d(L)c(L)^{-1}=\[d(L)m(L)\]\[c(L)m(L)\]^{-1}=\\bar d(L) \\bar c(L)^{-1}.
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0Ak%28L%29%3Dd%28L%29c%28L%29%5E%7B-1%7D%3D%5Bd%28L%29m%28L%29%5D%5Bc%28L%29m%28L%29%5D%5E%7B-1%7D%3D%5Cbar%20d%28L%29%20%5Cbar%20c%28L%29%5E%7B-1%7D.%0A "
k(L)=d(L)c(L)^{-1}=[d(L)m(L)][c(L)m(L)]^{-1}=\bar d(L) \bar c(L)^{-1}.
")

The problem is that the researcher cannot distinguish between
![d(L)c(L)^{-1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;d%28L%29c%28L%29%5E%7B-1%7D "d(L)c(L)^{-1}")
and
![\\bar d(L) \\bar c(L)^{-1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbar%20d%28L%29%20%5Cbar%20c%28L%29%5E%7B-1%7D "\bar d(L) \bar c(L)^{-1}")
from the first and second moments of the data, and therefore attempts at
drawing conclusions from the structural IRF
![k(z)H](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;k%28z%29H "k(z)H")
are futile without further assumptions on
![c(L)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;c%28L%29 "c(L)")
and
![d(L)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;d%28L%29 "d(L)").

Our insight is that
![k(L)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;k%28L%29 "k(L)")
can be identified, that is, the set of matrices
![m(L)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;m%28L%29 "m(L)")
can be narrowed down to
![I_q](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;I_q "I_q"),
using the identification restrictions that are standard in the
literature dealing with the identification of the vector autoregressive
moving average (VARMA) models (for a recent summary, see [Deistler and
Scherrer,
2019](https://www2.cirano.qc.ca/~dufourj/Web_Site/Vinod_Rao_2019_Bk_Elsevier_ConceptualEtxUsingR.pdf#page=162)).
In this model class, the IRF is given as
![k(L)=a(L)^{-1}b(L)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;k%28L%29%3Da%28L%29%5E%7B-1%7Db%28L%29 "k(L)=a(L)^{-1}b(L)")
for some AR and MA lag polynomials
![a(L)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;a%28L%29 "a(L)")
and
![b(L)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;b%28L%29 "b(L)")
of suitable dimensions, and here
![a(L)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;a%28L%29 "a(L)")
and
![b(L)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;b%28L%29 "b(L)")
can be pre-multiplied by some lag polynomial
![m(L)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;m%28L%29 "m(L)")
to obtain an observationally equivalent IRF:

![
k(L)=a(L)^{-1}b(L)=\[m(L)a(L)\]^{-1}\[m(L)b(L)\]=\\bar a(L)^{-1}\\bar b(L).
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0Ak%28L%29%3Da%28L%29%5E%7B-1%7Db%28L%29%3D%5Bm%28L%29a%28L%29%5D%5E%7B-1%7D%5Bm%28L%29b%28L%29%5D%3D%5Cbar%20a%28L%29%5E%7B-1%7D%5Cbar%20b%28L%29.%0A "
k(L)=a(L)^{-1}b(L)=[m(L)a(L)]^{-1}[m(L)b(L)]=\bar a(L)^{-1}\bar b(L).
")

One popular identification approach for the VARMA model is to use the
canonical echelon form parametrization, which involves zero and unity
restrictions in
![a(L)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;a%28L%29 "a(L)")
and
![b(L)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;b%28L%29 "b(L)").
By noting that the identification restrictions for the IRF in right
matrix fraction description (RMFD, for short) model,
i.e. ![k(L)=d(L)c(L)^{-1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;k%28L%29%3Dd%28L%29c%28L%29%5E%7B-1%7D "k(L)=d(L)c(L)^{-1}"),
are equivalent to those for the VARMA model after transposing, the
derivation is quite straightforward using the existing results for the
VARMA model (for details, see [Appendix
A](https://arxiv.org/pdf/2202.00310.pdf#section*.2)).

The structure of the RMFD model in echelon form is quite complicated,
and while the functions of this package `rmfd4dfm` code these
restrictions without the need for the user to trouble herself with the
implementation, one might still wonder what is the point of going
through the painstaking process of deriving and coding the restrictions
in the first place, as there are alternatives in the literature? To make
the point, let us introduce a popular alternative used in the literature
and compare it with our approach. One way to deal with the
identification and estimation of the DFM presented above is to stack the
dynamic factors into a
![r=q(s+1)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;r%3Dq%28s%2B1%29 "r=q(s+1)")-dimensional
vector

![
z_t = \\left(z_t^{\*'}, z\_{t-1}^{\*'}, \\cdots, z\_{t-s}^{\*'} \\right)',
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0Az_t%20%3D%20%5Cleft%28z_t%5E%7B%2A%27%7D%2C%20z_%7Bt-1%7D%5E%7B%2A%27%7D%2C%20%5Ccdots%2C%20z_%7Bt-s%7D%5E%7B%2A%27%7D%20%5Cright%29%27%2C%0A "
z_t = \left(z_t^{*'}, z_{t-1}^{*'}, \cdots, z_{t-s}^{*'} \right)',
")

and write Eq. (1) as

![
x_t = \\left(d_0,\\cdots,d_s\\right)z_t + \\xi_t = Dz_t + \\xi_t,
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0Ax_t%20%3D%20%5Cleft%28d_0%2C%5Ccdots%2Cd_s%5Cright%29z_t%20%2B%20%5Cxi_t%20%3D%20Dz_t%20%2B%20%5Cxi_t%2C%0A "
x_t = \left(d_0,\cdots,d_s\right)z_t + \xi_t = Dz_t + \xi_t,
")

which makes the model amenable to the principal components (PC) methods.
In the second step, the static factor process is modelled as a
VAR(![k](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;k "k"))
process

![
z_t - A_1 z\_{t-1} - \\cdots - A_k z\_{t-k} = A(L)z_t = B\\varepsilon_t,
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0Az_t%20-%20A_1%20z_%7Bt-1%7D%20-%20%5Ccdots%20-%20A_k%20z_%7Bt-k%7D%20%3D%20A%28L%29z_t%20%3D%20B%5Cvarepsilon_t%2C%0A "
z_t - A_1 z_{t-1} - \cdots - A_k z_{t-k} = A(L)z_t = B\varepsilon_t,
")

where
![A(L)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;A%28L%29 "A(L)")
is an
![r \\times r](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;r%20%5Ctimes%20r "r \times r")
VAR lag polynomial of order
![k](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;k "k"),
and
![B](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;B "B")
is an
![r\\times q](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;r%5Ctimes%20q "r\times q")
constant matrix. After the estimation of
![D](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;D "D")
and
![A(L)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;A%28L%29 "A(L)"),
the IRF
![k(L)=DA(L)^{-1}B](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;k%28L%29%3DDA%28L%29%5E%7B-1%7DB "k(L)=DA(L)^{-1}B")
can be constructed straightforwardly.

This alternative using PCs is straightforward from an estimation point
of view, but let us point two associated restrictive features. First,
note that the minimal number of identifying restrictions needed is
![r^2=q^2(s+1)^2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;r%5E2%3Dq%5E2%28s%2B1%29%5E2 "r^2=q^2(s+1)^2"),
which is higher than that of in our parametrization,
i.e. ![q^2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;q%5E2 "q^2"),
whenever
![s\>0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;s%3E0 "s>0")
[(Bai and Wang,
2012)](https://mpra.ub.uni-muenchen.de/38434/2/MPRA_paper_38434.pdf).
Second, we note that the reliable estimation of VAR on
![z_t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;z_t "z_t")
can be difficult if the lag order is misspecified, which is due to the
singularity of the covariance matrix of the innovation to
![z_t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;z_t "z_t"),
i.e.
![\\mathbb E (B\\varepsilon_t \\varepsilon_t'B')=B\\Sigma\_\\varepsilon B'](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbb%20E%20%28B%5Cvarepsilon_t%20%5Cvarepsilon_t%27B%27%29%3DB%5CSigma_%5Cvarepsilon%20B%27 "\mathbb E (B\varepsilon_t \varepsilon_t'B')=B\Sigma_\varepsilon B'")
has rank
![q](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;q "q"),
which is smaller than
![r=(s+1)q](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;r%3D%28s%2B1%29q "r=(s+1)q")
if
![s\>0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;s%3E0 "s>0")
[(Hörmann and Nisol,
2021)](https://onlinelibrary.wiley.com/doi/pdf/10.1111/jtsa.12568). On
the other hand, the law of motion for
![z_t^\*](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;z_t%5E%2A "z_t^*")
is standard non-singular VAR in our setup. These two examples make the
case that the “static” method is suboptimal whenever
![s\>0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;s%3E0 "s>0")
and warrant considering other alternatives, development of which is the
main goal of our paper.

## Data

The package makes available three different macroeconomic panels of data
measured at a monthly frequency. First is the data used by [Forni and
Gambetti
(2010)](http://pareto.uab.es/lgambetti/ForniGambettiJMEsecondRevSecondStage.pdf),
and the documentation can be accessed using command `?FG_data` in the R
console. Second and third data sets, `FRED_light` and `FRED_heavy`
originate from the FRED-MD data set, which is documented carefully in
[McCracken and Ng
(2015)](https://s3.amazonaws.com/real.stlouisfed.org/wp/2015/2015-012.pdf).
The difference between the data sets concerns the transformations to
obtain stationarity, while the time span is same in both data sets to
facilitate comparison to [Forni and Gambetti
(2010)](http://pareto.uab.es/lgambetti/ForniGambettiJMEsecondRevSecondStage.pdf).
The details can be accessed by `?FRED_light` or `?FRED_heavy`.

## Replication of the Monetary Policy Example

The empirical example compares the DFMs estimated by [Forni and Gambetti
(2010)](http://pareto.uab.es/lgambetti/ForniGambettiJMEsecondRevSecondStage.pdf)
to that developed in our paper. The idea is to assess the effects of
monetary policy to key macroeconomic variables via structural DFMs,
which guards against the omitted variable bias, to which the SVAR
methods can be more susceptible.

The following code snippets show how to replicate the results given in
the empirical section of our paper. First, we order the variables of
interest in the structural analysis variables first in the data matrix,
i.e. industrial production (INDPRO, the mnemonic in [McCracken and Ng,
2015](https://s3.amazonaws.com/real.stlouisfed.org/wp/2015/2015-012.pdf)),
consumer price index (CPIAUCSL), federal funds rate (FEDFUNDS), and real
Swiss/US exchange rate (EXSZUSx), in this order. In the structural
analysis with recursive identification, the ordering is important, as it
determines the short-run restrictions imposed on the structural impact
multiplier matrix. (For a brief and to the point summary of the
structural identification using DFMs, see [Section
9](http://www.barigozzi.eu/MB_DF_lecture_notes_online.pdf#page=58) of
the lecture notes written by [Matteo
Barigozzi](http://www.barigozzi.eu/Home.html).) Additionally, one must
include a vector `int_ix` in the data object `FRED_heavy` coding the
positions of the variables of interest, which in this case is simply a
sequence `1:4`.

``` r
library("rmfd4dfm")
# code the positions of the variables of interest
int_vars_fred <- c("INDPRO", "CPIAUCSL", "FEDFUNDS", "EXSZUSx")
int_ix <- sapply(int_vars_fred, function(x) which(names(FRED_heavy$df)==x))
# re-organize data matrix s.t. the variables of interest are ordered first
perm_ix <- c(int_ix, (1:ncol(FRED_heavy$df))[-int_ix])
FRED_heavy$df <- FRED_heavy$df[,perm_ix]
FRED_heavy$trans_ix <- FRED_heavy$trans_ix[perm_ix]
FRED_heavy$int_ix <- 1:4
```

Second, we need to determine the static factor dimension. This is needed
for an estimate of the lag structure in eqs. (1)–(2), while in the
“static” method it defines the dimension of
![z_t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;z_t "z_t").
To this end, we use functions `baingcriterion` and `abc_crit`, which
implement the tests developed in [Bai and Ng
(2002)](https://www.ssc.wisc.edu/~bhansen/718/BaiNg2002.pdf) and
[Alessi, Barigozzi and Capasso
(2009)](https://dipot.ulb.ac.be/dspace/bitstream/2013/54139/1/RePEc_eca_wpaper_2009_023.pdf),
respectively.

``` r
# Estimate the number of static factors ####
baing_cr <- baingcriterion(FRED_heavy$df, rmax = 25)$IC[1:2]
abc_cr <- replicate(50, abc_crit(FRED_heavy$df, kmax = 25))
get_mode <- function(x) unique(x)[which.max(colSums(sapply(unique(x), function(z) x %in% z)))]
r_hats <- c("ICp1" = baing_cr[1],
            "ICp2" = baing_cr[2],
            "ABCp1" = get_mode(unlist(abc_cr[1,])),
            "ABCp2" = get_mode(unlist(abc_cr[2,]))
            )
```

Upon fixing the state dimension of the state space representation of
eqs. (1)–(2), estimation of the DFM can be carried out. Here we fix the
static factor dimension to 8. In the RMFD-E parametrization of the DFM
and for a given state dimension, it is likely that one has many model
specifications available from which to choose “the best model”.
Therefore, the function `do_everything_rmfd` estimates all the model
alternatives that are consistent with the model selection criteria given
in [Section 3.3.](https://arxiv.org/pdf/2202.00310.pdf#subsection.3.3)
and chooses the one that minimizes Bayesian information criterion. For
this model, the function then estimates confidence intervals for the IRF
using block bootstrap. Finally, the function returns the impulse
responses of the variables to the third shock, which is standardized to
have an immediate impact of 0.5. For responses to other shocks and/or
size, one can include an index vector of length two, which specifies the
shock of interest and normalization constant, in the data object
specified as the first argument of `do_everything_rmfd` by
`df$shock_ix`. For a more flexible estimation setup, the user can
consult function `estim_wrap`, which produces an estimate of the
non-structural IRF
![k(L)=d(L)c(L)^{-1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;k%28L%29%3Dd%28L%29c%28L%29%5E%7B-1%7D "k(L)=d(L)c(L)^{-1}")
for a given model structure. The functions `do_everything_fglr` and
`do_everything_svar` perform the same task for the static DFM and
structural VAR, respectively.

``` r
# D-DFM
est_obj <- do_everything_rmfd(df = FRED_heavy,
                              r = 8,
                              h = 50,
                              nrep = 500,
                              conv_crit = 1e-5,
                              ci = 0.68,
                              init_rep = 1,
                              verbose = TRUE)
# S-DFM
est_fglr <- do_everything_fglr(df = FRED_heavy,
                               r = 8,
                               k = 2,
                               h = 50,
                               nrep = 500,
                               ci = 0.68)

# SVAR
svar_irf <- do_everything_svar(df = FRED_heavy,
                               p = 9,
                               nrep = 500,
                               h = 50,
                               ci = 0.68)
```

Finally, the `plot_irfs` function returns the IRFs along with the
possible error bands in a nice plot. To gather the IRF plots into one
figure, which can then be exported to pdf, for example, we use
`gridExtra::marrangeGrob`.

``` r
p1 <- plot_irfs(est_obj$irf, int_vars_fred, "D-DFM")
p2 <- plot_irfs(est_fglr, int_vars_fred, "S-DFM", label_y = FALSE)
p3 <- plot_irfs(svar_irf, int_vars_fred, "SVAR", label_y = FALSE)
plot1 <- gridExtra::marrangeGrob(c(p1,p2,p3), nrow = 4, ncol = 3, top = NULL)
```
