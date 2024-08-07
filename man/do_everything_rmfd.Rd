% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/do_everything.R
\name{do_everything_rmfd}
\alias{do_everything_rmfd}
\title{Estimate IRFs and associated bootstrap intervals for D-DFM}
\usage{
do_everything_rmfd(
  df,
  r,
  h,
  nrep,
  conv_crit = 1e-05,
  ci = 0.68,
  init_rep = 0,
  verbose = FALSE,
  mod_str = NULL
)
}
\arguments{
\item{df}{a list with the following elements:
\enumerate{
\item \code{df}, data in a \eqn{T x n} matrix;
\item \code{int_ix}, an index vector coding the columns corresponding to the variables of interest in the data;
\item \code{trans_ix}, an index vector giving the variable transformation codes as in McCracken and Ng (2016)
\item \code{shock_ix}, specify the shock of interest and the normalization at impact, defaults to \code{c(3,0.5)}
}}

\item{r}{state dimension of the D-DFM in state space model}

\item{h}{estimate \eqn{h}-step ahead IRFs}

\item{nrep}{number of replications in the block bootstrap procedure, set to zero for no bootstrapping}

\item{conv_crit}{convergence criterion of the EM algorithm, the default is 1e-5}

\item{ci}{confidence interval level for the IRFs, defaults to 0.68}

\item{init_rep}{bootstrap replications for the initial value procedure,
for details, see \code{\link{boot_init}}}

\item{verbose}{logical, print the estimation process?}

\item{mod_str}{optional list for a fixed model structure, containing at least the Kronecker index vector
and also possibly the orders of \eqn{c(z)} and \eqn{d(z)} as a second element}
}
\value{
a list with elements
\item{irf}{an array corresponding to the recursively identified IRF}
\item{conv_crit}{convergence criteria for the model minimizing BIC}
\item{model_selection}{table containing the model selection criteria, returned only if
the number of estimated models is larger than one}
}
\description{
\code{do_everything_rmfd} is a user-friendly function for implementing the estimation
and identification methodology of the D-DFM model presented in sections 2.2.
and 3, as well as replicating the empirical exercise of Section 5 in \emph{Estimation
of Impulse-Response Functions with Dynamic Factor Models: A New Parametrization},
available at \url{https://arxiv.org/pdf/2202.00310.pdf}.
}
\details{
The structure of the function is as follows:
\enumerate{
\item Given the state dimension of the D-DFM in state space format, which can be estimated
using functions \code{\link{abc_crit}} and \code{\link{baingcriterion}}, for example,
the function returns the set of model structures consistent with the model selection
criteria 3--5 given in Section 3.3 of the paper via the function \code{\link{admissible_mods}}.
Alternatively, the user can estimate a fixed model structure using the argument \code{mod_str} without
running the estimation algorithm across a set of model alternatives.
\item The function estimates a set of feasible model(s) using the function \code{\link{estim_wrap}}
and returns the statistics of the model selection criteria. The best model candidate is chosen
according to the BIC statistic.
\item The function estimates the bootstrapped confidence intervals for IRFs of the recursively identified
model with the smallest BIC value. In the case of no bootstrapping (\code{nrep=0}), the function returns
all the estimated IRFs, while with bootstrapping the function returns only the point estimates of the
IRF corresponding to the best model candidate determined using BIC.
}

For obtaining the qualified models and an estimation example, we refer the reader
to \code{\link{admissible_mods}} and \code{\link{estim_wrap}}. The running time of the algorithm
is determined by the size of the data matrix, complexity of the estimated model,
level of the convergence criterion, and the number of bootstrap draws.
For a standard data set (\eqn{n~100}, \eqn{T~300}), simple model structure,
a lenient convergence criterion (\code{conv_crit=1e-3}) and a small number of
bootstrap draws (\code{nrep=100}), the function should be finished in approximately 60 minutes.

The estimated models are identified recursively using Cholesky decomposition of the residual covariance
matrix. As a default, the function returns the impulse responses to the third shock, the size of which
is normalized to 0.5 on impact. For changing the shock of interest (i.e. not the third),
and the normalization constant, the user should include store a vector of length 2 as \code{df$shock_ix}
with the position of the shock of interest as the first element and the
normalization constant as the second element. For a more flexible setup, such as identifying
more than one shock or obtaining responses to many variables, the user can use function
\code{\link{estim_wrap}}, which returns the non-structural IRF.
}
\examples{
\dontrun{
# the following example takes around 60 minutes to finish
FRED_heavy$int_ix <- c(5,96,77,101) # remember to index the variables of interest
est_obj <- do_everything_rmfd(FRED_heavy, r = 8, h = 48, nrep = 100,
  conv_crit = 1e-3, ci = .68, init_rep = 0, verbose = TRUE)}

}
\seealso{
McCracken, M. W., & Ng, S. (2016). FRED-MD: A monthly database for macroeconomic research.
Journal of Business & Economic Statistics, 34(4), 574-589.
}
