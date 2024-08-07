% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/em_algo.R
\name{smoothed_moments8}
\alias{smoothed_moments8}
\title{Obtain Moment Matrices for the EM Algorithm Using Kalman Filter/Smoother}
\usage{
smoothed_moments8(stsp_mod, Sigma, data_wide, only_ll = FALSE)
}
\arguments{
\item{stsp_mod}{\code{stspmod} object containing the state space model used for smoothing and the idiosyncratic noise component}

\item{Sigma}{Matrix of dimension \code{dim_in x dim_in}. All \code{dim_in^2} elements are saved, symmetry is not taken into account.}

\item{data_wide}{Matrix of dimension \code{dim_out x n_obs}}

\item{only_ll}{return only the log-likelihood value calculated using KF}
}
\description{
For details, see section 3.2. and Appendix C in \emph{Estimation of Impulse-Response
Functions with Dynamic Factor Models: A New Parametrization}
available at \url{https://arxiv.org/pdf/2202.00310.pdf}.
}
\seealso{
Watson, M. W., & Engle, R. F. (1983). Alternative algorithms for the
estimation of dynamic factor, mimic and varying coefficient regression models.
Journal of Econometrics, 23(3), 385-400.
}
\keyword{internal}
