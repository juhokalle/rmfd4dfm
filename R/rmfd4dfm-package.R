#' @title Estimation of Impulse-Response Functions with Dynamic Factor Models: A New Parametrization
#'
#' @description
#' This is the associated R-package to the paper with the above title available at
#' \url{https://arxiv.org/pdf/2202.00310.pdf}.
#'
#' Abstract:
#'
#' We propose a new parametrization for the estimation and identification
#' of the impulse-response functions (IRFs) of dynamic factor models
#' (DFMs). The theoretical contribution of this paper concerns the problem
#' of observational equivalence between different IRFs, which implies
#' non-identification of the IRF parameters without further restrictions.
#' We show how the minimal identification conditions proposed by Bai and Wang (2015)
#' are nested in the proposed framework and can be further augmented
#' with overidentifying restrictions leading to efficiency gains. The
#' current standard practice for the IRF estimation of DFMs is based
#' on principal components, compared to which the new parametrization
#' is less restrictive and allows for modelling richer dynamics. As the
#' empirical contribution of the paper, we develop an estimation method
#' based on the EM algorithm, which incorporates the proposed identification
#' restrictions. In the empirical application, we use a standard high-dimensional
#' macroeconomic dataset to estimate the effects of a monetary policy
#' shock. We estimate a strong reaction of the macroeconomic variables,
#' while the benchmark models appear to give qualitatively counterintuitive
#' results. The estimation methods are implemented in the accompanying
#' R package.
#'
#' @references Bai, J., & Wang, P. (2015). Identification and Bayesian estimation
#' of dynamic factor models. Journal of Business & Economic Statistics, 33(2), 221-240.
#'
#' Koistinen, J., & Funovits, B. (2022). Estimation of Impulse-Response Functions with
#' Dynamic Factor Models: A New Parametrization. arXiv preprint arXiv:2202.00310.
#'
#' @section Dependencies:
#'
#' This package depends heavily on the \strong{rationalmatrices} and the \strong{RLDM} packages
#' developed by Wolfgang Scherrer and Bernd Funovits.
#' Since both packages may be subject to breaking changes, all necessary functionalities are
#' extracted in order to provide stable code with minimal external dependencies.
#' Importantly, please reference the \strong{rationalmatrices} and the \strong{RLDM}
#' packages should you use their functionalities and not this package.
#'
#' @section Author(s):
#'
#' Bernd Funovits, Juho Koistinen, and Wolfgang Scherrer
#'
#' Maintainer: <juho.koistinen@@helsinki.com>
#'
#'
#' @docType package
#' @name rmfd4dfm
#' @importFrom magrittr %>%
#' @importFrom corpcor pseudoinverse
#' @useDynLib rmfd4dfm, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' @title Pipe operator
#'
#' @description
#' Re-export pipe operator \code{\%>\%} to turn function composition into a series of imperative statements.
#' For more extensive description, see function \code{`\%>\%`} in package \emph{magrittr}.
#'
#' @importFrom magrittr %>%
#' @name %>%
#' @rdname pipe
#' @export
#' @param lhs First argument of the function of the right-hand-side of the pipe operator.
#' @param rhs Function whose first argument is given by the left-hand-side argument \code{lhs} of the pipe operator.
#'
NULL

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))
