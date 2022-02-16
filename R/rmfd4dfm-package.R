#' @title Estimation of Impulse-Response Functions with Dynamic Factor Models
#'
#' @description The \code{rmfd4dfm} package contains functions to estimate
#' impulse-response functions with two specifications of a dynamic factor model,
#' the so-called dynamic and static form of the DFM. This \code{R} package is
#' an accompanying resource for replication purposes of the results given in
#' the paper \emph{Estimation of Impulse-Response Functions with Dynamic Factor
#' Models: A New Parametrization} available at
#' \url{https://arxiv.org/pdf/2202.00310.pdf}.
#'
#' @author Juho Koistinen <juho.koistinen@@helsinki.fi>
#'
#' @docType package
#' @name rmfd4dfm
#' @importFrom magrittr %>%
#' @importFrom corpcor pseudoinverse
#'
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
