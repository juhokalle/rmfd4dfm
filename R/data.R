#' @title High-dimensional macroeconomic data set used in Forni and Gambetti (2010)
#'
#' @description
#' This data is used by Forni and Gambetti (2010) in their monetary policy
#' study using dynamic factor models. The data spans the period between
#' March 1973 and November 2007 and includes 112 macroeconomic variables.
#' The data is included in this package for purposes of comparison to the
#' main data set used in the empirical exercise, FRED-MD. The method of obtaining
#' this data set from the raw data is described in the script \code{data-raw/DATASET.R}.
#'
#' @format A list with the following elements
#' \describe{
#'   \item{df}{\eqn{T x n} data frame where \eqn{T} is the time-series dimension,
#'   and \eqn{n} is the cross-sectional dimension}
#'   \item{data}{\code{Date} variable giving the time span of the data}
#'   \item{trans_ix}{Variable transformation codes with the interpretation as in
#'   McCracken and Ng (2016)}
#' }
#'
#'
#' @source \url{http://pareto.uab.es/lgambetti/ReplicaForniGambettiJME.zip}
#'
#' @references Forni, M., & Gambetti, L. (2010). The dynamic effects of monetary policy:
#' A structural factor model approach. Journal of Monetary Economics, 57(2), 203-216.
#'
#' McCracken, M. W., & Ng, S. (2016). FRED-MD: A monthly database for macroeconomic research.
#' Journal of Business & Economic Statistics, 34(4), 574-589.
#'
"FG_data"

#' @title FRED-MD (December 2021 vintage) data set with heavy transformations
#'
#' @description
#' This data is taken from FRED-MD database, with the sample span corresponding to
#' that of \code{FG_data}, i.e. March 1973 to November 2007. The database consists
#' of 127 variables, out of which two are discarded due to many missing variables, and
#' the corresponding variable names can be accessed by \code{FRED_heavy$discarded}.
#' The data is transformed stationary according to the methodology given in
#' McCracken and Ng. The method of obtaining this data set from the raw data
#' is described in the script \code{data-raw/DATASET.R}.
#'
#' @format A list with the following elements
#' \describe{
#'   \item{discarded}{data frame listing the discarded variables and the number of
#'   missing values under the time span under consideration}
#'   \item{df}{\eqn{T x n} data frame where \eqn{T} is the time-series dimension,
#'   and \eqn{n} is the cross-sectional dimension}
#'   \item{data}{\code{Date} variable giving the time span of the data}
#'   \item{trans_ix}{Variable transformation codes with the interpretation as in
#'   McCracken and Ng (2016)}
#' }
#'
#'
#' @source \url{https://files.stlouisfed.org/files/htdocs/fred-md/monthly/2021-12.csv}
#'
#' @references McCracken, M. W., & Ng, S. (2016). FRED-MD: A monthly database for macroeconomic research.
#' Journal of Business & Economic Statistics, 34(4), 574-589.
#'
"FRED_heavy"

#' @title FRED-MD (December 2021 vintage) data set with light transformations
#'
#' @description
#' This data is taken from FRED-MD database, with the sample span corresponding to
#' that of \code{FG_data}, i.e. March 1973 to November 2007. The database consists
#' of 127 variables, out of which two are discarded due to many missing variables, and
#' the corresponding variable names can be accessed by \code{FRED_light$discarded}.
#' The data is transformed differenly to the dataset \code{FRED_heavy} such that
#' the transformations would correspond more closely to those given in Forni and Gambetti (2010).
#' The method of obtaining this data set from the raw data and the transformation
#' codes are described in the script \code{data-raw/DATASET.R}.
#'
#' @format A list with the following elements
#' \describe{
#'   \item{discarded}{data frame listing the discarded variables and the number of
#'   missing values under the time span under consideration}
#'   \item{df}{\eqn{T x n} data frame where \eqn{T} is the time-series dimension,
#'   and \eqn{n} is the cross-sectional dimension}
#'   \item{data}{\code{Date} variable giving the time span of the data}
#'   \item{trans_ix}{Variable transformation codes with the interpretation as in
#'   McCracken and Ng (2016)}
#' }
#'
#'
#' @source \url{https://files.stlouisfed.org/files/htdocs/fred-md/monthly/2021-12.csv}
#'
#' @references Forni, M., & Gambetti, L. (2010). The dynamic effects of monetary policy:
#' A structural factor model approach. Journal of Monetary Economics, 57(2), 203-216.
#'
#' McCracken, M. W., & Ng, S. (2016). FRED-MD: A monthly database for macroeconomic research.
#' Journal of Business & Economic Statistics, 34(4), 574-589.
#'
"FRED_light"
