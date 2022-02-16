#' Append file name to path
#'
#' Path is in enclosing environment!
#' Helper for running scripts.
#'
#' @param path file path to append to \code{path}
#'
#' @return New file path
#' @export
#'
#' @examples
#' path = "../local_data_ukko2/"
#' pap = pap_factory(path)
#' pap("myfile.whatever")
pap_factory <- function(path){
  function(str){
    paste0(path, str)
  }
}

#' Obtain a balanced and possibly imputed and standardized panel of data
#'
#' This convenience function allows the user to manipulate the supplied data matrix
#' with respect to
#' \itemize{
#'   \item{the sample span;}
#'   \item{reordering of the variables;}
#'   \item{imputation of outliers; and}
#'   \item{standardization.}
#' }
#'
#' @param df \eqn{T x (N+1)} data panel with the first column specifying the time stamp,
#' \eqn{T} giving the number of observations and \eqn{N} the number of variables
#' @param int_vars vector of strings specifying the variables of interest, which are ordered first
#' @param start_date \link{ymd} object specifying the start date of the sample
#' @param end_date \link{ymd} object specifying the end date of the sample
#' @param impute Are the outliers removed and replaced with an imputed value
#' (for the imputation scheme see \link[fbi]{tw_apc})? Defaults to TRUE
#' @param sdize Are the variables standardized to have zero mean and unit standard deviation? Defaults to FALSE
#' @param trans_ix vector containing the transformation codes of the variables
#'
#' @return list of components
#' \item{df}{A \eqn{T x N} data frame corresponding transformed data}
#' \item{date}{vector of \link{ymd} components giving the sample span}
#' \item{discarded}{A list of discarded variables containing missing values}
#' \item{sd_mat}{If \code{sdize = TRUE}, a matrix containing standard deviations of the variables on the diagonal
#' is also supplied}
#' \item{trans_ix}{vector of reordered transformation codes if \code{int_vars} is supplied}
#'
#' @keywords internal
prune_data <- function(df, int_vars = NULL,
                       start_date = NULL, end_date = NULL,
                       impute = TRUE, sdize = TRUE,
                       trans_ix = NULL){

  nobs <- nrow(df)
  if(!sapply(df, class)[1]=="Date") stop("Specify the first column as Date variable")
  if(names(df)[1]!="date") names(df)[1] <- "date"
  if(is.null(start_date)) start_date <- df$date[1]
  if(is.null(end_date)) end_date <- df$date[nobs]
  dd <- df$date[2]-df$date[1] # define the sampling frequency
  frq <- if(dd < 32){"month"} else if(dd < 93){"quarter"} else{"year"}
  time_span <- seq(start_date, end_date, frq)
  # drop observations outside of the sample span
  df = df %>% filter(date %in% time_span)
  date_vec <- df$date # save date variable and drop it from the data matrix
  df = df %>% select(-date)
  ix_na <- apply(df, 2, function(x) sum(is.na(x)))
  var_names <- names(df)

  # balance panel by dropping any variable containing missing vars
  if(any(ix_na>0)){
    discarded <- data.frame(NA_count = ix_na[ix_na!=0])
    df <- df[,ix_na==0]
    var_names <- var_names[!var_names %in% rownames(discarded)]
    trans_ix = trans_ix[var_names]
    nvar <- ncol(df)
    out <- list(discarded = discarded)
  }
  # data imputation
  if(impute){
    # helper function: detects outliers á la McCracken & Ng (JBES, 2016)
    rm_outlr <- function(x){
      devx <- abs(x-median(x, na.rm = T))
      x[devx > 10*IQR(x, na.rm = T)] <- NA
      return(x)
    }
    df <- df %>% mutate_if(is.numeric, rm_outlr)
    df <- tw_apc(df, kmax = 20)$data %>% as_tibble() # data imputation, for details see package 'fbi'
    colnames(df) <- var_names # rename the resulting data matrix
  }
  # standardization
  if(sdize){
    out$sd_mat <- df %>% select_if(is.numeric) %>% apply(2, sd) %>% diag
    df <- df %>% mutate_if(is.numeric, ~ (.x-mean(.x))/sd(.x))
  }
  # reordering
  if(!is.null(int_vars)){
    int_ix <- sapply(int_vars, function(x) which(var_names==x))
    trans_ix <- trans_ix[c(int_ix, (1:nvar)[-int_ix])]
    df <- df[ , c(int_ix, (1:nvar)[-int_ix])]
  }

  out$df <- df
  out$date <- date_vec
  out$trans_ix <- trans_ix

  return(out)
}

#' Get sample quantiles of the autocorrelation coefficients of the variables
#'
#' This function prints the percent quantiles of the autocorrelation coefficients
#' of the variables included in the data set. This is a simple way to
#' summarize the degree of non-stationarity and persistence of the
#' variables. The number of lags at which autocorrelations are calculated
#' and the grid of probabilities for the percent quantiles can be defined.
#'
#'
#' @param df \eqn{T x N} data matrix
#' @param lag.acf maximum number of lags at which the autocorrelations are calculated
#' @param qnt_grid vector of probabilities determining the sample percent quantiles
#'
#' @return Quantile table giving the distribution of autocorrelation coefficients at different lags
#'
#' @keywords internal
get_acf_qnt <- function(df, lag.acf = 8, qnt_grid = c(.05,.25,.5,.75,.95)){

  acf_obj <- apply(df, 2, function(x) acf(x, plot = FALSE)$acf)
  qnt_mat <- apply(acf_obj[(1:lag.acf)+1,], 1,
                   function(x) quantile(abs(x), probs=qnt_grid))
  return(qnt_mat)
}
