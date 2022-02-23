#' Append file name to path
#'
#' Path is in enclosing environment!
#' Helper for running scripts.
#'
#' @param path file path to append to \code{path}
#'
#' @return New file path
#' @keywords internal
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

#' @title Obtain a block bootstrap data sample
#'
#' @description
#' \code{X_boot} performs a simple block bootstrap for the data matrix.
#' This is a direct copy of the code used in Forni and Gambetti (2010)
#' available at \url{http://pareto.uab.es/lgambetti/ReplicaForniGambettiJME.zip}.
#'
#' @param X data matrix
#' @param L determines the block size
#'
#' @return bootstrapped data matrix of the same dimensions as the original one
#'
#' @keywords internal
#'
X_boot <- function(X, L){

  X <- as.matrix(X)
  nobs <- nrow(X)
  nvar <- ncol(X)
  K <- floor(nobs/L)
  blks <- ceiling(stats::runif(20)*K)
  X_boot <- matrix(0, nrow = K*L, ncol = nvar)
  for(i in 1:K){
    X_boot[((i-1)*L + 1):(i*L),] <- X[((blks[i]-1)*L + 1):(blks[i]*L),]
  }
  return(X_boot)
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
#' @param start_date \link[lubridate]{ymd} object specifying the start date of the sample
#' @param end_date \link[lubridate]{ymd} object specifying the end date of the sample
#' @param impute Are the outliers removed and replaced with an imputed value
#' (for the imputation scheme see \link[fbi]{tw_apc})? Defaults to TRUE
#' @param sdize Are the variables standardized to have zero mean and unit standard deviation? Defaults to FALSE
#' @param trans_ix vector containing the transformation codes of the variables
#'
#' @return list of components
#' \item{df}{A \eqn{T x N} data frame corresponding transformed data}
#' \item{date}{vector of \link[lubridate]{ymd} components giving the sample span}
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
  df = df %>% dplyr::filter(date %in% time_span)
  date_vec <- df$date # save date variable and drop it from the data matrix
  df = df %>% dplyr::select(-date)
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
      devx <- abs(x-stats::median(x, na.rm = T))
      x[devx > 10*stats::IQR(x, na.rm = T)] <- NA
      return(x)
    }
    df <- df %>% dplyr::mutate_if(is.numeric, rm_outlr)
    df <- fbi::tw_apc(df, kmax = 20)$data %>% tibble::as_tibble() # data imputation, for details see package 'fbi'
    colnames(df) <- var_names # rename the resulting data matrix
  }
  # standardization
  if(sdize){
    out$sd_mat <- df %>% dplyr::select_if(is.numeric) %>% apply(2, stats::sd) %>% diag
    df <- df %>% dplyr::mutate_if(is.numeric, ~ (.x-mean(.x))/stats::sd(.x))
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

  acf_obj <- apply(df, 2, function(x) stats::acf(x, plot = FALSE)$acf)
  qnt_mat <- apply(acf_obj[(1:lag.acf)+1,], 1,
                   function(x) stats::quantile(abs(x), probs=qnt_grid))
  return(qnt_mat)
}

#' @title Left and right multiply IRFs with constant matrices
#'
#' @description
#' \code{irf_x} is a helper function rotating arrays in the third
#' dimension of the given array.
#'
#' @export
#'
#' @param irf \code{n x q x (h + 1)} array containing the impulse responses,
#' where \eqn{n, q} and \eqn{h} are the cross-sectional dimension, number of
#' dynamic factors, and IRF horizon, resepectively.
#' @param pre_mat matrix used in the left multiplication of the IRF
#' @param post_mat matrix used in the right multiplication of the IRF
#'
#' @return An array rotated in the third dimension.
#'
#' @examples
#' # create a random array corresponding to the impulse response function
#' irf_ex <- stats::rnorm(4*4*4) %>% array(dim = c(4,4,4))
#' # assume the recursive identification via cholesky is given as
#' chol_ex <- stats::rnorm(16) %>% matrix(4,4) %>% tcrossprod() %>% chol() %>% t()
#' # then the recursive identification of the IRF can be obtained as
#' (sirf_ex <- irf_x(irf_ex, post_mat = solve(irf_ex[,,1])%*%chol_ex) %>% zapsmall())
#'
irf_x <- function(irf, pre_mat = NULL, post_mat = NULL){

  if(is.null(pre_mat)) pre_mat <- diag(dim(irf)[1])
  if(is.null(post_mat)) post_mat <- diag(dim(irf)[2])
  new_irf_list <- purrr::map(1:dim(irf)[3], ~ pre_mat %*% irf[,,.x] %*% post_mat)
  new_irf_array <- purrr::reduce(new_irf_list, .f = function(x,y){dbind(d=3,x,y)})
  return(new_irf_array)
}

#' @title Cumulate IRFs for variables in differences
#'
#' @description
#' This function cumulates IRFs such that variables in first and second differences
#' are cumulated once and twice, respectively.
#'
#' @param irf_array \eqn{(N x q x (h+1))} array containing the impulse response function, with \eqn{N} denoting
#' the number of variables, \eqn{q} the number of shocks, \eqn{h} the IRF horizon
#' @param trans_ix vector of transformation codes according in McCracken & Ng (JME, 2016) format
#'
#' @return \eqn{(N x q x (h+1))} array containing the cumulated IRFs
#'
#' @keywords internal
#'
cum_irf <- function(irf_array, trans_ix){

  firstdiff <- trans_ix %in% c(2,3,5,6,7)
  seconddiff <- trans_ix %in% c(3,6,7)

  if(any(firstdiff)){
    irf_array[firstdiff, ,] <- irf_array[firstdiff,,,drop=FALSE] %>%
      apply(MARGIN = c(1,2), FUN = cumsum) %>%
      aperm(perm = c(2,3,1))
  }
  if(any(seconddiff)){
    irf_array[seconddiff, ,] <- irf_array[seconddiff,,,drop=FALSE] %>%
      apply(MARGIN = c(1,2), FUN = cumsum) %>%
      aperm(perm = c(2,3,1))
  }
  return(irf_array)
}

#' @title Finalize IRF for the plots
#'
#' @description
#' \code{finalize_irf} transforms the structural impulse responses such
#' that they reflect the responses of the endogenous variables
#' to an exogenous shock of custom size, determined by the user.
#' Additionally, it accumulates impulse responses of differenced
#' variables such that they reflect the responses in levels.
#' Finally, variables in logs, log-differences or in relative differences
#' are multiplied by 100 such that they are transformed to percentages.
#' This function is intended for a single shock identification,
#' i.e. the output is a matrix with the columns containing the responses
#' of the variables to the shock.
#'
#' @param irf_obj array containing the structural impulse responses
#' @param shock_size size of the shock
#' @param norm_id vector of length two, which code the variable corresponding
#' to the shock of interest and its ordering in the system
#' @param trans_ix vector containing the transformation codes
#' @param int_vars vector coding the interest variables
#'
#' @return A matrix containing the IRFs as columns
#'
#' @keywords internal
#'
finalize_irf <- function(irf_obj, shock_size, norm_id, trans_ix, int_vars){

  dim_in <- dim(irf_obj)[2]
  irf_obj <- irf_obj/irf_obj[norm_id[1],norm_id[2],1]*shock_size
  irf_obj <- cum_irf(irf_obj, trans_ix)
  irf_obj[trans_ix%in%c(4,5,6,7),,] <- irf_obj[trans_ix%in%c(4,5,6,7),,]*100
  irf_mat <- sapply(int_vars, function(x) irf_obj[x, norm_id[2],])

  return(irf_mat)
}

#' @title Plot impulse response functions
#'
#' @description
#' \code{plot_irfs} plots impulse response functions for a given array
#' containing the impulse responses and possibly also error bands. This function
#' returns a list of figures containing IRFs per variable.
#'
#' @export
#'
#' @param irf_arr Matrix or array of dimension \eqn{h+1 x q} (w/o error bands)
#' or \eqn{(h+1 x 3 x q)} (w/ error bands), where \eqn{h} is the number of periods
#' in IRF and \eqn{q} is the number of variables of interest
#' @param var_names Vector of variable names of the \eqn{q} impacted variables
#' @param plot_title Title of the figure
#' @param label_y Should the variable name be printed next to y-axis?
#'
#' @return list containing the figures
#' @examples
#' # Create a random IRF array: in case a matrix, rows denote the
#' # periods ahead while the columns are the reactions of the variables
#' h <- 24 # IRF horizon
#' q <- 4 # number of variables
#' irf_ex <- matrix(stats::rnorm((h+1)*q), h+1, q)
#' irf_ex <- apply(irf_ex, 2, cumsum)
#' p1 <- plot_irfs(irf_arr = irf_ex, var_names = letters[1:q])
#' # in case the supplied IRF element is array, the interpretation is as follows
#' # 1st dimension: length of the IRF horizon, 2nd dimension contain the
#' # lower error band, point estimate and upper error band estimates,
#' # 3rd is the of variables
#' irf_arr <- c(irf_ex-sd(irf_ex)*1.645, irf_ex, irf_ex+sd(irf_ex)*1.645) %>%
#'   array(dim= c(h+1, q, 3)) %>%
#'   aperm(c(1,3,2))
#' p2 <- plot_irfs(irf_arr, var_names = letters[1:4], label_y = FALSE)
#' \dontrun{
#' gridExtra::marrangeGrob(c(p1,p2), nrow = 4, ncol = 2)}
plot_irfs <- function(irf_arr, var_names, plot_title = "", label_y = TRUE){

  h <- dim(irf_arr)[1]
  if(length(dim(irf_arr))==2){
    q <- dim(irf_arr)[2]
    irf_arr <- c(rep(0, h*q), irf_arr, rep(0, h*q)) %>%
      array(dim= c(h, q, 3)) %>%
      aperm(c(1,3,2))
  }
  q <- dim(irf_arr)[3]
  dim(irf_arr) <- c(h, q*3)
  colnames(irf_arr) <- c("low", "point", "up") %>%
    expand.grid(var_names) %>%
    apply(1, function(x) paste(x, collapse="_"))
  irf_arr <- tibble::as_tibble(dplyr::bind_cols(irf_arr, "period" = 0:(h-1)))
  irf_plot <- vector(mode = "list", length = q)

  for(jj in 1:q){
    irf_arr %>%
      dplyr::select(dplyr::contains(var_names[jj]), "period") %>%
      ggplot2::ggplot(ggplot2::aes_string(x = "period",
                                          y = paste("point", var_names[jj], sep="_"))) +
      ggplot2::geom_line() +
      ggplot2::geom_ribbon(ggplot2::aes_string(ymin = paste("low", var_names[jj], sep="_"),
                                               ymax = paste("up", var_names[jj], sep="_")),
                           alpha = 0.5, fill= "lightblue") +
      ggplot2::scale_x_continuous(expand = c(0,0)) +
      ggplot2::labs(title = ifelse(jj==1, plot_title, ""),
                    x = "",
                    y = ifelse(label_y, var_names[jj], "")) +
      ggplot2::theme(legend.position = "none",
                     plot.title = ggplot2::element_text(hjust = 0.5),
                     plot.margin = ggplot2::unit(rep(0.1,4),"cm"),
                     panel.background = ggplot2::element_rect(fill='transparent'), #transparent panel bg
                     plot.background = ggplot2::element_rect(fill='transparent', color=NA), #transparent plot bg
                     panel.grid.major = ggplot2::element_blank(), #remove major gridlines
                     panel.grid.minor = ggplot2::element_blank(), #remove minor gridlines
                     panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=.5)
      ) -> irf_plot[[jj]]
  }
  return(irf_plot)
}
