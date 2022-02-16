#' @title Two-step estimation method of S-DFM using VAR
#'
#' @description
#' \code{DfmRawImp} estimates non-structural IRFs for S-DFM according to section 2.3.1 in
#' \emph{Estimation of Impulse-Response Functions with Dynamic Factor Models:
#' A New Parametrization} available at \url{https://arxiv.org/pdf/2202.00310.pdf}.
#'
#' @param X \eqn{T x n} data matrix, where \eqn{T} and \eqn{n}
#'  are the time series and cross-sectional dimensions, respectively
#' @param q The dynamic factor dimension
#' @param r Number of static factors, must be larger than or equal to \eqn{q}
#' @param k Number of lags in the VAR for static factors
#' @param h estimate IRF \eqn{h}-period ahead
#'
#' @return
#' \eqn{(n x q x h+1)}-dimensional array of impulse responses
#' @export
#'
#' @examples
#' # Consider the following monetary policy example for illustration
#' X <- FRED_heavy$df
#' # positions of IP, CPI, FFR, and US/Swiss exchange rate in the data
#' int_ix <- c(6, 105, 77, 95)
#' # Estimate the static factor dimension r
#' (PC_r <- baingcriterion(X, rmax = 25)$PC)
#' (IC_r <- baingcriterion(X, rmax = 25)$IC)
#' # As seen from this example, PC criterion estimates r higher
#' # and suggests also the maximum value as an estimate of the static factor
#' # dimension. Choose a grid of values of r to compare results:
#' r <- c(4,8,16)
#' # We fix dynamic factor dimension q to 4 and the VAR order k to 2
#' est_irfs <- lapply(r, function(x) DfmRawImp(X = X, q = 4, r = x, k = 2, h = 50))
#' # Recursive identification of the model
#' for(j in 1:length(est_irfs)){
#'   # extract the lag zero polynomial corresponding to the variables of interest
#'   b0 <- est_irfs[[j]] %>% unclass %>% .[int_ix,,1]
#'   # cholesky decomposition of the covariance matrix
#'   C <- b0 %>% tcrossprod() %>% chol %>% t
#'   # recursive identification of the IRFs
#'   chol_irf <- irf_x(est_irfs[[j]], post_mat = solve(b0)%*%C)
#'   # normalize to 50bp shock to FFR on impact
#'   chol_irf <- chol_irf/chol_irf[int_ix[3],3,1]*0.5
#'   # cumulate IRFS for IP, inflation, and Exchange rate
#'   chol_irf[int_ix[c(1,2,4)], ,] <- chol_irf[int_ix[c(1,2,4)], ,] %>%
#'     apply(MARGIN = c(1,2), FUN = cumsum) %>%
#'     aperm(perm = c(2,3,1))
#'   # produce IRF figures
#'   est_irfs[[j]] <- plot_irfs(t(chol_irf[int_ix, 3,]), c("Ind_prod", "Inflation", "FF_rate", "ExRate"))
#' }
#' # plot the IRFs into a single figure
#' ## Not run:
#' gridExtra::marrangeGrob(c(est_irfs[[1]], est_irfs[[2]], est_irfs[[3]]),
#'                         ncol = 3, nrow = 4)
#' ## End(Not run)
DfmRawImp <- function(X, q, r, k, h){
  # save standard deviations for later and standardize data
  X <- as.matrix(X)
  WW <- apply(X, 2, sd)
  mu <- matrix(1, nrow = nrow(X), ncol = 1)%*%t(colMeans(X))
  x <- (X-mu)%*%diag(WW^(-1))
  # step 1: PCA
  Gamma0 = stats::cov(x)
  eig_gamma <- eigen(Gamma0)
  W <- eig_gamma$vectors[,1:r]
  # static factors, common component and idiosyncratic noise
  stat_fac <- x%*%W
  Chi <- stat_fac %*% t(W) %*% diag(WW) + mu
  xi <- X - Chi
  sigma_noise <- var(xi)
  # step 2: VAR on static factors
  var_fac <- est_ar_ols(stat_fac, p.max = k, p.min = k, mean_estimate = "intercept")
  var_sq <- lmfd(a = polm(dbind(3, diag(r), -var_fac$a)), b = diag(r))
  # step 3: low-rank approximation of the error term
  eig_eps <- eigen(cov(var_fac$residuals[-c(1:k),]))
  KM <- eig_eps$vectors[,1:q]%*%diag(sqrt(eig_eps$values[1:q]))
  # IRF of the DFM model
  var_irf <- pseries(var_sq, lag.max = h) %>% unclass
  dfm_irf <- irf_x(irf = var_irf, pre_mat = diag(WW) %*% W, post_mat = KM)
  return(dfm_irf)
}

#' @title Estimation of the D-DFM model for given data and model specification
#'
#' @description
#' \code{estim_wrap} estimates the model for a given data matrix and
#' a vector of Kronecker indices. The user must also specify the
#' maximum number of iterations used for the estimation as well as
#' the convergence criterion.
#'
#' @export
#'
#' @param df \eqn{T x n} data matrix, where \eqn{T} and \eqn{n}
#'  are the time series and cross-sectional dimensions, respectively
#' @param nu \eqn{q}-dimensional vector of Kronecker indices
#' @param degs Optional, vector of length 2, defining the lag orders of \eqn{c(z)} and \eqn{d(z)}
#' @param maxit Maximum number of iterations in the estimation procedure, defaults to 250
#' @param conv_crit Convergence criterion in the estimation, defaults to 1e-3
#' @param init0 For the estimation of initial values, specify either
#' 1) a positive integer, with \code{init0 = 1}: "no bootstrap"; or
#' 2) a positive integer, with \code{init0 > 1}: "number of bootstrap draws"; or
#' 3) a list of elements containing pre-specified initial values of
#' the SSM and \eqn{(q x q)} covariance matrix of error term in state equation
#' @param verbose Should the estimation process be printed?
#' @param h estimate IRF \eqn{h}-period ahead
#'
#' @return a list of elements
#' \item{npar}{Number of estimated parameters in \eqn{c(z)} and \eqn{d(z)}}
#' \item{ll_val}{The model log-likelihood value}
#' \item{aic}{Akaike IC for the model}
#' \item{bic}{Bayesian IC for the model}
#' \item{hqic}{Hannan-Quinn IC for the model}
#' \item{irf}{Impulse-responses in an array corresponding to the model}
#'
#' @examples
#' df <- FRED_light$df
#' # small Kronecker index dimension, should converge fast
#' nu <- c(1,1)
#' est_obj <- estim_wrap(df = df, nu = nu)
#' # graphical analysis of the convergence of the EM algorithm
#' \dontrun{
#' opar <- par() # save default plot setup
#' d <- data.frame("iter" = 1:length(est_obj$conv_stat$ll_value),
#'                 "loglik" = est_obj$conv_stat$ll_value,
#'                 "conv.log10" = est_obj$conv_stat$conv_cr %>% log10)
#' par(mar = c(5,5,2,5))
#' with(d, plot(iter, conv.log10, type="l", col="red3",
#'              ylab=expression(log[10](italic(Delta))),
#'              ylim=c(-3,0)))
#' par(new = T)
#' with(d, plot(iter, loglik, pch=16, axes=F, xlab=NA, ylab=NA, cex=1.2))
#' axis(side = 4)
#' mtext(side = 4, line = 3, 'model log likelihood value')
#' legend("topleft",
#'        legend=c(expression(log[10](italic(Delta))), "log-likelihood"),
#'        lty=c(1,0), pch=c(NA, 16), col=c("red3", "black"))
#' par(opar) # reset to default
#' }
estim_wrap <- function(df, nu, degs = NULL,
                       maxit = 250, conv_crit = 1e-3,
                       init0 = 1,
                       verbose = TRUE, h = 50){

  df <- as.matrix(df)
  if(is.null(degs)) degs <- rep(max(nu), 2)
  if(max(degs)!=max(nu)) stop("Polynomial degrees not compatible with the Kronecker indices")
  # scale data and save the standard deviations s.t. they can be multiplied to final IRFs
  sd_vec <- apply(df, 2, sd)
  mu <- outer(rep(1, nrow(df)), colMeans(df))
  df <- (df-mu)%*%diag(sd_vec^-1)
  # create the model template
  tmpl_ss <- tmpl_rmfd_echelon_ss(dim_out = ncol(df),
                                  nu = nu,
                                  degs = degs)
  # save this large matrix for the extraction of deep prms
  # such that no need to calculate it in every instance
  tmpl_ss$tmpl$xpx <- crossprod(tmpl_ss$tmpl$H)
  # initial values
  if(is.list(init0)){
    # this is for supplied initial values
    params0$params0 <- init0$params0
    params0$sigma <- init0$sigma
  } else {
    params0 <- boot_init(data = df,
                         nu = nu,
                         degs = degs,
                         tmpl_mod = tmpl_ss$tmpl,
                         nrep = init0)
  }

  # estimation
  rmfd_em <- update_em2(params_deep_init = params0$params0,
                        sigma_init = params0$sigma,
                        data_wide = t(df),
                        tmpl_ss = tmpl_ss,
                        MAXIT = maxit,
                        VERBOSE = verbose,
                        conv_crit = conv_crit)

  # calculate model selection criteria
  rmfd_em$npar <- length(params0$params0)-1 # consider only system parameters
  pen <- rmfd_em$npar/nrow(df)
  rmfd_em$aic <- -2*rmfd_em$ll_val + 2*pen
  rmfd_em$bic <- -2*rmfd_em$ll_val + log(nrow(df))*pen
  rmfd_em$hqic <- -2*rmfd_em$ll_val + 2*log(log(nrow(df)))*pen

  # transform the resulting stsp system into RMFD model
  rmfd_em$rmfd_final <- stsp2rmfd(stspsys =  rmfd_em$ssm_final$sys,
                                  dim_out = dim(df)[2],
                                  nu = nu,
                                  degs = degs)
  # raw IRFs
  rmfd_em$irf <- rmfd_em$rmfd_final %>% pseries(lag.max = h) %>% unclass
  rmfd_em$irf <- irf_x(rmfd_em$irf, pre_mat = diag(sd_vec))

  return(rmfd_em)
}

#' @title Obtain a bootstrap distribution of the IRFs
#'
#' @param df \code{(T x N)}-dimensional data matrix
#' @param nu \code{q}-dimensional vector of Kronecker indices
#' @param degs optional, vector of length 2, defining the orders of c(z) and d(z)
#' @param maxit Maximum number of iterations in the estimation procedure, defaults to 250
#' @param conv_crit Convergence criterion in the estimation, defaults to 1e-3
#' @param verbose Should the estimation process be printed?
#' @param h the IRF horizon \eqn{h}
#' @param nrep number of bootstrap draws
#' @param int_vars \code{q}-dimensional vector coding the positions of the variables of interest
#' @param save_file save the IRF array at the end?
#'
#' @return
#' An array of dimension \eqn{n x q h+1 x nrep} corresponding to the recursively identified
#' structural impulse-response function
#'
#' @keywords internal
#'
boot_irfs <- function(df, nu,
                      degs = NULL, maxit = 250, conv_crit = 1e-3,
                      verbose, h, nrep,
                      int_vars, save_file = FALSE){

  mc_ix <- 1
  dim_out <- ncol(df)
  dim_in <- length(nu)
  if(is.null(degs)) degs <- rep(max(nu), 2)
  rmfd_mc <- array(0, dim = c(dim_out, dim_in, h+1, nrep))
  while(mc_ix<=nrep){
    cat("Bootstrap draw no.", mc_ix, "\n")
    arg_list <- list(df = X_boot(df, L = 52),
                     nu = nu,
                     degs = degs,
                     verbose = verbose,
                     init0 = 1,
                     conv_crit = conv_crit,
                     maxit = maxit)

    est_obj <- do.call(estim_wrap, args = arg_list) %>% try()

    if(!inherits(est_obj, 'try-error')){
      # check that Cholesky does not fail
      is_sigma_pd <- min(eigen(est_obj$Sigma[1:dim_in,1:dim_in])$value) > sqrt(.Machine$double.eps)
      if(is_sigma_pd){
        rmfd_mc[,,,mc_ix] <- est_obj %>% chol_ident(int_vars = int_vars,
                                                    est_type = "rmfd")
        mc_ix <- mc_ix + 1
      }
    }
  }
  if(save_file) saveRDS(object = out, file = "rmfd_mc.Rds")
  return(rmfd_mc)
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
#' in IRF and \eqn{q} is the number of structural shocks
#' @param var_names vector of variable names of the \eqn{q} impacted variables
#' @param plot_title Title of the figure
#' @param label_y Should the variable name be printed next to y-axis?
#'
#' @return list containing the figures
#' @examples x <- rnorm(10)
#'
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
  irf_arr <- tibble::as_tibble(dplyr::bind_cols(irf_arr, lag = 0:(h-1)))
  irf_plot <- vector(mode = "list", length = q)

  for(jj in 1:q){
    irf_arr %>%
      dplyr::select(contains(var_names[jj]), lag) %>%
      ggplot2::ggplot(ggplot2::aes_string(x = "lag",
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
#' irf_ex <- rnorm(4*4*4) %>% array(dim = c(4,4,4))
#' # assume the recursive identification via cholesky is given as
#' chol_ex <- rnorm(16) %>% matrix(4,4) %>% tcrossprod() %>% chol() %>% t()
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

#' @title Cholesky identification of the IRF
#'
#' @description
#' \code{chol_ident} carries out the recursive identification using the Cholesky decomposition
#' of the \eqn{q x q} innovation covariance matrix. The resulting IRF has then
#' a structural interpretation depending on the ordering of the variables. The user
#' must specify whether the estimation method of the DFM is the S-DFM (\code{\link{DfmRawImp}})
#' or D-DFM (\code{\link{estim_wrap}}).
#'
#' @param est_obj \code{\link{estim_wrap}} or \code{\link{DfmRawImp}}
#' object containing the raw IRF and the matrix
#' used to obtain the lower triangular Cholesky factor
#' @param int_vars vector denoting the variables of interest and the respective ordering
#' @param est_type specify the DFM estimation method
#'
#' @return An IRF array of the same size as the raw IRF used as an input
#'
#' @keywords internal
#'
chol_ident <- function(est_obj, int_vars, est_type = c("rmfd", "fglr")){

  if(!est_type %in% c("rmfd", "fglr")) stop("Specify the estimation of IRFs properly")
  dim_in <- length(int_vars)
  if(est_type=="rmfd"){
    d0 <- est_obj$rmfd_final$d %>% unclass %>% .[int_vars,,1]
    # make sure the covariance matrix is p.d. for cholesky
    sigma <- (est_obj$Sigma[1:dim_in,1:dim_in]+t(est_obj$Sigma[1:dim_in,1:dim_in]))/2
    chol_mat <- t(chol(sigma[1:dim_in,1:dim_in]))
    chol_irf <- irf_x(est_obj$irf, post_mat = solve(d0)%*%chol_mat)
  } else {
    b0 <- est_obj$dfm_irf %>% unclass %>% .[int_vars,,1]
    C <- b0 %>% tcrossprod() %>% chol %>% t
    chol_irf <- irf_x(est_obj$dfm_irf, post_mat = solve(b0)%*%C)
  }

  return(chol_irf)
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

#' @title Obtain a block bootstrap data sample
#'
#' @description
#' \code{X_boot} performs a simple block bootstrap for the data matrix.
#' This is a direct copy of FG code used in Forni and Gambetti (2010)
#' available at \url{https://ars.els-cdn.com/content/image/1-s2.0-S0304393209001597-mmc1.zip}.
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
  blks <- ceiling(runif(20)*K)
  X_boot <- matrix(0, nrow = K*L, ncol = nvar)
  for(i in 1:K){
    X_boot[((i-1)*L + 1):(i*L),] <- X[((blks[i]-1)*L + 1):(blks[i]*L),]
  }
  return(X_boot)
}
