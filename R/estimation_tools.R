#' @title Two-step estimation method of S-DFM using VAR
#'
#' @description
#' \code{DfmRawImp} estimates non-structural IRFs for S-DFM according to section 2.3.1 in
#' \emph{Estimation of Impulse-Response Functions with Dynamic Factor Models:
#' A New Parametrization} available at \url{https://arxiv.org/pdf/2202.00310.pdf}.
#' The code used here follows closely to that of given in the supplementary
#' material to the Forni and Gambetti (2010), available at
#' \url{http://pareto.uab.es/lgambetti/ReplicaForniGambettiJME.zip}.
#'
#' @param X \eqn{T x n} data matrix, where \eqn{T} and \eqn{n}
#'  are the time series and cross-sectional dimensions, respectively
#' @param q The dynamic factor dimension
#' @param r Number of static factors, must be larger than or equal to \eqn{q}
#' @param k The VAR order for static factors
#' @param h Estimate IRF \eqn{h}-period ahead
#'
#' @return
#' \eqn{(n x q x h+1)}-dimensional array of impulse responses
#' @export
#'
#' @seealso Forni, M., & Gambetti, L. (2010). The dynamic effects of monetary policy:
#' A structural factor model approach. Journal of Monetary Economics, 57(2), 203-216.
#'
#' @examples
#' # Consider the following monetary policy example for illustration
#' # taken from Forni and Gambetti (2010)
#' X <- FG_data$df
#' # positions of IP, CPI, FFR, and US/Swiss exchange rate in the data
#' int_ix <- c(5, 96, 75, 106)
#' # Estimate the static factor dimension r
#' (PC_r <- baingcriterion(X, rmax = 25)$PC)
#' (IC_r <- baingcriterion(X, rmax = 25)$IC)
#' # As seen from this example, the criteria are not in consensus about the
#' # static factor dimension and suggest also the maximum value as an
#' # estimate for r dimension. Choose a grid of values of r to compare results:
#' r <- c(4,10,16)
#' # We fix dynamic factor dimension q to 4 and the VAR order k to 2
#' est_irfs <- lapply(r, function(x) DfmRawImp(X = X, q = 4, r = x, k = 2, h = 50))
#' # Recursive identification of the model
#' for(j in 1:length(est_irfs)){
#'   # extract the lag zero polynomial corresponding to the variables of interest
#'   b0 <- unclass(est_irfs[[j]])[int_ix, , 1]
#'   # cholesky decomposition of the residual covariance matrix
#'   C <- b0 %>% tcrossprod() %>% chol() %>% t()
#'   # recursive identification of the IRFs
#'   chol_irf <- irf_x(est_irfs[[j]], post_mat = solve(b0)%*%C)
#'   # normalize to 50bp shock to FFR on impact
#'   chol_irf <- chol_irf/chol_irf[int_ix[3],3,1]*0.5
#'   # cumulate IRFS for IP and inflation
#'   chol_irf[int_ix[c(1,2)], ,] <- chol_irf[int_ix[c(1,2)], ,] %>%
#'     apply(MARGIN = c(1,2), FUN = cumsum) %>%
#'     aperm(perm = c(2,3,1))
#'   # produce IRF figures
#'   est_irfs[[j]] <- plot_irfs(t(chol_irf[int_ix, 3,]),
#'   c("Ind_prod", "Inflation", "FF_rate", "ExRate"),
#'   plot_title = paste("DFM with ", r[j], " static factors", sep=""))
#' }
#' # plot the IRFs into a single figure
#' \dontrun{
#' gridExtra::marrangeGrob(c(est_irfs[[1]], est_irfs[[2]], est_irfs[[3]]),
#'                         ncol = 3, nrow = 4, top = NULL)}
DfmRawImp <- function(X, q, r, k, h){
  # save standard deviations for later and standardize data
  X <- as.matrix(X)
  WW <- apply(X, 2, stats::sd)
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
  sigma_noise <- stats::var(xi)
  # step 2: VAR on static factors
  var_fac <- est_ar_ols(stat_fac, p.max = k, p.min = k, mean_estimate = "intercept")
  var_sq <- lmfd(a = polm(dbind(3, diag(r), -var_fac$a)), b = diag(r))
  # step 3: low-rank approximation of the error term
  eig_eps <- eigen(stats::cov(var_fac$residuals[-c(1:k),]))
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
#' @param h The IRF horizon
#'
#' @return a list of elements
#' \item{Sigma}{The \eqn{(q x q)} reduced form residual covariance matrix}
#' \item{ll_val}{The log-likelihood value of the model}
#' \item{conv_stat}{Convergence criterion and the model log-likelihood value
#' for every iteration of the EM algorithm}
#' \item{npar}{The number of estimated parameters in \eqn{c(z)} and \eqn{d(z)}}
#' \item{aic}{Akaike IC for the model}
#' \item{bic}{Bayesian IC for the model}
#' \item{hqic}{Hannan-Quinn IC for the model}
#' \item{rmfd_final}{The final model in \code{rmfd} format}
#' \item{irf}{The estimated raw impulse response function in an array}
#'
#' @examples
#' df <- FRED_light$df
#' # small Kronecker index dimension, should converge fast
#' nu <- c(1,1)
#' est_obj <- estim_wrap(df = df, nu = nu)
#' # graphical analysis of the convergence of the EM algorithm
#' \dontrun{
#' # create data frame containing the number of iterations, log-likelihood
#' # value of the model for given iteration and the convergence
#' # criterion for plotting purposes
#' d <- data.frame("iter" = 1:length(est_obj$conv_stat$ll_value),
#'                 "loglik" = est_obj$conv_stat$ll_value,
#'                 "conv.log10" = est_obj$conv_stat$conv_cr %>% log10)
#' opar <- par()$mar # save default plot setup
#' par(mar = c(5,5,2,5))
#' with(d, plot(iter, conv.log10, type="l", col="red3",
#'              ylab=expression(log[10](italic(Delta))),
#'              ylim=c(-3,0)))
#' par(new = T)
#' with(d, plot(iter, loglik, pch=1, axes=F, xlab=NA, ylab=NA, cex=1.2))
#' axis(side = 4)
#' mtext(side = 4, line = 3, 'model log likelihood value')
#' legend("topleft",
#'        legend=c(expression(log[10](italic(Delta))), "log-likelihood"),
#'        lty=c(1,0), pch=c(NA, 1), col=c("red3", "black"))
#' # reset graphical parameters to default
#' par(new = FALSE, mar = opar)}
estim_wrap <- function(df, nu, degs = NULL,
                       maxit = 250, conv_crit = 1e-3,
                       init0 = 1,
                       verbose = TRUE, h = 50){

  df <- as.matrix(df)
  nobs <- dim(df)[1] # time series dim
  nvar <- dim(df)[2] # cross-section dim
  if(is.null(degs)) degs <- rep(max(nu), 2)
  if(max(degs)!=max(nu)) stop("Polynomial degrees not compatible with the Kronecker indices")
  # scale data and save the standard deviations s.t. they can be multiplied to final IRFs
  sd_vec <- apply(df, 2, stats::sd)
  mu <- outer(rep(1, nobs), colMeans(df))
  df <- (df-mu)%*%diag(sd_vec^-1)
  # create the model template
  tmpl_ss <- tmpl_rmfd_echelon_ss(dim_out = nvar,
                                  nu = nu,
                                  degs = degs)
  # save this large matrix for the extraction of deep prms
  # such that no need to calculate it in every iteration
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
                         tmpl_mod = tmpl_ss,
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
  pen <- rmfd_em$npar/nobs
  rmfd_em$aic <- -2*rmfd_em$ll_val + 2*pen
  rmfd_em$bic <- -2*rmfd_em$ll_val + log(nobs)*pen
  rmfd_em$hqic <- -2*rmfd_em$ll_val + 2*log(log(nobs))*pen

  # transform the resulting stsp system into RMFD model
  rmfd_em$rmfd_final <- stsp2rmfd(stspsys =  rmfd_em$ssm_final$sys,
                                  dim_out = nvar,
                                  nu = nu,
                                  degs = degs)
  rmfd_em <- rmfd_em[-1] # drop the redundant element from the output list

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

    est_obj <- do.call(estim_wrap, args = arg_list) %>% try(silent = TRUE)

    if(!inherits(est_obj, 'try-error')){
      # check that Cholesky does not fail
      is_sigma_pd <- min(eigen(est_obj$Sigma[1:dim_in,1:dim_in])$value) > sqrt(.Machine$double.eps)
      if(is_sigma_pd){
        rmfd_mc[,,,mc_ix] <- est_obj %>% chol_ident(int_vars = int_vars,
                                                    est_type = "rmfd")
        mc_ix <- mc_ix + 1
      }
    } else{
      cat("\nFailure to converge, try again...\n")
    }
  }
  if(save_file) saveRDS(object = rmfd_mc, file = "rmfd_mc.Rds")
  return(rmfd_mc)
}

#' @title Cholesky identification of the IRF
#'
#' @description
#' \code{chol_ident} carries out the recursive identification using the Cholesky decomposition
#' of the \eqn{(q x q)} innovation covariance matrix. The resulting IRF has then
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
