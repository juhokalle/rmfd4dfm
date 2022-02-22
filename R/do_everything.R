#' @title Estimate IRFs and associated bootstrap intervals for D-DFM
#'
#' @description
#' \code{do_everything_rmfd} is a user-friendly function for implementing the estimation
#' and identification methodology of the D-DFM model presented in sections 2.2.
#' and 3, as well as replicating the empirical exercise of Section 5 in \emph{Estimation
#' of Impulse-Response Functions with Dynamic Factor Models: A New Parametrization},
#' available at \url{https://arxiv.org/pdf/2202.00310.pdf}.
#'
#' @export
#'
#' @param df a list with the following elements:
#' 1) \code{df}, data in a \eqn{T x n} matrix;
#' 2) \code{int_ix}, an index vector coding the columns corresponding to the variables of interest in the data;
#' 3) \code{trans_ix}, an index vector giving the variable transformation codes as in McCracken and Ng (2016)
#' 4) \code{shock_ix}, specify the shock of interest and the normalization at impact, defaults to \code{c(3,0.5)}
#' @param r state dimension of the D-DFM in state space model
#' @param h estimate \eqn{h}-step ahead IRFs
#' @param nrep number of replications in the block bootstrap procedure, set to zero for no bootstrapping
#' @param conv_crit convergence criterion of the EM algorithm, the default is 1e-5
#' @param ci confidence interval level for the IRFs, defaults to 0.68
#' @param init_rep bootstrap replications for the initial value procedure,
#' for details, see \code{\link{boot_init}}
#' @param verbose logical, print the estimation process?
#' @param mod_str optional list for a fixed model structure, containing at least the Kronecker index vector
#' and also possibly the orders of \eqn{c(z)} and \eqn{d(z)} as a second element
#'
#' @details
#' The structure of the function is as follows:
#' 1) Given the state dimension of the D-DFM in state space format, which can be estimated
#' using functions \code{\link{abc_crit}} and \code{\link{baingcriterion}}, for example,
#' the function returns the set of model structures consistent with the model selection
#' criteria 3--5 given in Section 3.3 of the paper via the function \code{\link{admissible_mods}}.
#' Alternatively, the user can estimate a fixed model structure using the argument \code{mod_str} without
#' running the estimation algorithm across a set of model alternatives.
#' 2) The function estimates a set of feasible model(s) using the function \code{\link{estim_wrap}}
#' and returns the statistics of the model selection criteria. The best model candidate is chosen
#' according to the BIC statistic.
#' 3) The function estimates the bootstrapped confidence intervals for IRFs of the recursively identified
#' model with the smallest BIC value. In the case of no bootstrapping (\code{nrep=0}), the function returns
#' all the estimated IRFs, while with bootstrapping the function returns only the point estimates of the
#' IRF corresponding to the best model candidate determined using BIC.
#'
#' For obtaining the qualified models and an estimation example, we refer the reader
#' to \code{\link{admissible_mods}} and \code{\link{estim_wrap}}. The running time of the algorithm
#' is determined by the size of the data matrix, complexity of the estimated model,
#' level of the convergence criterion, and the number of bootstrap draws.
#' For a standard data set (\eqn{n~100}, \eqn{T~300}), simple model structure,
#' a lenient convergence criterion (\code{conv_crit=1e-3}) and a small number of
#' bootstrap draws (\code{nrep=100}), the function should be finished in approximately 60 minutes.
#'
#' The estimated models are identified recursively using Cholesky decomposition of the residual covariance
#' matrix. As a default, the function returns the impulse responses to the third shock, the size of which
#' is normalized to 0.5 on impact. For changing the shock of interest (i.e. not the third),
#' and the normalization constant, the user should include store a vector of length 2 as \code{df$shock_ix}
#' with the position of the shock of interest as the first element and the
#' normalization constant as the second element. For a more flexible setup, such as identifying
#' more than one shock or obtaining responses to many variables, the user can use function
#' \code{\link{estim_wrap}}, which returns the non-structural IRF.
#'
#' @return a list with elements
#' \item{irf}{an array corresponding to the recursively identified IRF}
#' \item{conv_crit}{convergence criteria for the model minimizing BIC}
#' \item{model_selection}{table containing the model selection criteria, returned only if
#' the number of estimated models is larger than one}
#'
#' @seealso McCracken, M. W., & Ng, S. (2016). FRED-MD: A monthly database for macroeconomic research.
#' Journal of Business & Economic Statistics, 34(4), 574-589.
#'
#' @examples
#' \dontrun{
#' # the following example takes around 60 minutes to finish
#' FRED_heavy$int_ix <- c(5,96,77,101) # remember to index the variables of interest
#' est_obj <- do_everything_rmfd(FRED_heavy, r = 8, h = 48, nrep = 100,
#'   conv_crit = 1e-3, ci = .68, init_rep = 0, verbose = TRUE)}
#'
do_everything_rmfd <- function(df, r, h, nrep,
                               conv_crit = 1e-5, ci = .68,
                               init_rep = 0, verbose = FALSE,
                               mod_str = NULL){
  df$df <- as.matrix(df$df)
  q <- length(df$int_ix)
  max_nu <- floor(r/q) # maximum Kronecker index
  # fix the shock of interest and the normalization
  if(is.null(df$shock_ix)) df$shock_ix <- c(3, 0.5)
  # Model estimation ####
  # get admissible models and fix the state dimension if max(Kronecker index) > 1
  if(is.null(mod_str)){
    int_list <- admissible_mods(q = q, max_nu = max_nu, fix_r = max_nu>1)
  } else if(is.list(mod_str)){
    if(length(mod_str)==1) mod_str[[2]] <- rep(max(mod_str[[1]]), 2)
    int_list <- list(nus = list(mod_str[[1]]), degs = list(mod_str[[2]]))
  } else if(!is.list(int_list)){
    stop("Specify the model correctly.")
  }
  n_mods <- length(int_list[[1]]) # number of candidate models
  est_list <- vector(mode = "list", length = n_mods) # preallocate result list

  for(jj in 1:n_mods){
    cat("\nmodel no. ", jj, "/", n_mods, " idx: ", int_list$nus[[jj]],
        ", degs: (", int_list$degs[[jj]], ")", "\n", sep = "")
    arg_list <- list(df = df$df,
                     nu = int_list$nus[[jj]],
                     degs = int_list$degs[[jj]],
                     maxit = 500,
                     conv_crit = conv_crit,
                     init0 = init_rep,
                     verbose = verbose,
                     h=h)
    est_list[[jj]] <- do.call(estim_wrap, args = arg_list) %>% try()
    while(inherits(est_list[[jj]], 'try-error')){
      est_list[[jj]] <- do.call(estim_wrap, args = arg_list) %>% try()
    }
  }

  # Model selection ####
  if(n_mods>1){
    mod_sel <- sapply(c("ll_val", "aic", "bic", "hqic", "npar"),
                      function(z) sapply(est_list, function(x) x[[z]]))
    rownames(mod_sel) <- lapply(int_list$nus, function(x) paste("(", paste(x, collapse = ","), ")", sep=""))
    mod_sel <- cbind(round(mod_sel, 2),
                     "p,s" = sapply(1:n_mods, function(x) paste(int_list$degs[[x]], collapse = "," )))
  }
  min_ix <- ifelse(n_mods>1, which.min(sapply(est_list, function(x) x$bic)), 1)

  # Bootstrap confidence intervals ####
  if(nrep>0){
    rmfd_mc <- boot_irfs(df = df$df,
                         nu = int_list$nus[[min_ix]],
                         degs = int_list$degs[[min_ix]],
                         maxit = 250,
                         conv_crit = conv_crit,
                         h = h,
                         verbose = verbose,
                         nrep = nrep,
                         int_vars = df$int_ix,
                         save_file = TRUE)
    rmfd_irf <- apply(rmfd_mc$irfs, c(1,2,3), function(x) stats::quantile(x, probs = c((1-ci)/2, .5, 1-(1-ci)/2)))
    rmfd_irf <- sapply(1:3, function(x) rmfd_irf[x,,,] %>% finalize_irf(shock_size = df$shock_ix[2],
                                                                        norm_id = c(df$int_ix[df$shock_ix[1]],
                                                                                    df$shock_ix[1]),
                                                                        trans_ix = df$trans_ix,
                                                                        int_vars = df$int_ix),
                       simplify = "array") %>% aperm(c(1,3,2))

    # add point estimates to the IRF array
    rmfd_irf[, 2, ] <- est_list[[min_ix]] %>%
      chol_ident(int_vars = df$int_ix, est_type = "rmfd") %>%
      finalize_irf(shock_size = 0.5,
                   norm_id = c(df$int_ix[3],3),
                   trans_ix = df$trans_ix,
                   int_vars = df$int_ix)
  } else{
    # return all irfs if no bootstrapping
    rmfd_irf <- sapply(1:n_mods,
                       function(x) est_list[[x]] %>%
                         chol_ident(int_vars = df$int_ix, est_type = "rmfd") %>%
                         finalize_irf(shock_size = df$shock_ix[2],
                                      norm_id = c(df$int_ix[df$shock_ix[1]],
                                                  df$shock_ix[1]),
                                      trans_ix = df$trans_ix,
                                      int_vars = df$int_ix),
                       simplify = "array") %>% aperm(c(1,3,2))
  }

  if(n_mods > 1){
    out <- list("irf" = rmfd_irf,
                "conv_stat" = est_list[[min_ix]]$conv_stat,
                "model_selection" = mod_sel)
  } else{
    out <- list("irf" = rmfd_irf,
                "conv_stat" = est_list[[min_ix]]$conv_stat)
  }
  return(out)
}

#' @title Estimate IRFs and associated bootstrap intervals for S-DFM
#'
#' @description
#' \code{do_everything_fglr} is a user-friendly wrapper estimating
#' the recursively identified S-DFM model and associated bootstrap
#' intervals for IRFs within one function. The function allows the user
#' to specify the model similarly to \code{\link{DfmRawImp}}, to which
#' the user is referred for examples.
#'
#' @export
#'
#' @param df list containing items:
#' 1) \code{df}, data in a \eqn{T x n} matrix;
#' 2) \code{int_ix}, an index vector coding the columns corresponding to the variables of interest in the data;
#' 3) \code{trans_ix}, an index vector giving the variable transformation codes as in McCracken and Ng (2016);
#' 4) \code{shock_ix}, specify the shock of interest and the normalization at impact, defaults to \code{c(3,0.5)}
#' @param r static factor dimension
#' @param k VAR degree on the static factors
#' @param h h-step ahead IRFs
#' @param nrep number of replications in the bootstrap
#' @param ci confidence intervals for the IRFs
#'
#' The estimated models are identified recursively using Cholesky decomposition of the residual covariance
#' matrix. As a default, the function returns the impulse responses to the third shock, the size of which
#' is normalized to 0.5 on impact. For changing the shock of interest (i.e. not the third),
#' and the normalization constant, the user should include store a vector of length 2 as \code{df$shock_ix}
#' with the position of the shock of interest as the first element and the
#' normalization constant as the second element. For a more flexible setup, such as identifying
#' more than one shock or obtaining responses to many variables, the user can use function
#' \code{\link{DfmRawImp}}, which returns the non-structural IRF.
#'
#' @seealso McCracken, M. W., & Ng, S. (2016). FRED-MD: A monthly database for macroeconomic research.
#' Journal of Business & Economic Statistics, 34(4), 574-589.
#'
#' @return Structural impulse response function as an array
#'
do_everything_fglr <- function(df, r, k, h, nrep = 500, ci = 0.8){

  q <- length(df$int_ix)
  # fix the shock of interest and the normalization
  if(is.null(df$shock_ix)) df$shock_ix <- c(3, 0.5)
  # Confidence intervals
  fglr_mc <- replicate(nrep, DfmRawImp(X = X_boot(df$df, 52), q = q, r = r, k = k, h = h))
  fglr_irf <- sapply(1:nrep,
                     function(x) chol_ident(fglr_mc[, , ,x],
                                            int_vars = df$int_ix,
                                            est_type = "fglr") %>%
                       finalize_irf(shock_size = df$shock_ix[2],
                                    norm_id = c(df$int_ix[df$shock_ix[1]],
                                                df$shock_ix[1]),
                                    trans_ix = df$trans_ix,
                                    int_vars = df$int_ix),
                     simplify = "array")
  fglr_irf <- apply(fglr_irf, c(1,2),
                    function(x) stats::quantile(x, probs = c((1-ci)/2, .5, 1-(1-ci)/2))) %>%
    aperm(perm = c(2,1,3))

  # Add point estimates
  fglr_point <- DfmRawImp(X = df$df, q = q, r = r, k = k, h = h)
  fglr_irf[,2,] <- fglr_point %>%
    chol_ident(int_vars = df$int_ix, est_type = "fglr") %>%
    finalize_irf(shock_size = df$shock_ix[2],
                 norm_id = c(df$int_ix[df$shock_ix[1]],
                             df$shock_ix[1]),
                 trans_ix = df$trans_ix,
                 int_vars = df$int_ix)
  return(fglr_irf)
}

#' @title Estimate IRFs and associated bootstrap intervals for SVAR
#'
#' @description
#' Similarly to other \code{do_everything_*} functions, \code{do_everything_svar}
#' estimates the benchmark SVAR model and the associated IRFs with error bands.
#' Note that the data object is similar as in the other \code{do_everything_*} functions,
#' i.e. the data matrix contains all the variables and those included in the
#' SVAR are indexed with a vector (see below).
#'
#' @export
#'
#' @param df list containing items:
#' 1) \code{df}, data in a \eqn{T x n} matrix;
#' 2) \code{int_ix}, an index vector coding the columns corresponding to the variables included in the SVAR;
#' 3) \code{trans_ix}, an index vector giving the variable transformation codes as in McCracken and Ng (2016);
#' 4) \code{shock_ix}, specify the shock of interest and the normalization at impact, defaults to \code{c(3,0.5)}
#' @param k VAR order
#' @param h h-step ahead IRFs
#' @param nrep number of replications in the bootstrap
#' @param ci confidence interval level for the IRFs, defaults to 0.68
#'
#' The estimated models are identified recursively using Cholesky decomposition of the residual covariance
#' matrix. As a default, the function returns the impulse responses to the third shock, the size of which
#' is normalized to 0.5 on impact. For changing the shock of interest (i.e. not the third),
#' and the normalization constant, the user should include a vector of length 2 in the data object
#' with name "shock_ix" with the position of the shock of interest as the first and the
#' normalization constant as the second element.
#'
#' @references McCracken, M. W., & Ng, S. (2016). FRED-MD: A monthly database for macroeconomic research.
#' Journal of Business & Economic Statistics, 34(4), 574-589.
#'
#' @return An \eqn{(h+1 x 3 x q)} array corresponding to the structurally identified IRFs
#'
#' @examples
#' # example from the paper: the data object is similar to the
#' # other do_everything_* function, i.e. the data matrix
#' # is the whole data
#' FRED_light$int_ix <- c(6, 105, 77, 95)
#' est_obj <- do_everything_svar(FRED_light, k = 9, nrep = 500, h = 50, ci = .68)
do_everything_svar <- function(df, k, nrep = 500, h = 50, ci = 0.68){

  q <- length(df$int_ix)
  # fix the shock of interest and the normalization
  if(is.null(df$shock_ix)) df$shock_ix <- c(3, 0.5)
  # pre-allocate
  irf_arr <- array(0, dim = c(q, q, h+1, nrep))
  # obtain first the bootstrap sample and subsequently the point estimates
  for(j in 1:nrep){
    var_boot <- est_ar(obj = X_boot(df$df[,df$int_ix], 52),
                       p.max = k,
                       method = "ols",
                       ic = "max",
                       mean_estimate = 'intercept')
    irf_arr[, , ,j] <- var_boot$model$sys %>%
      pseries(lag.max = h) %>%
      unclass() %>%
      irf_x(post_mat = var_boot$model$sigma_L)
  }

  svar_irf <- apply(irf_arr, c(1,2,3), function(x) stats::quantile(x, probs = c((1-ci)/2, .5, 1-(1-ci)/2)))
  svar_irf <- sapply(1:3, function(x) svar_irf[x,,,] %>% finalize_irf(shock_size = df$shock_ix[2],
                                                                      norm_id = rep(df$shock_ix[1],2),
                                                                      trans_ix = df$trans_ix[df$int_ix],
                                                                      int_vars = 1:q),
                     simplify = "array") %>% aperm(c(1,3,2))
  # add point estimates to the IRF array
  var_point <- est_ar(obj = df$df[,df$int_ix],
                      p.max = k,
                      method = "ols",
                      ic = "max",
                      mean_estimate = 'intercept')
  svar_irf[,2,] <- var_point$model$sys %>%
    pseries(lag.max = h) %>%
    unclass() %>%
    irf_x(post_mat = var_point$model$sigma_L) %>%
    finalize_irf(shock_size = df$shock_ix[2],
                 norm_id = rep(df$shock_ix[1],2),
                 trans_ix = df$trans_ix[df$int_ix],
                 int_vars = 1:q)
  return(svar_irf)
}
