#' @title Estimate IRFs and associated bootstrap intervals for D-DFM
#'
#' @export
#'
#' @param df list containing elements:
#' 1) \code{df}, data in a \code{T x n} matrix;
#' 2) \code{int_ix}, index coding the column positions of the variables of interest in the data;
#' 3) \code{trans_ix}, index giving the variable transformation codes as in McCracken and Ng (2016)
#' @param r static factor dimension
#' @param h h-step ahead IRFs
#' @param nrep number of replications in the block bootstrap procedure
#' @param conv_crit convergence criterion of the EM algorithm
#' @param ci confidence intervals for the IRFs
#' @param init_rep replications for the starting value procedure, for details, see \code{\link{boot_init}}
#' @param verbose logical, print the estimation process?
#' @param mod_str optional list for a fixed model structure, containing at least the Kronecker index vector
#' and also possibly the orders of c(z) and d(z) as a second element
#'
#' @return
#'
#' @seealso McCracken, M. W., & Ng, S. (2016). FRED-MD: A monthly database for macroeconomic research.
#' Journal of Business & Economic Statistics, 34(4), 574-589.
#'
#' @examples
do_everything_rmfd <- function(df, r, h = 50, nrep = 0,
                               conv_crit = 1e-5, ci = .8,
                               init_rep = 0, verbose = FALSE,
                               mod_str = NULL){
  df$df <- as.matrix(df$df)
  q <- length(df$int_ix)
  max_nu <- floor(r/q) # maximum Kronecker index

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
                     verbose = verbose,
                     init0 = init_rep,
                     conv_crit = conv_crit,
                     maxit = 500)
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
    rmfd_irf <- apply(rmfd_mc$irfs, c(1,2,3), function(x) quantile(x, probs = c((1-ci)/2, .5, 1-(1-ci)/2)))
    rmfd_irf <- sapply(1:3, function(x) rmfd_irf[x,,,] %>% finalize_irf(shock_size = 0.5,
                                                                        norm_id = c(df$int_ix[3],3),
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
                         finalize_irf(shock_size = 0.5,
                                      norm_id = c(df$int_ix[3],3),
                                      trans_ix = df$trans_ix,
                                      int_vars = df$int_ix),
                       simplify = "array") %>% aperm(c(1,3,2))
  }

  if(n_mods > 1){
    out <- list("irf" = rmfd_irf,
                "rirf" = est_list[[min_ix]]$irf,
                "model_selection" = mod_sel,
                "conv_stat" = est_list[[min_ix]]$conv_stat)
  } else{
    out <- list("irf" = rmfd_irf,
                "rirf" = est_list[[min_ix]]$irf,
                "conv_stat" = est_list[[min_ix]]$conv_stat)
  }
  return(out)
}

#' @title Estimate IRFs and associated bootstrap intervals for S-DFM
#'
#' @export
#'
#' @param df list containing items:
#' 1) \code{df}, data in a tall matrix;
#' 2) \code{int_ix}, index coding the rows of the interest variables in the data;
#' 3) \code{trans_ix}, index giving the variable transformation codes
#' @param r static factor dimension
#' @param k VAR degree on the static factors
#' @param h h-step ahead IRFs
#' @param nrep number of replications in the bootstrap
#' @param ci confidence intervals for the IRFs
#'
#' @return
#'
#' @examples
do_everything_fglr <- function(df, r, k, h, nrep = 500, ci = 0.8){

  q <- length(df$int_ix)
  # Confidence intervals
  fglr_mc <- replicate(nrep, DfmRawImp(X = X_boot(df$df, 52), q = q, r = r, k = k, h = h))
  fglr_irf <- sapply(1:nrep,
                     function(x) chol_ident(fglr_mc["dfm_irf", x],
                                            int_vars = df$int_ix,
                                            est_type = "fglr") %>%
                       finalize_irf(shock_size = 0.5,
                                    norm_id = c(df$int_ix[3],3),
                                    int_vars = df$int_ix,
                                    trans_ix = df$trans_ix),
                     simplify = "array")

  fglr_irf <- apply(fglr_irf, c(1,2),
                    function(x) stats::quantile(x, probs = c((1-ci)/2, .5, 1-(1-ci)/2))) %>%
    aperm(perm = c(2,1,3))

  # Add point estimates
  fglr_point <- DfmRawImp(X = df$df, q = q, r = r, k = k, h = h)
  fglr_irf[,2,] <- fglr_point %>%
    chol_ident(int_vars = df$int_ix, est_type = "fglr") %>%
    finalize_irf(shock_size = 0.5,
                 norm_id = c(df$int_ix[3],3),
                 int_vars = df$int_ix,
                 trans_ix = df$trans_ix)

  return(list("irf" = fglr_irf))
}

#' Estimate IRFs and associated bootstrap intervals for SVAR
#'
#' @export
#'
#' @param df list containing items:
#' 1) \code{df}, data in a tall matrix;
#' 2) \code{int_ix}, index coding the rows of the interest variables in the data;
#' 3) \code{trans_ix}, index giving the variable transformation codes
#' @param k VAR order
#' @param h h-step ahead IRFs
#' @param nrep number of replications in the bootstrap
#' @param ci confidence intervals for the IRFs
#'
#' @return
#'
#' @examples
do_everything_svar <- function(df, k, nrep = 500, h = 50, ci = 0.8){

  q <- length(df$int_ix)
  irf_arr <- array(0, dim = c(q, q, h+1, nrep))
  for(j in 1:500){
    var_boot <- est_ar_ols(y = X_boot(df$df[,df$int_ix], 52),
                           p.max = k,
                           p.min = k,
                           mean_estimate = 'intercept')
    irf_arr[,,,j] <- var_boot$model$sys %>%
      pseries(lag.max = h) %>%
      unclass() %>%
      irf_x(post_mat = var_boot$model$sigma_L)
  }

  svar_irf <- apply(irf_arr, c(1,2,3), function(x) quantile(x, probs = c((1-ci)/2, .5, 1-(1-ci)/2)))
  svar_irf <- sapply(1:3, function(x) svar_irf[x,,,] %>% finalize_irf(shock_size = 0.5,
                                                                      norm_id = c(3,3),
                                                                      trans_ix = df$trans_ix[df$int_ix],
                                                                      int_vars = 1:q),
                     simplify = "array") %>% aperm(c(1,3,2))
  # add point estimates to the IRF array
  var_point <- est_ar_ols(y = df$df[,df$int_ix],
                          p.max = k,
                          p.min = k,
                          mean_estimate = 'intercept')
  svar_irf[,2,] <- var_point$model$sys %>%
    pseries(lag.max = h) %>%
    unclass() %>%
    irf_x(post_mat = var_point$model$sigma_L) %>%
    finalize_irf(shock_size = 0.5,
                 norm_id = c(3,3),
                 trans_ix = df$trans_ix[df$int_ix],
                 int_vars = 1:q)
  return(svar_irf)
}
