#' @title Get Tall IRFs
#'
#' @description
#' \code{lr_approx} calculates a tall IRF from a square
#' state space model using low-rank approximation
#'
#' @param ssmod square state space model object
#' @param dim_in the dynamic factor dimension of the DFM
#' @param lag.max number of lags, default is 30
#'
#' @return a list of elements:
#' \item{irf}{new IRF array}
#' \item{noise_L}{n x q matrix used to obtain tall IRF}
#'
#' @keywords internal
lr_approx <- function(ssmod, dim_in, lag.max = 30){

  # low rank approximation to obtain (n x dim_in) IRF corresponding to "ssmod"
  sigma <- tcrossprod(ssmod$sigma_L)
  sigma_evd <- sigma %>% eigen(symmetric=TRUE)
  sigma_L_red <- sigma_evd$vectors[, 1:dim_in, drop = FALSE] %*% diag(sqrt(sigma_evd$values[1:dim_in]))
  irf_new <- pseries(ssmod$sys, lag.max = lag.max) %>% unclass()
  out <- irf_x(irf_new, post_mat = sigma_L_red)
  return(list(irf = out, noise_L = sigma_L_red))
}

#' @title Normalize zero-lag coefficient of the IRF
#'
#' @description
#' Depending on the template/identification scheme, the first \code{dim_in} rows of
#' the zero-lag coefficient of the transfer function are either normalized to the
#' identity matrix or lower triangular.
#'
#' @param pseries_obj IRF object, should have sufficiently many lags
#' @param type Character string. Either \code{identity} or \code{lq_decomp}.
#'   In the former case the first submatrix of the tall transfer function will be transformed to the identity matrix.
#'   In the latter case, this submatrix will be made lower triangular by an orthogonal transformation.
#'
#' @return a list of elements
#' \item{pseries_array}{the normalised transfer function as array}
#' \item{noise_L}{the inverse of the matrix used for normalising}
#' @keywords internal
normalize_irf = function(pseries_obj, type = c("identity", "lq_decomp")){

  type = match.arg(type)

  pseries_array = pseries_obj %>% unclass
  dim_pseries = dim(pseries_array)
  dim_in = dim_pseries[2]

  stopifnot("normalize_irf(): Input argument *type* must be a string equal to either *identity* or *lq_decomp*" = length(type) == 1)

  if (type == "identity"){
    noise_L = pseries_array[1:dim_in, 1:dim_in, 1]
    post_mat = noise_L %>% solve()
  } else if (type == "lq_decomp"){
    tmp = qr(t(pseries_array[1:dim_in, 1:dim_in, 1]))
    noise_L = t(qr.Q(tmp))
    post_mat = noise_L %>% t()
  } else {
    stop("normalize_irf(): Something wrong with input argument *type*.")
  }

  pseries_array = irf_x(pseries_array, post_mat = post_mat)

  return(list(pseries_array = pseries_array,
              noise_L = noise_L))
}

#' @title Initial values for RMFD model
#'
#' @description
#' \code{init_vals} gives the initial values of the RMFD model in state space format.
#' The initial values are calculated from the input IRF that is used to generate the Hankel matrix,
#' which are used in the regressions to determine parameter matrices associated with the lag
#' polynomials c(z) and d(z).
#'
#' @param irf_obj impulse response function object of dimension (\eqn{n, q,} \code{lag.max}+1)
#' @param nu vector of  Kronecker indices
#' @param degs vector of length 2, where one can optionally specify the degrees of c(z) and d(z),
#' either of the parameters must be equal to \code{max(nu)}
#' @param tmpl_mod model template for extracting the initial values
#'
#' @return a list of elements
#' \item{init_val}{a vector of initial values}
#' \item{sigma}{error covariance matrix}
#'
#' @keywords internal
init_vals <- function(irf_obj, nu, degs = NULL, tmpl_mod){

  # extract integer valued prms
  if(is.null(degs)) degs <- rep(max(nu), 2)
  deg_c <- degs[1]
  deg_d <- degs[2]
  dim_out <- dim(irf_obj[[1]])[1]
  dim_in <- dim(irf_obj[[1]])[2]
  dim_state <- max(nu, deg_d+1)*dim_in
  # transform the input IRF series into RMFD
  irf_normalised <- irf_obj$pseries_array
  rmfd_sys <- pseries2rmfd(obj = irf_normalised, Hsize = max(nu+1), mu=nu)$Xr
  # extract and restrict the degree of c(z) and d(z) according to the prespecified values
  # the polynomials remain unchanged if no restrictions to lag orders are made
  rmfd_c <- unclass(rmfd_sys$c)[, , 1:max(nu+1) %in% 1:(deg_c+1)]
  rmfd_d <- unclass(rmfd_sys$d)[, , 1:max(nu+1) %in% 1:(deg_d+1)]
  rmfd_sys <- rmfd(c = rmfd_c, d = rmfd_d)
  ss_init <- rmfd2stsp(rmfd_sys, nu = nu)

  # obtain the vector of initial values
  prm_vec <- c(as.vector(ss_init), as.vector(diag(dim_out))) # set sigma_noise to I_n
  prms_out <- solve(tmpl_mod$xpx, crossprod(tmpl_mod$H, prm_vec - tmpl_mod$h))

  return(list(prms_out = prms_out,
              sigma = tcrossprod(irf_obj$noise_L)))
}

#' @title Obtain initial values corresponding to a given model specification
#'
#' @description
#' \code{init_wrap} estimates the initial values for a given model specification, this
#' function is a helper to the function \code{boot_init}
#'
#' @param data matrix of dimension \eqn{(T x N)}, where\eqn{T} is the number of observations and
#' \eqn{N} is the number of variables
#' @param nu vector of right Kronecker indices
#' @param degs vector of length 2, where one can optionally specify the degrees of c(z) and d(z), either of the parameters must be equal to max(nu)
#' @param tmpl_mod model template, see \code{rmfd4dfm:::tmpl_rmfd_echelon_ss}
#' @param reg_prm regularization parameter governing the std deviation of the randomly drawn noise vector,
#' see also \link[rmfd4dfm]{estim_wrap}
#'
#' @return a list of elements
#' \item{params0}{a vector of initial parameters}
#' \item{sigma}{\code{q x q} initial error covariance matrix}
#'
#' @keywords internal
init_wrap = function(data, nu, degs = NULL, tmpl_mod, reg_prm = NULL){

  n_obs <- dim(data)[1]
  dim_out <- dim(data)[2]
  dim_in <- length(nu)
  dim_state <- max(nu)*dim_in # minimal state dimension
  regu_mat <- matrix(0, nrow = n_obs, ncol = dim_out)
  iter <- 0
  while(iter == 0 || inherits(ssm_sq, 'try-error')){
    # fill regularization matrix if so specified
    if(!is.null(reg_prm)) regu_mat[] <- stats::rnorm(n_obs*dim_out, 0, 10^(iter-reg_prm))
    ssm_sq <- est_stsp_ss(obj = data + regu_mat, method = 'cca',
                          s.max = dim_state, mean_estimate = "zero",
                          estorder = estorder_max) %>% try(silent=TRUE)
    if(is.null(reg_prm)) reg_prm <- 8
    iter <- iter + 1
  }
  tall_irf <- lr_approx(ssm_sq$model, dim_in)$irf
  sigma_noise <- tcrossprod(ssm_sq$model$sigma_L)
  # normalize such that the top q rows of k_0 are an identity matrix
  irf_norm <- normalize_irf(tall_irf, "identity")
  # calculate initial values corresponding to the normalized irf
  inval <- init_vals(irf_obj = irf_norm, nu = nu, degs = degs, tmpl_mod = tmpl_mod)
  # drop the last value of the initial values as it contains the sigma_noise prm,
  # which is calculated separately depending on the estimation scheme
  params0 <- c(utils::head(inval$prms_out, -1), mean(diag(sigma_noise)))
  return(list(params0 = params0, sigma = inval$sigma))
}

#' @title Obtain initial values to start the EM algorithm
#'
#' @description
#' \code{boot_init} implements the estimation scheme of the initial parameter
#' values given in Appendix B of \emph{Estimation of Impulse-Response
#' Functions with Dynamic Factor Models: A New Parametrization} available at
#' \url{https://arxiv.org/pdf/2202.00310.pdf}.
#'
#' @export
#'
#' @param data matrix of dimension \eqn{(T x N)}, where \eqn{T} is the number of observations and
#' \eqn{N} is the number of variables
#' @param nu a vector of Kronecker indices
#' @param degs optional, vector of length 2, where one can specify the degrees of c(z) and d(z),
#' either of the parameters must be equal to \code{max(nu)}
#' @param reg_prm integer valued parameter allowing the user to directly estimate initial values from
#' regularized data such that the std deviation of the randomly drawn noise vector of normally
#' distributed values is \code{10^-reg_prm}
#' @param nrep integer, \code{nrep=1} for no bootstrapping of initial values and \code{nrep>1}
#' for \code{nrep} bootstraps of initial values
#' @param tmpl_mod optional, supply the template as an input, using this speeds up the function,
#' enabled for convenience of the function \code{\link{do_everything_rmfd}}
#'
#' @return a list of elements
#' \item{params0}{vector of initial values}
#' \item{sigma}{estimate of initial error covariance matrix}
#' \item{llval}{the log-likelihood value of the initial state space model, returned only if \code{nrep>1}}
#'
#' @examples
#' # compare the difference in computing time
#' # w/ and w/o bootstrapping
#' # scale data first
#' Y <- scale(FRED_light$df)
#' tmp <- Sys.time()
#' arg_list <- list(data = Y, nu = c(1,1,1,1), reg_prm = NULL, nrep = 1)
#' est0 <- do.call(boot_init, arg_list)
#' time_elapsed <- Sys.time() - tmp
#' arg_list$nrep <- 10
#' tmp <- Sys.time()
#' est1 <- do.call(boot_init, arg_list)
#' time_elapsed1 <- Sys.time() - tmp
#' # the computation slows down roughly at a rate 1 sec/1 bootstrap round, and is
#' # caused by the evaluation of the log-likelihood of the boostrapped models
#' print(time_elapsed1 - time_elapsed)
#'
boot_init <- function(data, nu, degs = NULL, reg_prm = 6, nrep = 1, tmpl_mod = NULL){

  data <- as.matrix(data)
  if(is.null(degs)) degs <- rep(max(nu), 2)
  if(is.null(tmpl_mod)){
    tmpl_mod <- tmpl_rmfd_echelon_ss(dim_out = ncol(data),
                                     nu = nu,
                                     degs = degs)
    tmpl_mod$tmpl$xpx <- crossprod(tmpl_mod$tmpl$H)
  }

  boot_list <- replicate(nrep, init_wrap(data = data,
                                         nu = nu,
                                         degs = degs,
                                         tmpl_mod = tmpl_mod$tmpl,
                                         reg_prm = if(nrep>1) reg_prm else NULL))
  if(nrep>1){
    ll_vec <- sapply(1:nrep, function(x) smoothed_moments8(fill_template(boot_list[,x]$params0,
                                                                         tmpl_mod$tmpl),
                                                           Sigma = boot_list[,x]$sigma,
                                                           data_wide = t(data),
                                                           only_ll = TRUE))
    return(list(params0 = boot_list[, which.max(ll_vec)]$params0,
                sigma = boot_list[, which.max(ll_vec)]$sigma,
                llval = max(ll_vec)))
  } else{
    return(list(params0 = boot_list[, 1]$params0,
                sigma = boot_list[, 1]$sigma))
  }
}
