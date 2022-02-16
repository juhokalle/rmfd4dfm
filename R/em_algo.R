#' Estimation of the D-DFM model using the EM algorithm
#'
#' This function carries out the EM estimation of RMFD model by iterating between
#' E and M step (\code{smoothed_moments8} and \code{max_step}, resp.) until convergence
#'
#' @param params_deep_init initial system and idiosyncratic variance parameter values
#' @param sigma_init initial value of the \eqn{dim_in x dim_in} innovation covariance matrix
#' @param data_wide \eqn{n x T} matrix containing the data
#' @param tmpl_ss \code{tmpl_rmfd_echelon_ss} object specifying the model to be estimated
#' @param MAXIT maximum number of iterations
#' @param VERBOSE print estimation progress
#'
#' @keywords internal
update_em2 = function(params_deep_init, sigma_init,
                      data_wide, tmpl_ss,
                      MAXIT = 100, VERBOSE = FALSE, conv_crit = 1e-2){

  dim_in <- dim(sigma_init)[1]

  # Construct a model corresponding to the initial values
  model_this = fill_template(params_deep_init, tmpl_ss$tmpl)

  # Start with E-Step and exit if unstable
  out_e = smoothed_moments8(stsp_mod = model_this,
                            Sigma = sigma_init,
                            data_wide = data_wide)
  ll_this <- out_e$loglik
  if(!is.finite(ll_this)) stop("Initial estimation results in unstable system, the algorithm stops.")

  # Prepare while loop
  flag_converged = FALSE
  iter = 1
  conv_stat <- vector("list")

  while(!flag_converged && iter<MAXIT){

    # M-step
    out_m <- max_step(out_smooth = out_e, tmpl_ss = tmpl_ss)

    # E-step
    out_e = smoothed_moments8(stsp_mod = out_m$model_new,
                              Sigma = out_e$Q[1:dim_in,1:dim_in],
                              data_wide = data_wide)
    # check convergence
    ll_new = out_e$loglik
    if(!is.finite(ll_new)) stop(paste("Unstable system, estimation stopped at ", iter, ". iteration.", sep = ""))
    delta_loglik = abs(ll_new - ll_this)
    avg_ll = 0.5*(abs(ll_new + ll_this) + 1e-6)
    diff_ll = delta_loglik / avg_ll
    ll_this <- ll_new

    if(diff_ll < conv_crit) flag_converged = TRUE

    if(VERBOSE){
      if(flag_converged){
        cat("\nThe EM algorithm converged after ", iter, " iterations.\n", sep ="")
      } else if(iter==1){
        cat("estimation ongoing.")
      } else if(iter%%50==0){
        cat("\nestimation ongoing")
      } else if(iter==MAXIT){
        cat("\n Maximum no. iterations reached.")
      } else{
        cat(".")
      }
      conv_stat$ll_value[iter] <- ll_this
      conv_stat$conv_cr[iter] <- diff_ll
    }
    iter = iter + 1
  }

  # finalize with M-step to get parameter estimates corresponding to the last iteration
  out_m <- max_step(out_smooth = out_e, tmpl_ss = tmpl_ss)

  return(list(ssm_final = out_m$model_new,
              Sigma = out_e$Q,
              ll_val = ll_this,
              conv_stat = conv_stat))
}

#' @title Obtaining Moment Matrices for EM Algorithm Using Kalman Filter/Smoother
#'
#' @description
#' For details see section 3.2. and Appendix C in \emph{Estimation of Impulse-Response
#' Functions with Dynamic Factor Models: A New Parametrization}
#' available at \url{https://arxiv.org/pdf/2202.00310.pdf}.
#'
#' @param stsp_mod \link{stspmod} object containing the state space model used for smoothing and the idiosyncratic noise component
#' @param Sigma Matrix of dimension \code{dim_in x dim_in}. All \code{dim_in^2} elements are saved, symmetry is not taken into account.
#' @param data_wide Matrix of dimension \code{dim_out x n_obs}
#' @param only_ll return only the log-likelihood value calculated using KF
#'
#' @keywords internal
smoothed_moments8 = function(stsp_mod, Sigma, data_wide, only_ll = FALSE){

  ABCD <- stsp_mod$sys %>% unclass
  sigma_noise <- stsp_mod$sigma_L %>% unclass
  # Integer-valued params
  dim_in = dim(Sigma)[1]
  n_obs = dim(data_wide)[2]
  dim_out = attr(stsp_mod$sys, "order")[1]
  dim_state = attr(stsp_mod$sys, "order")[3]

  # Extract the state space matrices
  A = ABCD[1:dim_state, 1:dim_state]
  B = ABCD[1:dim_state, (dim_state+1):(dim_state+dim_in)]
  C = ABCD[(dim_state+1):(dim_state+dim_out), 1:dim_state]
  D = ABCD[(dim_state+1):(dim_state+dim_out), (dim_state+1):(dim_state+dim_out)]

  # Derived quantities
  Q = B%*%Sigma%*%t(B) # Sigma is \Sigma
  R = sigma_noise

  # call Kalman smoother
  # return early if only the log-likelihood value is of interest
  if(only_ll){
    return(kfks_cpp(A, C, R, Q,
                    data_wide, lyapunov(A, Q), dim_state,
                    only_ll = TRUE)$loglik)
  } else{
    out_kfks = kfks_cpp(A, C, R, Q,
                        data_wide, lyapunov(A, Q), dim_state,
                        only_ll = FALSE)
  }
  P_t_T = out_kfks$P_t_T
  C_t_T = out_kfks$C_t_T
  s_t_T = out_kfks$s_t_T
  e_t_T = out_kfks$e_t_T
  v_t_T = out_kfks$v_t_T

  # Calculate moments
  mean_P_t_T = apply(P_t_T[, , 2:n_obs, drop = FALSE],
                     MARGIN = c(1,2),
                     mean) %>% zapsmall()
  mean_P_l1_T = apply(P_t_T[, , 1:(n_obs-1), drop = FALSE],
                      MARGIN = c(1,2),
                      mean) %>% zapsmall()
  mean_C_t_T = apply(C_t_T[, , (2:n_obs), drop = FALSE],
                     MARGIN = c(1,2),
                     mean) %>% zapsmall()

  # eq. 19 WE'83 JoE
  R = 1/n_obs * tcrossprod(e_t_T) + C %*% mean_P_t_T %*% t(C)
  # we model the idiosyncratic variance homoskedastic and serially uncorrelated
  sigma_noise <- diag(mean(diag(R)), dim_out, dim_out)

  # eq. 20 WE'83 JoE
  Q = 1/n_obs * tcrossprod(v_t_T) + mean_P_t_T + A %*% mean_P_l1_T %*% t(A) - A %*% mean_C_t_T - t(mean_C_t_T) %*% t(A)

  # sample variance of lagged smoothed states, i.e. s_{t-1|T}
  var_s_l1_T <- 1/(n_obs-1) * s_t_T[,1:(n_obs-1)] %*% t(s_t_T[,1:(n_obs-1)])
  # sample covariance of lagged and current smoothed states, i.e s_{t-1|T} and s_{t|T}
  cov_s_l1_t_T <- 1/(n_obs-1) * s_t_T[,1:(n_obs-1)]%*%t(s_t_T[,2:n_obs])

  # eq. 17 WE'83 JoE: conditional variance of s_{t-1} given data and theta
  var_ss_l1 <- var_s_l1_T + mean_P_l1_T
  # eq. 18 WE'83 JoE: conditional covariance of s_t and s_{t-1} given data and theta
  cov_ss_l1_t <- cov_s_l1_t_T + mean_C_t_T

  # sample variance of smoothed states, i.e. s_{t|T}
  var_s_t_T <- 1/(n_obs-1) * s_t_T[,2:n_obs]%*%t(s_t_T[,2:n_obs])
  # eq. 14 WE'83 JoE: sample covariance of smoothed states s_{t|T} and data y_t:
  cov_s_t_T_yt <- 1/n_obs * s_t_T%*%t(data_wide)
  cov_s_tl1_T_yt <- 1/(n_obs-1) * s_t_T[,1:(n_obs-1)] %*% t(data_wide[,2:(n_obs)])

  # eq. 16 WE'83 JoE: conditional variance of s_t given data and theta
  var_ss_t <- var_s_t_T + mean_P_t_T

  return(list(Q = Q,
              R = R,
              sigma_noise = sigma_noise,
              var_ss_l1 = var_ss_l1,
              cov_ss_l1_t = cov_ss_l1_t,
              var_ss_t = var_ss_t,
              cov_s_t_T_yt = cov_s_t_T_yt,
              loglik = out_kfks$loglik,
              dim_in = dim_in,
              dim_out = dim_out,
              dim_state = dim_state))
}

#' The maximization step associated with the EM algorithm
#'
#' For details see section 3.2. in \emph{Estimation of Impulse-Response
#' Functions with Dynamic Factor Models: A New Parametrization}
#' available at \url{https://arxiv.org/pdf/2202.00310.pdf}.
#'
#' @param out_smooth the output of \code{out_smooth}
#' @param tmpl_ss model template
#'
#' @return a list of elements
#' \item{model_new}{a state space model corresponding to the new estimated parameters}
#' \item{params_new}{a vector of deep parameters}
#'
max_step <- function(out_smooth, tmpl_ss){

  dim_in <- out_smooth$dim_in
  dim_out <- out_smooth$dim_out
  dim_state <- out_smooth$dim_state

  H_A = tmpl_ss$A$H_A
  h_A = tmpl_ss$A$h_A
  H_C = tmpl_ss$C$H_C
  h_C = tmpl_ss$C$h_C
  H = tmpl_ss$tmpl$H
  h = tmpl_ss$tmpl$h
  xpx = tmpl_ss$tmpl$xpx

  # GLS regression for deep parameters in A
  # restrict all the other prms of Q to zero apart from those in the upper left dim_in x dim_in block
  Q_inv = diag(1, dim_state, dim_in) %*% solve(out_smooth$Q[1:dim_in,1:dim_in]) %*% diag(1, dim_in, dim_state)
  var_ss_l1 = out_smooth$var_ss_l1
  cov_ss_l1_t = out_smooth$cov_ss_l1_t

  # setting the tolerance value in the generalized inverse function below helps to avoid
  # numerical instabilities associated with the ill-conditioned matrix subject to inversion
  deep_A1 <- corpcor::pseudoinverse(t(H_A) %*% (var_ss_l1 %x% Q_inv) %*% H_A, tol = 1e-2)
  deep_A2_1 <- t(H_A) %*% (cov_ss_l1_t %x% Q_inv) %*% c(diag(dim_state))
  deep_A2_2 <- t(H_A) %*% (var_ss_l1 %x% Q_inv) %*% h_A
  deep_A <- deep_A1 %*% (deep_A2_1 - deep_A2_2)

  A = matrix(h_A + H_A %*% deep_A, nrow = dim_state, ncol = dim_state)

  # GLS regression for deep parameters in C
  R_inv <- solve(out_smooth$sigma_noise)
  var_ss_t = out_smooth$var_ss_t
  cov_s_t_T_yt = out_smooth$cov_s_t_T_yt

  deep_C1 <- solve(t(H_C) %*% (var_ss_t %x% R_inv) %*% H_C)
  deep_C2_1 <- t(H_C) %*% (cov_s_t_T_yt %x% R_inv) %*% c(diag(dim_out))
  deep_C2_2 <- t(H_C) %*% (var_ss_t %x% R_inv) %*% h_C
  deep_C <- deep_C1 %*% (deep_C2_1-deep_C2_2)

  C = matrix(h_C + H_C %*% deep_C, nrow = dim_out, ncol = dim_state)

  # Extract deep parameters
  BD = rbind(cbind(diag(1, dim_state, dim_in), matrix(0, dim_state, dim_out-dim_in)),
             diag(dim_out))
  ABCD = cbind(rbind(A, C), BD)

  prm_vec <- c(as.vector(ABCD), as.vector(out_smooth$sigma_noise))
  params_new <- solve(xpx, crossprod(H, prm_vec - h)) # this line extracts the deep parameters

  model_new <- fill_template(params_new, tmpl_ss$tmpl)

  return(list(model_new = model_new,
              params_new = params_new))
}

