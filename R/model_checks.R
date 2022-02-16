#' @title Bai & Ng information criteria for determining the number of static factors
#'
#' @description
#' This function calculates two different information criteria for
#' determining the number of static factors, which were introduced in
#' Bai & Ng (2002). The function gives two different sets of
#' information criteria (IC and PC), which use three different
#' penalty functions to take into account the increased complexity of the
#' model.
#'
#' @export
#'
#' @seealso Bai, J., & Ng, S. (2002). Determining the number of factors in approximate factor models.
#' Econometrica, 70(1), 191-221.
#'
#' @param X \eqn{T x N} data matrix
#' @param rmax maximum number of static factors to be considered
#' @param rvar fix penalty function scaling term, defaults to \code{rvar=rmax}
#'
#' @return lists of IC and PC, where the values indicate the suggested number of factors
#'
#' @examples
#' (bn1 <- baingcriterion(X = FRED_heavy$df, 25))
#' (bn2 <- baingcriterion(X = FRED_light$df, 25))
baingcriterion <- function(X, rmax, rvar = NULL)
{
  X <- as.matrix(X)
  if(is.null(rvar)) rvar = rmax
  nobs = dim(X)[1]
  nvar = dim(X)[2]
  mu <- matrix(1, nrow=nobs,ncol=1)%*%apply(X, 2, mean)
  Winv <- diag(apply(X, 2, sd)^-1)
  x = (X-mu)%*%Winv

  Gamma0 = cov(x)
  H <- eigen(Gamma0)$vectors
  V <- vector()
  for(j in 1:rmax){
    I = x%*%(diag(nvar) - tcrossprod(H[,1:j]))
    V[j] = sum(I^2)/(nvar*nobs)
  }
  e = (nvar + nobs)/(nvar*nobs)

  penalty1 = (1:rmax)*e*log(1/e)
  penalty2 = (1:rmax)*e*log(min(nvar,nobs))
  penalty3 = (1:rmax)*(log(min(nvar,nobs))/min(nvar,nobs))

  IC = cbind(log(V) + penalty1, log(V) + penalty2, log(V) + penalty3)
  PC = cbind(V + V[rvar]*penalty1, V + V[rvar]*penalty2, V + V[rvar]*penalty3)
  IC.min <- apply(IC, 2, which.min)
  PC.min <- apply(PC, 2, which.min)
  return(list(IC = IC.min, PC = PC.min))
}

#' @title ABC criterion for choosing the number of static factors
#'
#' @description
#' The function estimates the number of static factors in the
#' dynamic factor model. These information criteria are based on
#' Alessi, Barigozzi and Capasso (2010). This function is a direct
#' translation of the MATLAB code given in \url{http://www.barigozzi.eu/Codes.html}
#' written by Matteo Barigozzi.
#'
#' @export
#'
#' @seealso Alessi, L., M. Barigozzi, and M. Capasso (2010).
#' Improved penalization for determining the number of factors
#' in approximate static factor models.
#' Statistics and Probability Letters 80, 1806-1813.
#'
#' @param X \eqn{T x N} data matrix
#' @param kmax maximum number of factors
#' @param nbck number of sublocks to be used, defaults to \code{floor(N/10)}
#' @param cmax maximum value for the penalty constant, default is 3
#'
#' @return
#' \item{rhat1}{determines the number of shocks using a large window}
#' \item{rhat2}{determines the number of shocks using a small window}
#'
#' @examples
#' (abc1 <- abc_crit(X = FRED_heavy$df, kmax = 25))
#' (abc2 <- abc_crit(X = FRED_light$df, kmax = 25))
abc_crit <- function(X, kmax, nbck = NULL, cmax = 3)
{
  npace <- 1
  step <- 500
  nobs <- dim(X)[1]
  nvar <- dim(X)[2]

  if(is.null(nbck)) nbck <- floor(nvar/10)
  x = scale(X, center = TRUE, scale = TRUE)

  # running the criterion
  s <- 0
  abc <- matrix(0, 0, floor(cmax*step))
  IC1 <- matrix(0, nrow = kmax + 1, ncol = 1)
  for(N in seq(from = nvar-nbck, to = nvar, by = npace)){
    s <- s+1
    Ns <- order(stats::runif(nvar))
    xs <- x[1:nobs, Ns[1:N]]
    xs <- scale(xs, center = TRUE, scale = TRUE)
    eigv <- eigen(cov(xs), only.values = TRUE)$values
    for(k in 1:(kmax+1)){
      IC1[k,1] <- sum(eigv[k:N])
    }
    p <- ((N+nobs)/(N*nobs))*log((N*nobs)/(N+nobs))
    T0 <- matrix(0:kmax, nrow = kmax + 1, ncol = 1)*p
    abc <- rbind(abc, matrix(0, 1, floor(cmax*step)))
    for(c in 1:floor(cmax*step)){
      cc <- c/step
      IC <- IC1/N + T0*cc
      rr <- apply(IC, 2, which.min)
      abc[s,c] <- rr-1
    }
  }
  # select automatically the number of shocks
  cr <- (1:floor(cmax*500))/500
  ABC <- matrix(0, length(cr), 4)

  for(ll in 1:2){
    # don't know what's the point of the loop here, just copied it
    ABC[1,1] <- kmax
    ABC[1,2:3] <- 0
  }
  sabc <- apply(abc, 2, sd)
  c1 <- 2
  for(ii in 1:length(cr)){
    if(sabc[ii]==0){
      if(abc[nrow(abc),ii]==ABC[c1-1,1]){
        ABC[c1-1,3] <- cr[ii]
      } else{
        ABC[c1,1] <- abc[nrow(abc), ii]
        ABC[c1, 2:3] <- cr[ii]
        c1 <- c1+1
      }
    }
  }
  ABC[,4] <- ABC[,3]-ABC[,2]
  q <- ABC[which(ABC[2:nrow(ABC),4]>0.05)+1,1]
  rhat1 <- q[1]
  q <- ABC[which(ABC[2:nrow(ABC),4]>0.01)+1,1]
  rhat2 <- q[1]
  return(list(rhat1 = rhat1, rhat2 = rhat2))
}

#' @title Admissible models
#'
#' @description
#' \code{admissible_mods} gives the qualified models according to
#' the model selection criteria in section 3.3 of \emph{Estimation of
#' Impulse-Response Functions with Dynamic Factor Models: A New Parametrization}
#' available at \url{https://arxiv.org/pdf/2202.00310.pdf}.
#'
#' @export
#'
#' @note This is a rather crude way of obtaining the set of models satisfying the
#' model selection criteria because the function simulates
#' the model parameters from a \emph{U(-1,1)}-distribution and then
#' runs the checks. This entails that the function may not always find the
#' maximum number of qualified models, and this being the case especially with
#' specifications where the dynamic factor dimension and maximum Kronecker index are
#' set to a high value. However, for a small values of these parameters the function
#' should work fine.
#'
#' @param q input dimension, number of shocks
#' @param max_nu maximum Kronecker index value
#' @param fix_r boolean, whether the state dimension of the state space
#' representation is fixed to \code{q*max_nu}
#'
#' @return a list of two lists of equal length, where each element in the list
#' correspond either to the Kronecker index specification or to the additional
#' possible restrictions imposed on the lag orders of \eqn{c(z)} and \eqn{d(z)}
#' for given DFM specification. The elements are saved in named lists
#' \item{nus}{list of Kronecker index vectors}
#' \item{degs}{list of polynomial degrees}
#'
#' @examples
#' tmp <- admissible_mods(4, 4)
#' sample(1:length(tmp$nus), 1) %>% sapply(function(x)
#'  cat("The qualified model no. ", x,
#'      "has a Kronecker index structure ", tmp$nus[[x]], "with the
#'  orders of the lag polynomials restricted to deg(c(z)) = ", tmp$degs[[x]][1],
#'      " and deg(d(z)) = ", tmp$degs[[x]][2], "\n")) -> tmp
admissible_mods <- function(q, max_nu, fix_r = TRUE){

  nus <- replicate(q, list(0:max_nu)) %>%
    expand.grid() %>%
    apply(1, sort) %>%
    as.data.frame() %>%
    as.list() %>%
    unique()
  degs <- replicate(2, list(1:max_nu)) %>%
    expand.grid() %>%
    t() %>%
    as.data.frame() %>%
    as.list()
  int_list <- expand.grid(nus, degs)
  int_ix <- sapply(int_list[[1]], max)==sapply(int_list[[2]], max)
  int_list <- int_list[int_ix, ]
  flag <- rep(NA, nrow(int_list))
  for(i in 1:nrow(int_list)){

    tmpl_i <- tmpl_rmfd_echelon_ss(dim_out = q+1,
                                   nu = unlist(int_list[i,]$Var1),
                                   degs = unlist(int_list[i,]$Var2))
    mod_i <- fill_template(th = runif(tmpl_i$tmpl$n.par, -1, 1),
                           template = tmpl_i$tmpl)
    ctr_mat <- ctr_matrix(mod_i$sys)
    obs_mat <- obs_matrix(mod_i$sys)
    flag[i] <- min(dim(ctr_mat))==qr(ctr_mat)$rank && min(dim(obs_mat))==qr(obs_mat)$rank
    if(fix_r) flag[i] <- flag[i] && dim(mod_i$sys)[3] == q*max_nu
  }
  out_list <- list(nus = int_list[flag,]$Var1,
                   degs = int_list[flag,]$Var2)
  return(out_list)
}
