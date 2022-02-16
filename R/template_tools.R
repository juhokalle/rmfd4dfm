#' @title Template for echelon RMFD model in state space format
#'
#' @description
#' In this template one obtains a pure echelon form by supplying the vector
#' of Kronecker indices and leaving the degrees of c(z) and d(z) unspecified.
#' On the other hand, the (P,Q)-form is obtained by supplying a vector
#' of Kronecker indices where all the values are equal to the maximum between
#' the degrees of c(z) and d(z), and then specifying the polynomial degrees
#' c(z) and d(z). Finally, a hybrid version of these two forms is obtained
#' by setting the Kronecker indices and polynomial degrees freely with
#' the restriction that \eqn{max(nu) = max(P, Q)}.
#'
#' @param dim_out output dimension, i.e. number of variables
#' @param nu vector of Kronecker indices
#' @param degs vector of length 2, where one can optionally specify the degrees of
#'   \code{c(z)} and \code{d(z)}, either of the parameters must be equal to max(nu)
#'
#' @keywords internal
tmpl_rmfd_echelon_ss = function(dim_out, nu, degs = NULL) {

  if(is.null(degs)) degs <- rep(max(nu), 2)
  if(!max(degs)==max(nu)) stop("Polynomial degrees not compatible with the Kronecker indices")
  dim_in = length(nu)
  deg_c <- degs[1]
  deg_d <- degs[2]
  m = max(nu)

  # code the position of the basis columns of the Hankel matrix
  basis = nu2basis(nu)

  # coefficients of c(z) in reverse order!!!
  # c = [c[p]',...,c[0]']'
  c_ech = matrix(0, nrow = dim_in*(m+1), ncol = dim_in)
  # d = [d[0]',...,d[p]']'
  d_ech = matrix(0, nrow = dim_out*(m+1), ncol = dim_in)

  # code free entries with NA's
  for (i in (1:dim_in)) {
    shift = (m-nu[i])*dim_in
    k = nu[i]*dim_in + i
    basis_i = basis[basis < k]
    # cat(i, shift, k, basis_i, '\n')
    c_ech[k + shift, i] = 1
    c_ech[basis_i + shift, i] = NA_real_
    d_ech[iseq(dim_out + 1, (nu[i] + 1)*dim_out), i] = NA_real_   # i-th column has degree nu[i]
    d_ech[iseq(dim_in + 1, dim_out), i] = NA_real_               # the last (m-n) rows of b[0] are free
  }

  if(deg_c < max(nu)) c_ech[1:(dim_in*(m-deg_c)),] <- 0 # c(z) coeffs in reverse order
  if(deg_d < max(nu)) d_ech[(dim_out*(deg_d+1)+1):(dim_out*(m+1)),] <- 0

  # state space implementation
  m <- max(nu, deg_d+1)
  A <- matrix(0, nrow = dim_in, ncol = m*dim_in)
  C <- matrix(0, nrow = dim_out, ncol = m*dim_in)

  for(jj in 1:m){
    col_ix <- (1+(jj-1)*dim_in):(jj*dim_in)
    row_ix <- ((m-jj)*dim_in+1):((m+1-jj)*dim_in) # c[i] stored in reverse order
    ix <- (1+(jj-1)*dim_out):(jj*dim_out)
    C[,col_ix] <- d_ech[ix,]
    if(jj>1 || deg_c > deg_d){
      # skip c[0] and fix it into C after loop
      if(deg_c > deg_d) A[,col_ix] <- c_ech[row_ix,] else A[,col_ix-dim_in] <- c_ech[row_ix,]
    }
  }

  # c[0] = d[0][1:dim_in,1:dim_in]
  C[1:dim_in, 1:dim_in] <- utils::tail(c_ech, dim_in) # this is vague...
  A <- rbind(A, diag(1, nrow = (m-1)*dim_in, ncol = m*dim_in))
  B = cbind(diag(x = 1, nrow = m * dim_in, ncol = dim_in),
            matrix(0, nrow = m * dim_in, ncol = dim_out-dim_in))
  D = diag(dim_out)
  ABCD <- cbind(rbind(A, C), rbind(B, D))
  sys = structure(ABCD,
                  order = c(dim_out, dim_out, dim_in*m), class = c('stsp','ratm'))

  # create a helper model
  model = list(sys = sys, sigma_L = diag(dim_out), names = NULL, label = NULL)
  model = structure(model, class = c('stspmod', 'rldm'))
  tmpl <- model2template(model)

  # change Sigma to \sigma^2*I_n
  ix <- (length(ABCD)+1):(length(ABCD)+dim_out^2)
  tmpl$H <- bdiag(tmpl$H[-ix, ], matrix(diag(dim_out), nrow = dim_out^2, ncol = 1))
  tmpl$h[ix] <- 0
  tmpl$n.par <- tmpl$n.par + 1

  # Affine transformation: A
  h_A = as.vector(A)
  ix = which(!is.finite(h_A))
  h_A[ix] = 0
  n_A = length(ix)

  H_A = matrix(0, nrow = length(h_A), ncol = n_A)
  H_A[ix,] = diag(n_A)

  # Affine transformation: C
  h_C <- as.vector(C)
  ix <- which(!is.finite(h_C))
  h_C[ix] <- 0
  n_C <- length(ix)

  H_C <- matrix(0, nrow = length(h_C), ncol = n_C)
  H_C[ix,] <- diag(n_C)

  return(list(A = list(H_A = H_A,
                       h_A = h_A,
                       n_A = n_A),
              C = list(H_C = H_C,
                       h_C = h_C,
                       n_C = n_C),
              tmpl = tmpl)
  )
}

#' @title Represent RMFD model in state space format
#'
#' @param rmfdsys \code{rmfd} object
#' @param nu vector of Kronecker indices
#'
#' @return \code{ABCD}, matrix of dimensions \eqn{r+n x r+n} of the corresponding state space model,
#' where state and output dimensions are denoted by \eqn{r} and \eqn{n}.
#'
#' @keywords internal
rmfd2stsp = function(rmfdsys, nu){

  # Extract integer valued parameters
  dim_out = attr(rmfdsys, "order")[1]
  dim_in = attr(rmfdsys, "order")[2]
  deg_c = attr(rmfdsys, "order")[3]
  deg_d = attr(rmfdsys, "order")[4]
  if(!max(deg_c,deg_d)==max(nu)) stop("Polynomial degrees not compatible with the Kronecker indices")
  m = max(nu, deg_d+1)

  # stsp matrices
  if(deg_d>=deg_c){
    A = dbind(3, unclass(rmfdsys$c), array(0, dim = c(dim_in, dim_in, (deg_d-deg_c+1)))) %>% companion_matrix()
    C <- rmfdsys$d %>% unclass
    dim(C) <- c(dim_out, m*dim_in)
  } else{
    A = dbind(3, unclass(rmfdsys$c)) %>% companion_matrix()
    C <- rmfdsys$d %>% unclass
    dim(C) <- c(dim_out, (deg_d+1)*dim_in)
    C <- cbind(C, matrix(0, dim_out, dim_in*(deg_c-deg_d-1)))
  }


  B = cbind(diag(x = 1, ncol = dim_in, nrow = m * dim_in),
            matrix(0, ncol = dim_out-dim_in, nrow = m * dim_in))
  D = diag(dim_out)

  ABCD = cbind(rbind(A,C), rbind(B,D))

  return(ABCD)
}

#' @title Obtain RMFD object from a state space system
#'
#' @param stspsys matrix corresponding to ABCD state space system
#' @param dim_out output dimension, i.e. number of variables
#' @param nu vector of Kronecker indices
#' @param degs vector of length 2, where one can optionally specify the degrees of \code{c(z)}
#' and \code{d(z)}, either of the parameters must be equal to \code{max(nu)}
#'
#' @return \code{rmfd} object corresponding to the state space setup
#'
#' @keywords internal
stsp2rmfd = function(stspsys, dim_out, nu, degs = NULL){

  if(is.null(degs)) degs <- rep(max(nu), 2)
  if(!max(degs)==max(nu)) stop("Polynomial degrees not compatible with the Kronecker indices")
  dim_in = length(nu)
  deg_c <- degs[1]
  deg_d <- degs[2]
  m = max(nu,deg_d+1)
  dim_state <- m*dim_in
  # get polynomial matrices c(z) and d(z) by arraying stsp matrices A and C
  deep_C <- stspsys %>% unclass %>% .[1:dim_in, 1:(dim_in*deg_c)]
  polm_C <- c(diag(dim_in), -deep_C) %>%
    array(dim=c(dim_in,dim_in,deg_c+1)) %>%
    polm()
  deep_D <- stspsys %>%
    unclass %>%
    .[(dim_state+1):(dim_out+dim_state), 1:(dim_in*(deg_d+1))] %>% c()
  polm_D <- deep_D %>% array(dim=c(dim_out, dim_in, deg_d+1)) %>% polm()
  # make sure d[0][1:dim_in, 1:dim_in]=c[0]
  d_0 <- polm_D %>% unclass %>% .[1:dim_in,1:dim_in,1]
  rmfd_out <- rmfd(d_0%r%polm_C, polm_D)
  return(rmfd_out)
}
