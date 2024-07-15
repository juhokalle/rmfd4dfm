# as_methods.R ##################################################
#

as.stspmod = function(obj, ...){
  UseMethod("as.stspmod", obj)
}


as.stspmod.armamod = function(obj, ...){
  obj = unclass(obj)
  obj$sys = as.stsp(obj$sys)
  class(obj) = c('stspmod', ' rldm')
  return(obj)
}


# autocov_methods.R ######################################
# defines the main methods for the 'autocov' class


autocov = function(obj, type, ...) {
  UseMethod("autocov", obj)
}

autocov.default = function(obj, type=c('covariance','correlation','partial'), lag.max = NULL,
                           na.action = stats::na.fail, demean = TRUE, ...) {

  if (!is.null(lag.max)) {
    lag.max = as.integer(lag.max)[1]
    if (lag.max < 0) stop('negative lag.max')
  }
  type = match.arg(type)

  out_acf = try(stats::acf(obj, lag.max = lag.max, type = 'covariance', plot = FALSE,
                           na.action = na.action, demean = demean))

  if (inherits(out_acf, 'try-error')) stop('stats:acf failed, the input "obj" may not be supported.')

  gamma = aperm(out_acf$acf, c(2,3,1))
  n.obs = out_acf$n.used

  if (type == 'covariance') {
    acf = gamma
  }
  if (type == 'correlation') {
    # compute correlations
    m = dim(gamma)[1]
    # standard deviations
    s = sqrt(diag(matrix(gamma[,,1], nrow = m, ncol = m)))
    s = matrix(s, nrow = m, ncol = m)
    s2 = s * t(s)
    acf = gamma
    for (i in (1: dim(gamma)[3])) acf[,,i] = acf[,,i] / s2
  }
  if (type == 'partial') {
    acf = est_ar_dlw(gamma, penalty = -1)$partial
  }
  acf = structure(acf, class = c('pseries', 'ratm') )

  out = structure( list(acf = acf, type = type, gamma = gamma, names = NULL, label = NULL,  n.obs = n.obs),
                   class = c('autocov', 'rldm') )

  return(out)
}


autocov.autocov = function(obj, type=c('covariance','correlation','partial'), ...) {
  type = match.arg(type)
  out = unclass(obj)
  gamma = out$gamma

  out$type = type

  if (type == 'covariance') {
    acf = gamma
  }
  if (type == 'correlation') {
    # compute correlations
    m = dim(gamma)[1]
    # standard deviations
    s = sqrt(diag(matrix(gamma[,,1], nrow = m, ncol = m)))
    s = matrix(s, nrow = m, ncol = m)
    s2 = s * t(s)
    acf = gamma
    for (i in (1: dim(gamma)[3])) acf[,,i] = acf[,,i] / s2
  }
  if (type == 'partial') {
    acf = est_ar_dlw(gamma, penalty = -1)$partial
  }
  acf = structure(acf, class = c('pseries', 'ratm') )
  out$acf = acf

  out = structure( out, class = c('autocov', 'rldm') )

  return(out)
}


autocov.armamod = function(obj, type=c('covariance','correlation','partial'), lag.max = 12, ...) {
    lag.max = as.integer(lag.max[1])
    if (lag.max<0) stop('negative lag.max')
    type = match.arg(type)

    d = unname(dim(obj$sys))
    m = d[1]
    n = d[2]
    p = d[3]
    q = d[4]

    a = unclass(obj$sys$a)
    b = unclass(obj$sys$b)
    sigma = obj$sigma_L
    sigma = sigma %*% t(sigma)
    k = unclass(pseries(obj$sys, lag.max = q)) # impulse response

    # the ACF is computed via the Yule Walker equations (for j = 0, ... , p)
    # a[0] gamma[j] + a[1] gamma[j-1] + ... + a[p] gamma[j-p] =
    #         = E (b[0] eps[t] + b[1] eps[t-q] + ... + b[q] eps[t-q]) y[t-j]'

    # compute "right hand side" of the YW equations
    # E (b[0] eps[t] + b[1] eps[t-q] + ... + b[q] eps[t-q]) y[t-j]' =
    #    b[j] sigma k[0]' + ... + b[q] sigma k[q-j]'
    Ewy = array(0, dim = c(m, m, (max(p , q, lag.max)+1)) )
    for (j in (0:q)) {
      for (i in (j:q)) Ewy[,,j+1] = Ewy[,,j+1] + b[,,i+1] %*% sigma %*% t(k[,,i-j+1])
    }

    # construct the Yule-Walker equations in vectorized form:
    # note: vec(a[i] * gamma[j]) = (I x a[i]) vec(gamma[j]), where x stands for the Kronecker product
    # and:  vec(a[i] * gamma[-j]) = (I x a[i]) vec(gamma[-j]) = (I x a[i]) P vec(gamma[j])
    #       where P is a permutation matrix!

    # construct permutation matrix P
    junk = as.vector(t(matrix(1:(m^2), nrow = m, ncol = m)))
    P = numeric(m^2)
    P[junk] = 1:(m^2)

    # construct an array with the Kronecker products (I x a[i])
    A = array(0, dim = c(m^2, m^2, p+1))
    for (i in iseq(0, p)) A[,,i+1] = diag(1, nrow = m) %x% a[ , , i + 1]

    # left hand side of equation system
    L = array(0, dim = c(m^2, m^2, p+1, p+1))
    # right hand side of equation system
    R = array(0, dim = c(m, m, p+1))
    for (j in (0:p)) {
      # E y[t] t(y[t-j]) = ...
      R[ , , j+1] = Ewy[ , , j+1]

      L[ , , j+1, j+1] = A[ , , 1]
      for (i in (0:p)) {
        lag = j - i
        #      cat(j, i, lag,'\n')
        if (lag >= 0) {
          L[ , , j+1, lag+1] = A[ , , i+1]
        } else {
          L[ , , j+1, 1-lag] = L[ , , j+1, 1-lag] + A[, P, i+1]
        }
      }
    }
    # print(L)
    # print(R)
    L = bmatrix(L, rows = c(1,3), cols = c(2,4))
    R = as.vector(R)
    # print(cbind(L,R))

    gamma0 = solve(L, R)
    #  print(gamma0)
    gamma0 = array(gamma0, dim=c(m, m, p+1))

    if (p >= lag.max) {
      gamma = gamma0[ , , 1:(lag.max+1), drop=FALSE]
    } else {
      # extend acf for lags p+1,...lag.max
      # simply by two nested for loops!
      # however we have to transform the AR coefficients
      # a[0] y_t + a[1] y[t-1] + ...  ==> y[t] = (-a[0]^{-1} a[1])y[t-1] + ...
      for (i in iseq(1,p)) {
        a[,,i+1] = solve(a[,,1], -a[,,i+1])
      }

      gamma = dbind(d = 3, gamma0, Ewy[ , , (p+2):(lag.max+1), drop = FALSE])
      for (lag in ((p+1):lag.max)) {
        for (i in iseq(1,p)) {
          gamma[ , , lag+1] = gamma[ , , lag+1] + a[ , ,i+1] %*% gamma[ , , lag+1-i]
        }
      }
    }
    #  print(gamma)

    if (type == 'covariance') {
      acf = gamma
    }
    if (type == 'correlation') {
      # compute correlations
      m = dim(gamma)[1]
      # standard deviations
      s = sqrt(diag(matrix(gamma[,,1], nrow = m, ncol = m)))
      s = matrix(s, nrow = m, ncol = m)
      s2 = s * t(s)
      acf = gamma
      for (i in (1: dim(gamma)[3])) acf[,,i] = acf[,,i] / s2
    }
    if (type == 'partial') {
      acf = est_ar_dlw(gamma, penalty = -1)$partial
    }
    acf = structure(acf, class = c('pseries', 'ratm') )

    out = structure( list(acf = acf, type = type, gamma = gamma, names = NULL, label = NULL,  n.obs = NULL),
                     class = c('autocov', 'rldm') )
    return(out)
}

autocov.stspmod = function(obj, type=c('covariance','correlation','partial'), lag.max = 12, ...) {

  lag.max = as.integer(lag.max)[1]
  if (lag.max < 0) stop('negative lag.max')
  type = match.arg(type)

  A = obj$sys$A
  B = obj$sys$B
  C = obj$sys$C
  D = obj$sys$D
  m = nrow(C)
  n = ncol(B)
  s = ncol(A)

  sigma = obj$sigma_L %*% t(obj$sigma_L)

  gamma = array(0, dim = c(m, m, lag.max+1))
  if (s == 0) {
    gamma[,,1] = D %*% sigma %*% t(D)
  } else {
    P = lyapunov(A, B %*% sigma %*% t(B), non_stable = 'stop')
    M = A %*% P %*% t(C) + B %*% sigma %*% t(D)

    gamma[,,1] = C %*% P %*% t(C) + D %*% sigma %*% t(D)
    for (i in iseq(1, lag.max)) {
      gamma[,,i+1] = C %*% M
      M = A %*% M
    }
  }

  if (type == 'covariance') {
    acf = gamma
  }
  if (type == 'correlation') {
    # compute correlations
    m = dim(gamma)[1]
    # standard deviations
    s = sqrt(diag(matrix(gamma[,,1], nrow = m, ncol = m)))
    s = matrix(s, nrow = m, ncol = m)
    s2 = s * t(s)
    acf = gamma
    for (i in (1: dim(gamma)[3])) acf[,,i] = acf[,,i] / s2
  }
  if (type == 'partial') {
    acf = est_ar_dlw(gamma, penalty = -1)$partial
  }
  acf = structure(acf, class = c('pseries', 'ratm') )

  out = structure( list(acf = acf, type = type, gamma = gamma, names = NULL, label = NULL,  n.obs = NULL),
                   class = c('autocov', 'rldm') )
  return(out)
}



# classes.R ######################################
# defines the main classes


armamod = function(sys, sigma_L = NULL, names = NULL, label = NULL) {
  if (!inherits(sys, 'lmfd')) stop('"sys" must be an lmfd object')

  d = dim(sys)
  m = d[1]
  n = d[2]

  if (is.null(sigma_L)) sigma_L = diag(n)
  if (!is.numeric(sigma_L)) stop('parameter sigma_L is not numeric')
  if ( is.vector(sigma_L) ) {
    if (length(sigma_L) == n) sigma_L = diag(sigma_L, nrow = n, ncol = n)
    if (length(sigma_L) == (n^2)) sigma_L = matrix(sigma_L, nrow = n, ncol = n)
  }
  if ( (!is.matrix(sigma_L)) || any(dim(sigma_L) != n) ) {
    stop('"sigma_L" is not compatible')
  }

  x = list(sys = sys, sigma_L = sigma_L, names = names, label = label)
  x = structure(x, class = c('armamod', 'rldm'))
  return(x)
}

rmfdmod = function(sys, sigma_L = NULL, names = NULL, label = NULL) {
  # Check sys input
  if (!inherits(sys, 'rmfd')) stop('"sys" must be an rmfd object')

  d = dim(sys)
  m = d[1]
  n = d[2]

  # Parameterisation of sigma through left-factor: Different cases
  if (is.null(sigma_L)) sigma_L = diag(n)
  if (!is.numeric(sigma_L)) stop('parameter sigma_L is not numeric')
  if ( is.vector(sigma_L) ) {
    if (length(sigma_L) == n) sigma_L = diag(sigma_L, nrow = n, ncol = n)
    if (length(sigma_L) == (n^2)) sigma_L = matrix(sigma_L, nrow = n, ncol = n)
  }
  if ( (!is.matrix(sigma_L)) || any(dim(sigma_L) != n) ) {
    stop('"sigma_L" is not compatible')
  }

  x = list(sys = sys, sigma_L = sigma_L, names = names, label = label)
  x = structure(x, class = c('rmfdmod', 'rldm'))
  return(x)
}

stspmod = function(sys, sigma_L = NULL, names = NULL, label = NULL) {
  if (!inherits(sys, 'stsp')) stop('"sys" must be an stsp object')

  d = dim(sys)
  m = d[1]
  n = d[2]

  if (is.null(sigma_L)) sigma_L = diag(n)
  if (!is.numeric(sigma_L)) stop('parameter sigma_L is not numeric')
  if ( is.vector(sigma_L) ) {
    if (length(sigma_L) == n) sigma_L = diag(sigma_L, nrow = n, ncol = n)
    if (length(sigma_L) == (n^2)) sigma_L = matrix(sigma_L, nrow = n, ncol = n)
  }
  if ( (!is.matrix(sigma_L)) || any(dim(sigma_L) != n) ) {
    stop('"sigma_L" is not compatible')
  }

  x = list(sys = sys, sigma_L = sigma_L, names = names, label = label)
  x = structure(x, class = c('stspmod', 'rldm'))
  return(x)
}
est_ar = function(obj, p.max = NULL, penalty = NULL, ic = c('AIC','BIC','max'),
                  method = c('yule-walker', 'ols', 'durbin-levinson-whittle'),
                  mean_estimate = c('sample.mean', 'intercept','zero'), n.obs = NULL) {

  method = match.arg(method)
  ic = match.arg(ic)
  mean_estimate = match.arg(mean_estimate)

  if (inherits(obj, 'autocov')) {
    if (method == 'ols') {
      stop('for method="ols" the input parameter "obj" must be a matrix, time series or data-frame')
    }
    gamma = obj$gamma
    m = dim(gamma)[1] # output dimension
    lag.max = dim(gamma)[3] - 1
    if ( (m*(lag.max+1)) == 0 ) stop('"obj" contains no data')
    names = obj$names
    label = obj$label
    y.mean = rep(NA_real_, m)

    # check n.obs
    if (is.null(n.obs)) {
      n.obs = obj$n.obs
      if (is.null(n.obs)) n.obs = Inf
    }
    n.obs = as.numeric(n.obs)[1]
    if (is.finite(n.obs)) {
      n.obs = as.integer(n.obs)[1]
      if (n.obs <= 0) stop('the sample size "n.obs" must be non negative')
    } else {
      n.obs = Inf
    }

    if (is.null(p.max)) {
      p.max = max(0, min(12, lag.max, 10*log10(n.obs), floor((n.obs-1)/(m+1))))
    }
    p.max = as.integer(p.max)[1]
    if (p.max < 0) stop('p.max must be a non negative integer!')
  } else {
    y = try(as.matrix(obj))
    if ( inherits(y, 'try-error') || (!is.numeric(y)) || (!is.matrix(y)) ) {
      stop('input "obj" must be an "autocov" object or a "time series" ',
           'object which may be coerced to a matrix with "as.matrix(y)"')
    }
    m = ncol(y)      # number of "outputs"
    n.obs = nrow(y)  # sample size (optional parameter n.obs is ignored)
    if (m*n.obs == 0) stop('"obj" contains no data')

    names = colnames(y)
    label = NULL

    if (is.null(p.max)) {
      p.max = max(0, min(12, 10*log10(n.obs), floor((n.obs-1)/(m+1))))
    }
    p.max = as.integer(p.max)[1]
    if (p.max < 0) stop('p.max must be a non negative integer!')

    if (method != 'ols') {
      # compute ACF and mean
      demean = (mean_estimate != 'zero')
      gamma = autocov(y, lag.max = p.max, type = 'covariance',
                      na.action = stats::na.fail, demean = demean)$gamma
      if (demean) {
        y.mean = colMeans(y)
      } else {
        y.mean = double(m)
      }
    }
  }

  # set penalty
  if (is.null(penalty)) {
    if (ic == 'max') penalty = -1
    if (ic == 'AIC') penalty = 2/n.obs
    if (ic == 'BIC') penalty = log(n.obs)/n.obs
  }
  penalty = as.vector(penalty)[1]

  # now call the estimation methods
  if (method == 'ols') {
    out = est_ar_ols(y, p.max = p.max, mean_estimate = mean_estimate, penalty = penalty)
    y.mean = out$y.mean
  }
  if (method == 'yule-walker') {
    out = est_ar_yw(gamma, p.max = p.max, penalty = penalty)
  }
  if (method == 'durbin-levinson-whittle') {
    out = est_ar_dlw(gamma, p.max = p.max, penalty = penalty)
  }

  model = armamod(lmfd(a = dbind(d = 3, diag(m), -out$a)),
                  sigma_L = t(chol(out$sigma)), names = names, label = label)

  # log Likelihood
  p = out$p
  ll = unname(out$stats[p+1, 'lndetSigma'])
  if (is.finite(n.obs)) {
    ll = (-1 / 2) * (m*log(2*pi) + m + ll)
  } else {
    ll = NA_real_
  }

  return(list(model = model, p = p, stats = out$stats, y.mean = y.mean, ll = ll))
}


est_ar_yw = function(gamma, p.max = (dim(gamma)[3]-1), penalty = -1) {
  # no input check!

  m = dim(gamma)[1]   # number of "outputs"

  # G is the covariance matrix of (y[t-p.max]',y[t-p+1]',...,y[t]')'
  G = btoeplitz(C = gamma[,,1:(p.max+1),drop=FALSE])

  # cholesky decomposition of G
  R = chol(G)

  # vector of log likelihood values and information criteria
  # compute log(det(sigma_p)), where sigma_p is the noise covariance matrix of the AR(p) model
  # note that R[p+1,p+1] is the cholesky decomposition of sigma_p
  ldS = 2*apply(matrix(log(diag(R)), nrow = m, ncol = p.max+1), MARGIN = 2, FUN = sum)

  stats = matrix(NA_real_, nrow = p.max+1, 4)
  colnames(stats) = c('p', 'n.par', 'lndetSigma', 'ic')
  stats[, 'p'] = 0:p.max
  stats[, 'lndetSigma'] = ldS
  stats[, 'n.par'] = stats[,'p']*(m^2)
  stats[, 'ic'] = stats[, 'lndetSigma'] + stats[, 'n.par']*penalty

  # optimal order
  p = unname(which.min(stats[, 'ic']) - 1)

  # compute/estimate AR model of order p
  if (p > 0) {
    # AR coefficients a = (a[p],...,a[1]) are determined by R[1:p,1:p] t(a) = R[1:p,p+1]
    a = backsolve(R[1:(p*m), 1:(p*m), drop=FALSE], R[1:(m*p), (m*p+1):(m*(p+1)), drop = FALSE])
    # we have to reshuffle the coefficients
    a = t(a)
    dim(a) = c(m, m, p)
    a = a[ , , p:1, drop = FALSE]

    sigmaR = R[(p*m+1):((p+1)*m), (p*m+1):((p+1)*m), drop=FALSE]
    sigma = t(sigmaR) %*% sigmaR
  } else {
    a = array(0, dim = c(m, m, 0))
    sigma = matrix(gamma[,,1], nrow = m)
  }

  return(list(a = a, sigma = sigma, p = p, stats = stats))
}

est_ar_dlw = function(gamma, p.max = (dim(gamma)[3]-1), penalty = -1) {
  # no input checks!

  m = dim(gamma)[1]
  g0 = matrix(gamma[,,1], nrow = m, ncol = m)  # lag zero covariance

  # convert gamma to (m*p.max,m) matrix [gamma(1),...,gamma(p.max)]' = E [x[t-1]',...,x[t-p.max]']' x[t]'
  g.past = gamma[,,iseq(2, p.max+1),drop=FALSE]
  dim(g.past) = c(m, m*p.max)
  g.past = t(g.past)
  # print(g.past)

  # convert gamma to (m*p.max,m) matrix [gamma(1)',...,gamma(p.max)']' = E [x[t+1]',...,x[t+p.max]']' x[t]'
  # g.future = aperm(gamma[,,iseq(2,p.max+1),drop=FALSE],c(1,3,2))
  # dim(g.future) = c(m*p.max, m)
  # print(g.future)

  # initialize ( <=> order p=0 model)
  a = matrix(0, nrow = m, ncol = m*p.max) # matrix with the coefficients of the forward model
  a0 = a
  b = a            # matrix with the coefficients of the backward model
  sigma = g0       # covariance of the forecasting errors
  sigma0 = sigma
  omega = sigma    # covariance of the backcasting errors
  omegae0 = omega

  # partial autocorrelation coefficients
  partial = array(0, dim = c(m, m, p.max+1))
  su = matrix(sqrt(diag(sigma)), nrow = m, ncol = m)
  partial[,,1] = sigma / (su * t(su))

  stats = matrix(NA_real_, nrow = p.max+1, 4)
  colnames(stats) = c('p', 'n.par', 'lndetSigma', 'ic')
  stats[, 'p'] = 0:p.max
  stats[, 'n.par'] = stats[,'p']*(m^2)

  # optimal model so far
  stats[1, 'lndetSigma'] = log(det(sigma))
  stats[1, 'ic'] = stats[1, 'lndetSigma'] + stats[1, 'n.par']*penalty
  ic.opt = stats[1, 'ic']
  p.opt = 0
  a.opt = a
  sigma.opt = sigma

  # i is a matrix of indices, such that we can easily access the contents of g.future, g.past, a, b, ...
  # e.g. i[,1] = 1:m, i[,2] = (m+1):2*m, ...
  i = matrix(iseq(1,m*p.max), nrow = m, ncol = p.max)

  for (p in iseq(1, p.max)) {
    # compute order p model

    sigma0 = sigma
    omega0 = omega
    if (p > 1) a0[, i[,1:(p-1)]] = a[, i[,1:(p-1)]]

    # update forward model
    # E v[t-p] y[t]', where v[t-p] = y[t-p] - b[1] y[t-p+1] - ... - b[p-1] y[t-1]
    # is the backcasting error
    if (p > 1) {
      Evy = g.past[i[, p], , drop=FALSE] - b[,i[,1:(p-1)],drop=FALSE] %*% g.past[i[,(p-1):1],,drop=FALSE]
    } else {
      Evy = g.past[i[, p], , drop=FALSE]
    }
    aa = t(solve(omega, Evy))
    sigma = sigma - aa %*% Evy
    a[, i[, p]] = aa
    if (p>1) {
      a[, i[, 1:(p-1)]] = a[, i[, 1:(p-1)], drop=FALSE] - aa %*% b[, i[, (p-1):1], drop = FALSE]
    }

    # update backward model
    # E u[t] y[t-p]' = t( E v[t-p] y[t]' ), where u[t] = y[t] - a[1] y[t-1] - ... - a[p-1] y[t-p+1]
    # is the forecasting error
    # Euy = g.future[i[,p],,drop=FALSE] - a0[,i[,iseq(1,p-1)],drop=FALSE] %*%
    #       g.future[i[,rev(iseq(1,p-1))],,drop=FALSE]
    Euy = t(Evy)
    bb = t(solve(sigma0, Euy))
    omega = omega - bb %*% Euy
    b[, i[, p]] = bb
    if (p>1) {
      b[, i[, 1:(p-1)]] = b[, i[, 1:(p-1)], drop=FALSE] - bb %*% a0[, i[, (p-1):1], drop=FALSE]
    }

    # partial autocorrelations
    su = matrix(sqrt(diag(sigma0)), nrow=m, ncol=m)
    sv = matrix(sqrt(diag(omega0)), nrow=m, ncol=m,byrow = TRUE)
    # print(cbind(Euy,sigma0,omega0,su,sv))
    partial[,,p+1] = Euy / (su * sv)

    stats[p+1, 'lndetSigma'] = log(det(sigma))
    stats[p+1, 'ic'] = stats[p+1, 'lndetSigma'] + stats[p+1, 'n.par']*penalty

    # cat(p,'\n')
    # print(Evy)
    # print(omega)
    # print(eigen(omega)$values)
    # print(sigma)
    # print(eigen(sigma)$values)
    # print(ldS)
    # print(ic)
    if (stats[p+1, 'ic'] < ic.opt) {
      # update optimal model
      ic.opt = stats[p+1, 'ic']
      a.opt = a
      sigma.opt = sigma
      p.opt = p
    }
  }

  a = a.opt[, iseq(1,m*p.opt), drop = FALSE]
  dim(a) = c(m, m, p.opt)
  return(list(a = a, sigma = sigma.opt, p = p.opt, stats = stats, partial = partial))
}

est_ar_ols = function(y, p.max = NULL, penalty = -1,
                      mean_estimate = c('sample.mean', 'intercept','zero'), p.min = 0L) {
  # only some basic input checks
  mean_estimate = match.arg(mean_estimate)
  intercept = (mean_estimate == 'intercept')

  # coerce data objects to matrix
  y = try(as.matrix(y))
  if ( inherits(y, 'try-error') || (!is.numeric(y)) || (!is.matrix(y)) ) {
    stop('input "y" must be a data object which may be coerced to a matrix with "as.matrix(y)')
  }
  m = ncol(y)      # number of "outputs"
  n.obs = nrow(y)  # sample size
  if ( (m*n.obs == 0) ) stop('"y" contains no data')

  p.min = as.integer(p.min)[1]
  if (p.min < 0) stop('minimum order "p.min" must be a non negative integer')

  if (is.null(p.max)) {
    p.max = max(p.min, min(12, 10*log10(n.obs)/m, floor((n.obs-1)/(m+1)) ))
  }
  p.max = as.integer(p.max)[1]
  if ( (p.max < 0) ) stop('maximum order "p.max" must be a non negative integer!')
  if ( (p.max < p.min) ) stop('maximum order "p.max" is smaller than the minimum order "p.min"!')
  if ( (n.obs-p.max) < (m*p.max + intercept) ) {
    stop('sample size N is not sufficient for the desired maximum order "p.max"')
  }

  y.mean = double(m)
  if (mean_estimate == 'sample.mean') {
    y.mean = colMeans(y) # sample mean
    y = y - matrix(y.mean, nrow = n.obs, ncol = m, byrow = TRUE)
  }

  # create a "big" data matrix X, where the rows are of the form
  #    (y[t,], y[t-1,],...,y[t-p.max,])

  # use "btoeplitz"
  # dim(y) = c(n.obs,n,1)
  # y = aperm(y, c(3,2,1))
  # junk = array(0,dim=c(1,n,p.max+1))
  # junk[,,1] = y[,,1]
  # X = btoeplitz(R = junk, C = y)

  # use "for loop"
  X = matrix(NA_real_, nrow = n.obs, ncol = m*(p.max+1))
  for (i in (0:p.max)) {
    X[(1+i):n.obs, (m*i+1):(m*(i+1))] = y[1:(n.obs-i),]
  }
  # print(X[1:(n*(p.max+2)),])

  stats = matrix(NA_real_, nrow = p.max-p.min+1, 4)
  colnames(stats) = c('p', 'n.par', 'lndetSigma', 'ic')
  stats[, 'p'] = p.min:p.max
  ic.opt = Inf
  u.opt = matrix(NA_real_, nrow = n.obs, ncol = m)
  p.opt = NA_integer_

  # estimate AR models with order p = p.min, 1, ..., p.max
  for (p in (p.min:p.max)) {
    if (p == 0) {
      # AR(0) model
      a = matrix(0, nrow = m, ncol = 0)
      if (intercept) {
        d = colMeans(y)
        u = y - matrix(d, nrow = n.obs, ncol = m, byrow = TRUE)
      } else {
        d = double(m)
        u = y
      }
    } else {
      # estimate coefficients by OLS
      out = stats::lsfit(X[(1+p):n.obs, (m+1):(m*(p+1)), drop=FALSE],
                         X[(1+p):n.obs, 1:m, drop=FALSE], intercept = intercept)
      a = t(out$coef)
      if (intercept) {
        d = a[, 1]
        a = a[, -1, drop = FALSE]
      } else {
        d = NULL
      }
      u = out$residuals
    }
    sigma = crossprod(u)/(n.obs - p)
    # log(det(sigma)) may generate NaN's, suppress warning
    i = which(stats[,'p'] == p)
    stats[i, 'lndetSigma'] = suppressWarnings( log(det(sigma)) )
    stats[i, 'n.par'] = (m^2)*p + intercept*m
    stats[i, 'ic'] = stats[i, 'lndetSigma'] + stats[i, 'n.par']*penalty
    # take care of NaN's
    if ( is.finite(stats[i, 'ic']) && (stats[i, 'ic'] < ic.opt)) {
      # we have found a better model!
      a.opt = a
      d.opt = d
      sigma.opt = sigma
      p.opt = p
      ic.opt = stats[i, 'ic']
      u.opt[1:p, ] = NA_real_
      u.opt[(p+1):n.obs, ] = u
    }
  }

  # this should not happen
  if (is.na(p.opt)) stop('could not find an optimalt AR model')

  # coerce a.opt to 3-D array
  dim(a.opt) = c(m, m, p.opt)

  if (intercept) {
    # estimated mean from intercept mu = (I -a[1] - ... - a[p])^{-1} d
    if (p.opt == 0) {
      y.mean = d.opt
    } else {
      # a(1) = I - a[1] - ... - a[p]
      # print( apply(a.opt, MARGIN = c(1,2), FUN = sum) )
      a1 = diag(m) - matrix(apply(a.opt, MARGIN = c(1,2), FUN = sum), nrow = m, ncol = m)
      y.mean = try(solve(a1, d.opt))
      if (inherits(a1, 'try-error')) {
        # model has a unit root a(1) is singular!!!!
        y.mean = rep(NaN, m)
      }
    }
  }

  # # log likelihood
  # ll = unname((-(n.obs - p.opt)/2)*(m*log(2*pi) + m + stats[p.opt+1,'lndetSigma']))
  # # cat('est_ar_ols', (n.obs-p.opt)/2, m*log(2*pi), m, stats[p.opt+1, 'lndetSigma'], ll, '\n')

  # return the best model!
  return(list(a = a.opt, sigma = sigma.opt, p = p.opt,
              stats = stats, y.mean = y.mean, residuals = u.opt))
}

# freqresp_methods.R ######################################
# defines the main methods for the 'freqresp' class

dft_3D = function(a, n.f = dim(a)[3]) {
  n.f = as.integer(n.f)[1]
  if (n.f > 0) {
    z = exp((0:(n.f-1)) * (2*pi*complex(imaginary = -1)/n.f))
  } else {
    z = complex(0)
  }

  d = dim(a)
  if (min(c(d,n.f)) == 0) {
    a = array(complex(real = 0), dim = c(d[1], d[2], n.f))
    attr(a,'z') = z
    class(a) = c('zvalues','ratm')
    return(a)
  }

  dim(a) = c(d[1]*d[2], d[3])
  a = t(a)
  if (d[3] > n.f) {
    h = ceiling(d[3]/n.f)
    a = rbind(a, matrix(0, nrow = h*n.f - d[3], ncol = d[1]*d[2]))
    dim(a) = c(n.f, h, d[1]*d[2])
    a = apply(a, MARGIN = c(1,3), FUN = sum)
  }
  if (d[3] < n.f) {
    a = rbind(a, matrix(0, nrow = n.f - d[3], ncol = d[1]*d[2]))
  }
  # print(dim(a))
  a = stats::mvfft(a)
  # print(dim(a))
  a = t(a)
  dim(a) = c(d[1],d[2],n.f)
  attr(a,'z') = exp((0:(n.f-1)) * (2*pi*complex(imaginary = -1)/n.f))
  class(a) = c('zvalues','ratm')
  return(a)
}


freqresp = function(obj, n.f, ...) {
  UseMethod("freqresp", obj)
}

freqresp.armamod = function(obj, n.f = 128, ...) {
  n.f = as.integer(n.f)[1]
  if (n.f < 0) stop('the number of frequencies "n.f" must be a non negative integer')

  d = dim(obj$sys)
  m = d[1]
  n = d[2]
  p = d[3]
  q = d[4]
  if (n.f > 0) {
    z = exp((0:(n.f-1)) * (2*pi*complex(imaginary = -1)/n.f))
  } else {
    z = complex(0)
  }


  if ( min(c(m,n,q+1,n.f)) == 0) {
    frr = array(complex(real = 0), dim = c(m, n, n.f))
    attr(frr, 'z') = complex(0)
    class(frr) = c('zvalues', 'ratm')

    out = structure(list(frr = frr, sigma_L = obj$sigma_L, names = obj$names, label = obj$label),
                    class = c('freqresp', 'rldm'))
  }
  a = unclass(dft_3D(unclass(obj$sys$a), n.f))
  frr = unclass(dft_3D(unclass(obj$sys$b), n.f))
  if (m == 1) {
    a = aperm(array(a, dim = c(m, n.f, n)), c(1, 3, 2))
    frr = frr / a
  } else {
    for (i in (1:n.f)) {
      b = try(solve(a[,,i], frr[,,i]), silent = TRUE)
      if (!inherits(b, 'try-error')) {
        frr[,,i] = b
      } else {
        frr[,,i] = NA_complex_
      }
    }
  }
  attr(frr, 'z') = z
  class(frr) = c('zvalues', 'ratm')

  out = structure(list(frr = frr, sigma_L = obj$sigma_L, names = obj$names, label = obj$label),
                  class = c('freqresp', 'rldm'))
  return(out)
}


freqresp.stspmod = function(obj, n.f = 128, ...) {
  n.f = as.integer(n.f)[1]
  if (n.f < 0) stop('the number of frequencies "n.f" must be a non negative integer')

  frr = zvalues(obj$sys, n.f = n.f)
  # make sure that the frr object is of type 'complex'
  if (!is.complex(frr)) {
    frr = unclass(frr)
    z = attr(frr, 'z')
    d = dim(frr)
    frr = array(as.complex(frr), dim = d)
    frr = structure(frr, z = z, class = c('zvalues', 'ratm'))
  }

  out = structure(list(frr = frr, sigma_L = obj$sigma_L, names = obj$names, label = obj$label),
                  class = c('freqresp', 'rldm'))
  return(out)
}




freqresp.impresp = function(obj, n.f = 128, ...) {
  n.f = as.integer(n.f)[1]
  if (n.f < 0) stop('the number of frequencies "n.f" must be a non negative integer')

  irf = unclass(obj$irf)
  frr = dft_3D(irf, n.f)

  out = structure(list(frr = frr, sigma_L = obj$sigma_L, names = obj$names, label = obj$label),
                  class = c('freqresp', 'rldm'))
  return(out)
}


# impresp_methods.R ######################################
# defines the main methods for the 'impresp' class

impresp = function(obj, lag.max, H) {
  UseMethod("impresp", obj)
}

# internal helper function which computes the transformation matrix "H"
make_H = function(type = c('chol','eigen','sigma_L'), sigma_L) {
  n = nrow(sigma_L)
  type = match.arg(type)

  # H = sigma_L
  if ( type == 'sigma_L' ) {
    H = sigma_L
  }

  # cholesky decomposition
  if ( type == 'chol' ) {
    sigma = sigma_L %*% t(sigma_L)
    H = t(chol(sigma))
  }

  # eigenvalue decomposition
  if ( type == 'eigen' ) {
    sigma = sigma_L %*% t(sigma_L)
    ed = eigen(sigma, symmetric = TRUE)
    H = ed$vectors %*% diag(x = sqrt(ed$values), nrow = n) %*% t(ed$vectors)
  }

  return(H)
}


impresp.armamod = function(obj, lag.max = 12, H = NULL) {

  lag.max = as.integer(lag.max)[1]
  if (lag.max < 0) {
    stop('lag.max must be a non-negative integer.')
  }

  sys = obj$sys
  sigma_L = obj$sigma_L

  irf = pseries(sys, lag.max = lag.max)

  # Compute orthogonalized IRF #
  if ((!is.null(H)) && (nrow(sigma_L) > 0)) {
    if (is.character(H)) {
      H = make_H(H, sigma_L)
    }
    if ( (!is.numeric(H)) || (!is.matrix(H)) || (any(dim(H) != nrow(sigma_L))) ) {
      stop('H must be a numeric matrix of the same dimension as the noise.')
    }
    irf = irf %r% H
    sigma_L = solve(H, sigma_L)
  }

  out = structure(list(irf = irf, sigma_L = sigma_L, names = obj$names, label = obj$label),
                  class = c('impresp','rldm'))

  return(out)
}

impresp.rmfdmod = function(obj, lag.max = 12, H = NULL) {

  lag.max = as.integer(lag.max)[1]
  if (lag.max < 0) {
    stop('lag.max must be a non-negative integer.')
  }

  sys = obj$sys
  sigma_L = obj$sigma_L

  irf = pseries(sys, lag.max = lag.max)

  # Compute orthogonalized IRF #
  if ((!is.null(H)) && (nrow(sigma_L) > 0)) {
    if (is.character(H)) {
      H = make_H(H, sigma_L)
    }
    if ( (!is.numeric(H)) || (!is.matrix(H)) || (any(dim(H) != nrow(sigma_L))) ) {
      stop('H must be a numeric matrix of the same dimension as the noise.')
    }
    irf = irf %r% H
    sigma_L = solve(H, sigma_L)
  }

  out = structure(list(irf = irf, sigma_L = sigma_L, names = obj$names, label = obj$label),
                  class = c('impresp','rldm'))

  return(out)
}

impresp.stspmod = function(obj, lag.max = 12, H = NULL) {

  # code is completey identical to the code of 'impresp.armamod #

  lag.max = as.integer(lag.max)[1]
  if (lag.max < 0) {
    stop('lag.max must be a non-negative integer.')
  }

  sys = obj$sys
  sigma_L = obj$sigma_L

  irf = pseries(sys, lag.max = lag.max)

  # Compute orthogonalized IRF #
  if ((!is.null(H)) && (nrow(sigma_L) > 0)) {
    if (is.character(H)) {
      H = make_H(H, sigma_L)
    }
    if ( (!is.numeric(H)) || (!is.matrix(H)) || (any(dim(H) != nrow(sigma_L))) ) {
      stop('H must be a numeric matrix of the same dimension as the noise.')
    }
    irf = irf %r% H
    sigma_L = solve(H, sigma_L)
  }

  out = structure(list(irf = irf, sigma_L = sigma_L, names = obj$names, label = obj$label),
                  class = c('impresp','rldm'))

  return(out)
}


impresp.impresp = function(obj, lag.max = NULL, H = NULL) {
  # code is almost identical to the code of 'impresp.armamod #

  sigma_L = obj$sigma_L
  irf = obj$irf

  # Compute orthogonalized IRF #
  if ((!is.null(H)) && (nrow(sigma_L) > 0)) {
    if (is.character(H)) {
      H = make_H(H, sigma_L)
    }
    if ( (!is.numeric(H)) || (!is.matrix(H)) || (any(dim(H) != nrow(sigma_L))) ) {
      stop('H must be a numeric matrix of the same dimension as the noise.')
    }
    irf = irf %r% H
    sigma_L = solve(H, sigma_L)
  }

  out = structure(list(irf = irf, sigma_L = sigma_L, names = obj$names, label = obj$label),
                  class = c('impresp','rldm'))

  return(out)
}

# poles.___ and zeroes.___ method ############################################################
#
#
zeroes.armamod = function(x, tol = sqrt(.Machine$double.eps), print_message = TRUE, ...) {
  z = zeroes(x$sys, tol = tol, print_message = print_message)
  return(z)
}


poles.armamod = function(x, tol = sqrt(.Machine$double.eps), print_message = TRUE, ...) {
  z = poles(x$sys, tol = tol, print_message = print_message)
  return(z)
}

zeroes.rmfdmod = function(x, tol = sqrt(.Machine$double.eps), print_message = TRUE, ...) {
  z = zeroes(x$sys, tol = tol, print_message = print_message)
  return(z)
}


poles.rmfdmod = function(x, tol = sqrt(.Machine$double.eps), print_message = TRUE, ...) {
  z = poles(x$sys, tol = tol, print_message = print_message)
  return(z)
}

zeroes.stspmod = function(x, tol = sqrt(.Machine$double.eps), print_message = TRUE, ...) {
  z = zeroes(x$sys, tol = tol, print_message = print_message)
  return(z)
}


poles.stspmod = function(x, tol = sqrt(.Machine$double.eps), print_message = TRUE, ...) {
  z = poles(x$sys, tol = tol, print_message = print_message)
  return(z)
}
# RLDM - print.____() methods ##############################################################

NULL

print.armamod = function(x, digits = NULL,
                         format = c('i|jz', 'i|zj', 'iz|j', 'zi|j', 'i|j|z', 'character'), ...) {
  if (!is.null(digits)) {
    digits = as.vector(as.integer(digits))[1]
  }

  format = match.arg(format)
  if ((format == 'character') && (is.complex(unclass(x$sys)))) {
    stop('the format option "character" is only implemented for ARMA models with real coefficients.')
  }

  label = ''
  if (!is.null(x$label)) label = paste(x$label, ': ', sep = '')

  d = attr(x$sys, 'order')
  m = d[1]
  n = d[2]
  p = d[3]
  q = d[4]

  cat(label, 'ARMA model [',d[1],',',d[2],'] with orders p = ', d[3], ' and q = ', d[4], '\n', sep = '')

  if ((m*m*(p+1)) > 0) {
    cat('AR polynomial a(z):\n')

    a = unclass(x$sys$a)

    # use the function print_3D (contained in rationalmatrices)
    dimnames(a) = list(paste('[', 1:m, ',]', sep = ''),
                       paste('[,', 1:m, ']', sep = ''),
                       paste('z^',0:p, sep = ''))
    print_3D(a, digits, format)
  }

  if ((m*n*(q+1)) > 0) {
    cat('MA polynomial b(z):\n')

    a = unclass(x$sys$b)

    # use the above defined internal function print_3D
    dimnames(a) = list(paste('[', 1:m, ',]', sep = ''),
                       paste('[,', 1:n, ']', sep = ''),
                       paste('z^',0:q, sep = ''))
    print_3D(a, digits, format)
  }

  if (n > 0) {
    cat('Left square root of noise covariance Sigma:\n')

    a = x$sigma_L
    dimnames(a) = list(paste('u[',1:n,']',sep = ''),paste('u[',1:n,']',sep = ''))
    print(a)
  }

  return(invisible(x))
}

print.rmfdmod = function(x, digits = NULL,
                         format = c('i|jz', 'i|zj', 'iz|j', 'zi|j', 'i|j|z', 'character'), ...) {
  if (!is.null(digits)) {
    digits = as.vector(as.integer(digits))[1]
  }

  format = match.arg(format)
  if ((format == 'character') && (is.complex(unclass(x$sys)))) {
    stop('the format option "character" is only implemented for ARMA models with real coefficients.')
  }

  label = ''
  if (!is.null(x$label)) label = paste(x$label, ': ', sep = '')

  d = attr(x$sys, 'order')
  m = d[1]
  n = d[2]
  p = d[3]
  q = d[4]

  cat(label, 'RMFD model [',d[1],',',d[2],'] with orders p = ', d[3], ' and q = ', d[4], '\n', sep = '')

  if ((n*n*(p+1)) > 0) {
    cat('right factor polynomial c(z):\n')

    c = unclass(x$sys$c)

    # use the above defined internal function print_3D
    dimnames(c) = list(paste('[', 1:n, ',]', sep = ''),
                       paste('[,', 1:n, ']', sep = ''),
                       paste('z^',0:p, sep = ''))
    print_3D(c, digits, format)
  }

  if ((m*n*(q+1)) > 0) {
    cat('left factor polynomial d(z):\n')

    a = unclass(x$sys$d)

    # use the above defined internal function print_3D
    dimnames(a) = list(paste('[', 1:m, ',]', sep = ''),
                       paste('[,', 1:n, ']', sep = ''),
                       paste('z^',0:q, sep = ''))
    print_3D(a, digits, format)
  }

  if (n > 0) {
    cat('Left square root of noise covariance Sigma:\n')

    a = x$sigma_L
    dimnames(a) = list(paste('u[',1:n,']',sep = ''),paste('u[',1:n,']',sep = ''))
    print(a)
  }

  return(invisible(x))
}

print.stspmod = function(x, digits = NULL, ...) {
  if (!is.null(digits)) {
    digits = as.vector(as.integer(digits))[1]
  }

  d = attr(x$sys, 'order')
  m = d[1]
  n = d[2]
  s = d[3]

  label = ''
  if (!is.null(x$label)) label = paste(x$label, ': ', sep = '')

  cat(label, 'Statespace model [', m, ',', n, '] with s = ', s, ' states\n', sep = '')

  a = unclass(x$sys)
  attr(a, 'order') = NULL
  if (length(a) == 0) {
    return(invisible(x))
  }

  # rounding digits
  if (!is.null(digits)) {
    a = round(a, digits)
  }

  snames = character(s)
  if (s > 0) snames = paste('s[',1:s,']',sep = '')
  xnames = character(m)
  if (m > 0) {
    if ( !is.null(x$names) && is.character(x$names) && is.vector(x$names) && (length(x$names) == m) ) {
      xnames = x$names
    } else {
      xnames = paste('x[',1:m,']',sep = '')
    }
  }
  unames = character(n)
  if (n > 0) unames = paste('u[',1:n,']',sep = '')

  rownames(a) = c(snames, xnames)
  colnames(a) = c(snames, unames)
  print(a)

  if (n > 0) {
    cat('Left square root of noise covariance Sigma:\n')

    a = x$sigma_L
    dimnames(a) = list(paste('u[',1:n,']',sep = ''),paste('u[',1:n,']',sep = ''))
    print(a)
  }

  invisible(x)
}

print.impresp = function(x, digits = NULL,
                         format = c('i|jz', 'i|zj', 'iz|j', 'zi|j', 'i|j|z'), ...) {
  if (!is.null(digits)) {
    digits = as.vector(as.integer(digits))[1]
  }
  format = match.arg(format)

  d = dim(x$irf)
  m = d[1]
  n = d[2]
  lag.max = d[3]

  label = ''
  if (!is.null(x$label)) label = paste(x$label, ': ', sep = '')

  orth = FALSE
  if (n > 0) {
    sigma = x$sigma_L
    sigma = sigma %*% t(sigma)
    orth = isTRUE(all.equal(x$sigma_L, diag(n)))
  }

  if (orth) {
    cat(label, 'Orthogonalized impulse response [', m, ',', n, '] with ', lag.max, ' lags\n', sep = '')
  } else {
    cat(label, 'Impulse response [', m, ',', n, '] with ', lag.max, ' lags\n', sep = '')
  }

  if ((m*n*(lag.max+1)) > 0) {
    a = unclass(x$irf)

    # use the above defined internal function print_3D
    dimnames(a) = list(paste('[', 1:m, ',]', sep = ''),
                       paste('[,', 1:n, ']', sep = ''),
                       paste('lag=', 0:lag.max, sep = ''))
    print_3D(a, digits, format)
  }

  invisible(x)
}

print.autocov = function(x, digits = NULL,
                         format = c('i|jz', 'i|zj', 'iz|j', 'zi|j', 'i|j|z'), ...) {
  if (!is.null(digits)) {
    digits = as.vector(as.integer(digits))[1]
  }
  format = match.arg(format)

  d = dim(x$acf)
  m = d[1]
  lag.max = d[3]

  label = ''
  if (!is.null(x$label)) label = paste(x$label, ': ', sep = '')

  if (x$type == 'covariance') {
    label = paste(label, 'Autocovariance function [',d[1],',',d[2],'] with ', d[3], ' lags', sep = '')
  }
  if (x$type == 'correlation') {
    label = paste(label, 'Autocorrelation function [',d[1],',',d[2],'] with ', d[3], ' lags', sep = '')
  }
  if (x$type == 'partial') {
    label = paste(label, 'Partial autocorrelation function [',d[1],',',d[2],'] with ', d[3], ' lags', sep = '')
  }

  n.obs = x$n.obs
  if (!is.null(n.obs)) {
    n.obs = as.integer(n.obs)[1]
  } else {
    n.obs = Inf
  }
  if ( is.finite(n.obs) ) {
    label = paste(label, ', sample size is ', n.obs, sep = '')
  }

  cat(label, '\n')

  if ((m*(lag.max+1)) > 0) {
    a = unclass(x$acf)

    # use the above defined internal function print_3D
    dimnames(a) = list(paste('[', 1:m, ',]', sep = ''),
                       paste('[,', 1:m, ']', sep = ''),
                       paste('lag=', 0:lag.max, sep = ''))
    print_3D(a, digits, format)
  }
  invisible(x)
}

print.fevardec = function(x, digits = NULL,
                          format = c('i|jz', 'i|zj', 'iz|j', 'zi|j', 'i|j|z'), ...) {
  if (!is.null(digits)) {
    digits = as.vector(as.integer(digits))[1]
  }
  format = match.arg(format)

  vd = x$vd
  d = dim(vd)
  n = d[1]
  h.max = d[3]

  label = ''
  if (!is.null(x$label)) label = paste(x$label, ': ', sep = '')
  cat(label, 'Forecast error variance decompositon [', n, ',', n,'] maximum horizon = ', h.max, '\n', sep = '')

  if ((n*h.max)> 0) {
    names = x$names
    if (is.null(names)) {
      names = paste('y[', 1:n, ']', sep = '')
    }
    unames = paste('u[', 1:n, ']', sep = '')
    hnames = paste('h=', 1:h.max, sep='')

    dimnames(vd) = list(names, unames, hnames)
    print_3D(vd, digits, format)
  }

  invisible(x)
}

print.freqresp = function(x, digits = NULL,
                          format = c('i|jz', 'i|zj', 'iz|j', 'zi|j', 'i|j|z'), ...) {
  if (!is.null(digits)) {
    digits = as.vector(as.integer(digits))[1]
  }
  format = match.arg(format)

  d = dim(x$frr)
  m = d[1]
  n = d[2]
  n.f = d[3]

  label = ''
  if (!is.null(x$label)) label = paste(x$label, ': ', sep = '')

  cat(label, 'Frequency response [',d[1],',',d[2],'] with ', d[3], ' frequencies\n', sep = '')

  if ((m*n*n.f) > 0) {
    a = unclass(x$frr)
    attr(a, 'z') = NULL
    f = (0:(n.f-1))/n.f

    # use the above defined internal function print_3D
    if ((format == 'i|jz') || (format == 'i|zj')) {
      dimnames(a) = list(paste('[', 1:m, ',]', sep = ''),
                         paste('[,', 1:n, ']', sep = ''),
                         paste('f[',1:n.f, ']', sep = ''))
    } else {
      dimnames(a) = list(paste('[', 1:m, ',]', sep = ''),
                         paste('[,', 1:n, ']', sep = ''),
                         paste('f=', round(f,3), sep = ''))
    }
    print_3D(a, digits, format)
  }

  invisible(x)
}

print.spectrald = function(x, digits = NULL,
                           format = c('i|jz', 'i|zj', 'iz|j', 'zi|j', 'i|j|z'), ...) {
  if (!is.null(digits)) {
    digits = as.vector(as.integer(digits))[1]
  }
  format = match.arg(format)

  d = dim(x$spd)
  m = d[1]
  n = d[2]
  n.f = d[3]

  label = ''
  if (!is.null(x$label)) label = paste(x$label, ': ', sep = '')

  label = paste(label, 'Frequency response [',d[1],',',d[2],'] with ', d[3], ' frequencies', sep = '')

  n.obs = x$n.obs
  if (!is.null(n.obs)) {
    n.obs = as.integer(n.obs)[1]
  } else {
    n.obs = Inf
  }
  if ( is.finite(n.obs) ) {
    label = paste(label, ', sample size is ', n.obs, sep = '')
  }

  cat(label, '\n')

  if ((m*n.f) > 0) {
    a = unclass(x$spd)
    attr(a, 'z') = NULL
    f = (0:(n.f-1))/n.f

    # use the above defined internal function print_3D
    if ((format == 'i|jz') || (format == 'i|zj')) {
      dimnames(a) = list(paste('[', 1:m, ',]', sep = ''),
                         paste('[,', 1:n, ']', sep = ''),
                         paste('f[',1:n.f, ']', sep = ''))
    } else {
      dimnames(a) = list(paste('[', 1:m, ',]', sep = ''),
                         paste('[,', 1:n, ']', sep = ''),
                         paste('f=', round(f,3), sep = ''))
    }
    print_3D(a, digits, format)
  }

  invisible(x)
}

# spectrald_methods.R ######################################
# defines the main methods for the 'spectrald' class


spectrald = function(obj, n.f, ...) {
  UseMethod("spectrald", obj)
}

spectrald.armamod = function(obj, n.f = 128, ...) {
  n.f = as.integer(n.f)[1]
  if (n.f < 0) stop('the number of frequencies "n.f" must be a non negative integer')

  out = unclass(freqresp(obj, n.f))
  frr = out$frr
  out$frr = NULL
  frr = frr %r% (out$sigma_L / sqrt(2*pi))
  out$spd = frr %r% Ht(frr)

  out = out[c('spd','names','label')]
  class(out) = c('spectrald','rldm')
  return(out)
}

spectrald.stspmod = function(obj, n.f = 128, ...) {
  n.f = as.integer(n.f)[1]
  if (n.f < 0) stop('the number of frequencies "n.f" must be a non negative integer')

  out = unclass(freqresp(obj, n.f))
  frr = out$frr
  out$frr = NULL
  frr = frr %r% (out$sigma_L / sqrt(2*pi))
  out$spd = frr %r% Ht(frr)

  out = out[c('spd','names','label')]
  class(out) = c('spectrald','rldm')
  return(out)
}


spectrald.autocov = function(obj, n.f = 128, ...) {
  n.f = as.integer(n.f)[1]
  if (n.f < 0) stop('the number of frequencies "n.f" must be a non negative integer')

  gamma = obj$gamma
  gamma[,,1] = gamma[,,1] / 2
  gamma = gamma/(2*pi)

  spd = dft_3D(gamma, n.f)
  spd = spd + Ht(spd)
  out = list(spd = spd, names = obj$names, label = obj$label)
  if (!is.null(obj$n.obs)) out$n.obs = obj$n.obs
  class(out) = c('spectrald', 'rldm')
  return(out)
}


spectrald.impresp = function(obj, n.f = 128, ...) {
  n.f = as.integer(n.f)[1]
  if (n.f < 0) stop('the number of frequencies "n.f" must be a non negative integer')

  out = unclass(freqresp(obj, n.f))
  frr = out$frr
  out$frr = NULL
  frr = frr %r% (out$sigma_L / sqrt(2*pi))
  out$spd = frr %r% Ht(frr)

  out = out[c('spd','names','label')]
  class(out) = c('spectrald','rldm')
  return(out)
}


spectrald.default = function(obj, n.f = NULL, demean = TRUE, ...) {
  y = try(as.matrix(obj))
  if (inherits(y, 'try-error')) stop('could not coerce input parameter "obj" to matrix')

  n.obs = nrow(y)
  m = ncol(y)

  if (is.null(n.f)) n.f = n.obs
  n.f = as.integer(n.f)[1]
  if (n.f < 0) stop('the number of frequencies "n.f" must be a non negative integer')

  if (demean) {
    y.mean = apply(y, MARGIN = 2, FUN = mean)
    y = y - matrix(y.mean, nrow = n.obs, ncol = m, byrow = TRUE)
  }
  dim(y) = c(n.obs, 1, m)
  y = aperm(y, c(3,2,1))

  spd = dft_3D(y / sqrt(2*pi*n.obs), n.f)
  spd = spd %r% Ht(spd)

  out = list(spd = spd, names = colnames(y), label = NULL, n.obs = n.obs)
  class(out) = c('spectrald', 'rldm')

  return(out)
}

# str.____ methods ##############################################################

str.armamod = function(object, ...) {
  d = attr(object$sys, 'order')
  cat('ARMA model [',d[1],',',d[2],'] with orders p = ', d[3], ' and q = ', d[4], '\n', sep = '')
  return(invisible(NULL))
}

str.stspmod = function(object, ...) {
  d = attr(object$sys, 'order')
  cat('Statespace model [',d[1],',',d[2],'] with s = ', d[3], ' states\n', sep = '')
  return(invisible(NULL))
}

str.impresp = function(object, ...) {
  d = dim(object$irf)
  orth = FALSE
  if (d[2] > 0) {
    sigma = object$sigma_L
    sigma = sigma %*% t(sigma)
    orth = isTRUE(all.equal(sigma, diag(d[2])))
  }
  if (orth) {
    cat('Orthogonalized impulse response [',d[1],',',d[2],'] with ', d[3], ' lags\n', sep = '')
  } else {
    cat('Impulse response [',d[1],',',d[2],'] with ', d[3], ' lags\n', sep = '')
  }
  return(invisible(NULL))
}

str.autocov = function(object, ...) {
  d = dim(object$acf)
  type = object$type
  if (type == 'covariance') {
    cat('Autocovariance function [',d[1],',',d[2],'] with ', d[3], ' lags\n', sep = '')
  }
  if (type == 'correlation') {
    cat('Autocorrelation function [',d[1],',',d[2],'] with ', d[3], ' lags\n', sep = '')
  }
  if (type == 'partial') {
    cat('Partial autocorrelation function [',d[1],',',d[2],'] with ', d[3], ' lags\n', sep = '')
  }
  return(invisible(NULL))
}

str.fevardec = function(object, ...) {
  d = dim(object$vd)
  cat('Forecast error variance decompositon [',d[1],',',d[2],'] maximum horizon = ', d[3], '\n', sep = '')
  return(invisible(NULL))
}


str.freqresp = function(object, ...) {
  d = dim(object$frr)
  cat('Frequency response [',d[1],',',d[2],'] with ', d[3], ' frequencies\n', sep = '')
  return(invisible(NULL))
}


str.spectrald = function(object, ...) {
  d = dim(object$spd)
  cat('Spectral density [',d[1],',',d[2],'] with ', d[3], ' frequencies\n', sep = '')
  return(invisible(NULL))
}


# subspace estimates #############

estorder_SVC = function(s.max, Hsv=NULL, n.par, n.obs, Hsize, penalty = 'lnN', ...) {
  # cat('estorder_SVC start\n')
  # print(penalty)
  if (is.null(Hsv)) return(NULL)
  if (isTRUE(all.equal(penalty, 'lnN'))) penalty = log(n.obs)
  if (isTRUE(all.equal(penalty, 'fplnN'))) penalty = prod(Hsize)*log(n.obs)
  if (is.character(penalty)) stop('unknown option penalty=', penalty)
  penalty = as.numeric(penalty)[1]
  # print(Hsv)
  # print(n.par)
  # print(n.obs)
  # print(penalty)
  penalty = penalty / n.obs
  # take care of the case n.obs = Inf: log(Inf)/Inf = NaN
  if (!is.finite(penalty)) penalty = 0

  # Hsv should have at least s.max entries!
  if (length(Hsv) < (s.max+1)) Hsv = c(Hsv, 0)
  criterion = Hsv[1:(s.max+1)]^2 + n.par*penalty
  s = which.min(criterion) - 1
  # print(criterion)
  # print(s)
  # cat('estorder_SVC end\n')
  return(list(s = s, criterion = criterion))
}

estorder_IVC = function(s.max, lndetSigma=NULL, n.par, n.obs, penalty = 'BIC', ...) {
  # cat('estorder_IVC start\n')
  # print(penalty)
  if (is.null(lndetSigma)) return(NULL)
  if (isTRUE(all.equal(penalty, 'AIC'))) penalty = 2
  if (isTRUE(all.equal(penalty, 'BIC'))) penalty = log(n.obs)
  if (is.character(penalty)) stop('unknown option penalty=', penalty)
  penalty = as.numeric(penalty)[1]
  # print(lndetSigma)
  # print(n.par)
  # print(n.obs)
  # print(penalty)
  penalty = penalty / n.obs
  # take care of the case n.obs = Inf: log(Inf)/Inf = NaN
  if (!is.finite(penalty)) penalty = 0

  criterion = lndetSigma + n.par*penalty
  s = which.min(criterion) - 1
  # print(criterion)
  # print(s)
  # cat('estorder_IVC end\n')
  return(list(s = s, criterion = criterion))
}

estorder_max = function(s.max, ...) {
  return(list(s = s.max, criterion = c(rep(1, s.max), 0)))
}

estorder_rkH = function(s.max, Hsv = NULL, tol = sqrt(.Machine$double.eps), ...) {
  if (is.null(Hsv)) return(NULL)
  if (Hsv[1] <= .Machine$double.eps) {
    return(list(s = 0, criterion = numeric(s.max+1)))
  }
  criterion = as.integer(c(Hsv[1:s.max],0) > tol * Hsv[1])
  s = sum(criterion)
  return(list(s = s, criterion = criterion))
}

# MOE(n), which is implemented in the N4SID procedure of the system
# identifcation toolbox of MATLAB (Ljung, 1991):
# The idea here is to formalise the search for a "gap" in
# the singular values.

estorder_MOE = function(s.max, Hsv = NULL, ...) {
  if (is.null(Hsv)) return(NULL)
  criterion = Hsv
  if (Hsv[1] <= .Machine$double.eps) {
    return(list(s = 0, criterion = numeric(s.max+1)))
  }
  criterion = as.integer(log(Hsv) > (log(Hsv[1])+log(Hsv[length(Hsv)]))/2)
  criterion = c(criterion[1:s.max], 0)
  s = sum(criterion)
  return(list(s = s, criterion = criterion))
}


# core computations of the AOKI method #####################
aoki = function(m, s, svd.H, gamma0, Rp, Rf) {
  if (s == 0) {
    return(list(s = 0, A = matrix(0, nrow = 0, ncol = 0), B = matrix(0, nrow = 0, ncol = m),
                C = matrix(0, nrow = m, ncol = 0), D = diag(m),
                sigma = gamma0, sigma_L = t(chol(gamma0))))
  }

  sv2 = sqrt(svd.H$d[1:s])

  # (approximately) factorize H as H = U V
  U = svd.H$u[, 1:s, drop = FALSE] * matrix(sv2, nrow = nrow(svd.H$u), ncol = s, byrow = TRUE)
  V = t( svd.H$v[, 1:s, drop = FALSE] * matrix(sv2, nrow = nrow(svd.H$v), ncol = s, byrow = TRUE) )
  U = t(Rf) %*% U
  V = V %*% Rp

  C = U[1:m, , drop = FALSE]
  M = V[, 1:m, drop = FALSE]
  A = suppressWarnings(stats::lsfit(U[1:(nrow(U) - m), , drop = FALSE],
                                    U[(m+1):nrow(U), , drop = FALSE],
                                    intercept = FALSE)$coef)
  A = unname(A)
  A = matrix(A, nrow = s, ncol = s)

  # solve the Riccati equation to get the "B" matrix (Kalman gain)
  out = suppressWarnings(riccati(A, M, C, G = gamma0, only.X = FALSE))

  # construct state space model object
  sigma = out$sigma
  sigma_L = t(chol(sigma))
  return(list(s = s, A = A, B = out$B, C = C, D = diag(m), sigma = sigma, sigma_L = sigma_L))
}

est_stsp_aoki = function(gamma, s.max, p, estorder = estorder_SVC,
                         keep_models = FALSE, n.obs = NULL, ...) {

  # check gamma
  if ( (!is.numeric(gamma)) || (!is.array(gamma)) ||
       (length(dim(gamma)) != 3) || (dim(gamma)[1] != dim(gamma)[2]) ) {
    stop('input "gamma" is not a valid 3-D array.')
  }
  d = dim(gamma)
  m = d[1]
  if (m == 0) {
    stop('the autocovariance function "gamma" is empty.')
  }
  lag.max = d[3] - 1
  if (lag.max <= 2) {
    stop('the autocovariance function must contain at least 2 lags.',
         ' lag.max=', lag.max)
  }
  if (is.null(n.obs)) {
    n.obs = Inf
  }
  n.obs = as.numeric(n.obs)[1]
  if (is.finite(n.obs)) {
    n.obs = as.integer(n.obs)[1]
    if (n.obs <= 0) stop('the sample size "n.obs" must be non negative')
  } else {
    n.obs = Inf
  }

  # p <=> past
  p = as.integer(p)[1]
  if (p < 1) stop('the number of lags "p" must be positive')
  if (lag.max < (2*p)) {
    stop('the autocovariance function must contain at least 2p lags.',
         ' lag.max=', lag.max, ' < 2*p=', 2*p)
  }
  # f <=> future
  f = p + 1

  s.max = as.integer(s.max)[1]
  if (s.max < 0) stop('maximum state dimension "s.max" must be non negative.')
  if (s.max > (p*m)) {
    stop('the number of lags "p" is too small for the desired maximum state dimension:',
         ' p*m=', p*m, ' < s.max=', s.max)
  }

  gamma0 = matrix(gamma[,,1], nrow = m, ncol = m)

  # covariance between "future" (y[t]', y[t+1]',...,y[t+f-1]')' and
  #                    "past" (y[t-1]', y[t-2]',...,y[t-p]')'
  H = bhankel(gamma[,,-1,drop=FALSE], d = c(f,p))
  # H0 = H
  # junk = svd(H, nu = 0, nv = 0)$d
  # print(junk / junk[1])

  # covariance of the "past" (y[t-1]', y[t-2]',...,y[t-p]')'
  Gp = btoeplitz(R = gamma[,,1:p,drop = FALSE])
  # covariance matrix of the "future" (y[t]', y[t+1]',...,y[t+f-1]')'
  Gf = btoeplitz(C = gamma[,,1:f,drop = FALSE])

  # cholesky factors of covariance matrices
  Rp = chol(Gp)
  Rf = chol(Gf)

  # weighted Hankel matrix: inv(t(Rf)) * H * inv(Rp)
  # note: inv(t(Rf)) * H * inv(Rp) is the covariance between the
  # "standardized" future and the "standardized" past.
  # The singular value decomposition of this matrix defines the
  # canonical correlation coefficients between future and past.
  ## junk = H
  H = backsolve(Rf, H, transpose = TRUE)
  ## testthat::expect_equivalent(junk, t(Rf) %*% H)
  ## junk = H
  H = t(backsolve(Rp, t(H), transpose = TRUE))
  ## testthat::expect_equivalent(junk, H %*% Rp)

  # compute SVD of weighted Hankel matrix
  svd.H = svd(H)
  Hsv = svd.H$d
  # print(Hsv / Hsv[1])

  # construct a matrix to collect the selection criteria for the models with
  # state dimension s = 0:s.max
  stats = matrix(NA_real_, nrow = s.max+1, 5)
  colnames(stats) = c('s', 'n.par', 'Hsv', 'lndetSigma', 'criterion')
  stats[, 's'] = 0:s.max                # state dimension
  stats[, 'n.par'] = 2*stats[, 's']*m   # number of parameters
  # Hankel singular values
  # print(stats)
  # print(c(NA_real_, Hsv[1:s.max]))
  stats[, 'Hsv'] = c(NA_real_, Hsv[iseq(1,s.max)])

  info = list(m = m, Hsize = c(f,p), n.obs = n.obs, s.max = s.max, Hsv = Hsv)

  if (s.max == 0) {
    # maximum order = 0 ==> there is nothing to compute
    sigma = gamma0
    sigma_L = t(chol(sigma))
    model= list(s = 0, A = matrix(0, nrow = 0, ncol = 0), B = matrix(0, nrow = 0, ncol = m),
                C = matrix(0, nrow = m, ncol = 0), D = diag(m), sigma = sigma, sigma_L = sigma_L)
    stats[1, 'lndetSigma'] = 2*sum(log(diag(sigma_L)))
    if (keep_models) {
      models = list(model)
    } else {
      models = NULL
    }
    model = stspmod(sys = stsp(A = model$A, B = model$B, C = model$C, D = model$D),
                    sigma_L = model$sigma_L)
    return( list(model = model, models = models,
                 s = 0, stats = stats, info = info) )
  }

  # estimate order based on the Hankel singular values
  out.estorder = do.call(estorder, list(s.max = s.max, Hsv = Hsv,
                                        n.par = stats[, 'n.par'],
                                        m = m, n.obs = n.obs, Hsize=c(f,p), ...))

  s.hat = integer(0)
  if (!is.null(out.estorder)) {
    s.hat = out.estorder$s
    stats[, 'criterion'] = out.estorder$criterion
  }
  if (is.null(out.estorder) || (keep_models)) {
    s = (0:s.max)
  } else {
    s = s.hat
  }
  # cat('first round\n')
  # print(out.estorder)
  # cat(s, 's.hat', s.hat, '\n')

  models = vector('list', length(s))
  for (i in (1:length(s))) {
    si = s[i]
    j = which(stats[, 's'] == si)
    # estimate model with order s = s[i]
    out = try(aoki(m, si, svd.H, gamma0, Rp, Rf), silent = TRUE)
    if (!inherits(out, 'try-error')) {
      models[[i]] = out
      lndetSigma = 2*sum(log(diag(out$sigma_L)))
      stats[j, 'lndetSigma'] = lndetSigma
    } else {
      # 'aoki' failed!!!
      model[[i]] = list(s = si)
      stats[j, 'lndetSigma'] = NA_real_
    }
  }

  if (is.null(out.estorder)) {
    # estimate order based on the lndetSigma values
    out.estorder = do.call(estorder, list(s.max = s.max, Hsv = Hsv,
                                          lndetSigma = stats[, 'lndetSigma'],
                                          n.par = stats[, 'n.par'],
                                          m = m, n.obs = n.obs, Hsize=c(f,p), ...))
    if (is.null(out.estorder)) stop('could not estimate the model order?!')
    # cat('second round\n')
    # print(out.estorder)
    # cat(s, 's.hat', s.hat, '#models', length(models), '\n')
    s.hat = out.estorder$s
    stats[, 'criterion'] = out.estorder$criterion
  }


  if (length(s) == 1) {
    model = models[[1]]
  } else {
    model = models[[s.hat+1]]
  }
  if (is.null(model$A)) stop('AOKI failed for the selecetd model order!')
  model = stspmod( sys = stsp(A = model$A, B = model$B, C = model$C, D = model$D),
                   sigma_L = model$sigma_L )
  if (!keep_models) models = NULL

  return( list(model = model, models = models,
               s = s.hat, stats = stats, info = info ) )

}

est_stsp_cca = function(gamma, s.max, p, estorder = estorder_SVC,
                        keep_models = FALSE, n.obs = NULL, ...) {

  # check gamma
  if ( (!is.numeric(gamma)) || (!is.array(gamma)) ||
       (length(dim(gamma)) != 3) || (dim(gamma)[1] != dim(gamma)[2]) ) {
    stop('input "gamma" is not a valid 3-D array.')
  }
  d = dim(gamma)
  m = d[1]
  if (m == 0) {
    stop('the autocovariance function "gamma" is empty.')
  }
  lag.max = d[3] - 1
  if (lag.max <= 2) {
    stop('the autocovariance function must contain at least 2 lags.',
         ' lag.max=', lag.max)
  }
  if (is.null(n.obs)) {
    n.obs = Inf
  }
  n.obs = as.numeric(n.obs)[1]
  if (is.finite(n.obs)) {
    n.obs = as.integer(n.obs)[1]
    if (n.obs <= 0) stop('the sample size "n.obs" must be non negative')
  } else {
    n.obs = Inf
  }

  # p <=> past
  p = as.integer(p)[1]
  if (p < 1) stop('the number of lags "p" must be positive')
  if (lag.max < (2*p)) {
    stop('the autocovariance function must contain at least 2p lags.',
         ' lag.max=', lag.max, ' < 2*p=', 2*p)
  }
  # f <=> future
  f = p + 1

  s.max = as.integer(s.max)[1]
  if (s.max < 0) stop('maximum state dimension "s.max" must be non negative.')
  if (s.max > (p*m)) {
    stop('the number of lags "p" is too small for the desired maximum state dimension:',
         ' p*m=', p*m, ' < s.max=', s.max)
  }

  gamma0 = matrix(gamma[,,1], nrow = m, ncol = m)

  # covariance between "future" Y[t] = (y[t]', y[t+1]',...,y[t+f-1]')'
  #                  and "past" X[t] = (y[t-1]', y[t-2]',...,y[t-p]')'
  YX = bhankel(gamma[,,-1,drop=FALSE], d = c(f,p))

  # covariance matrix of "past" and "present" (y[t], y[t-1]', y[t-2]',...,y[t-p]')'
  X1X = btoeplitz(R = gamma[,,1:(p+1),drop = FALSE])
  # cat(dim(X1X), m,p, '\n')
  # covariance matrix of "past" X[t] = (y[t-1]', y[t-2]',...,y[t-p]')'
  XX = X1X[1:(m*p), 1:(m*p), drop = FALSE]
  # covariance between y[t] and "past" X[t] = (y[t-1]', y[t-2]',...,y[t-p]')'
  yX = X1X[1:m, (m+1):(m*(p+1)), drop = FALSE]
  # covariance between y[t] and the "next past" X[t+1] = (y[t]', y[t-1]',...,y[t-p+1]')'
  yX1 = XX[1:m, 1:(m*p), drop = FALSE]
  # covariance matrix of y[t] (=gamma[0])
  yy = XX[1:m, 1:m, drop = FALSE]
  # covariance of the "next past" X[t+1] = (y[t], y[t-1]', ...,y[t+1-p]')'
  #                and the "past" X[t] = (y[t-1]', y[t-2]',...,y[t-p]')'
  X1X = X1X[1:(m*p), (m+1):(m*(p+1)), drop = FALSE]

  # covariance matrix of the "future" Y[t] = (y[t]', y[t+1]',...,y[t+f-1]')'
  YY = btoeplitz(C = gamma[,,1:f,drop = FALSE])

  # cholesky factors of covariance matrices
  Rp = chol(XX)
  Rf = chol(YY)

  # weighted Hankel matrix: inv(t(Rf)) * YX * inv(Rp)
  # note: inv(t(Rf)) * YX * inv(Rp) is the covariance between the
  # "standardized" future and the "standardized" past.
  # The singular value decomposition of this matrix defines the
  # canonical correlation coefficients between future and past.
  YX = backsolve(Rf, YX, transpose = TRUE)
  YX = t(backsolve(Rp, t(YX), transpose = TRUE))

  # compute SVD of weighted Hankel matrix
  svd.YX = svd(YX, nu = 0)
  Hsv = svd.YX$d # Hankel singular values
  # print(Hsv / Hsv[1])

  # construct a matrix to collect the selection criteria for the models with
  # state dimension s = 0:s.max
  stats = matrix(NA_real_, nrow = s.max+1, 5)
  colnames(stats) = c('s', 'n.par', 'Hsv', 'lndetSigma', 'criterion')
  stats[, 's'] = 0:s.max                # state dimension
  stats[, 'n.par'] = 2*stats[, 's']*m   # number of parameters
  # Hankel singular values
  # print(stats)
  # print(c(NA_real_, Hsv[iseq(1,s.max)]))
  stats[, 'Hsv'] = c(NA_real_, Hsv[iseq(1,s.max)])

  info = list(m = m, Hsize = c(f,p), n.obs = n.obs, s.max = s.max, Hsv = Hsv)

  if (s.max == 0) {
    # maximum order = 0 ==> there is nothing to compute
    sigma = gamma0
    sigma_L = t(chol(sigma))
    model= list(s = 0, A = matrix(0, nrow = 0, ncol = 0), B = matrix(0, nrow = 0, ncol = m),
                C = matrix(0, nrow = m, ncol = 0), D = diag(m), sigma = sigma, sigma_L = sigma_L)
    stats[1, 'lndetSigma'] = 2*sum(log(diag(sigma_L)))
    if (keep_models) {
      models = list(model)
    } else {
      models = NULL
    }
    model = stspmod(sys = stsp(A = model$A, B = model$B, C = model$C, D = model$D),
                    sigma_L = model$sigma_L)
    return( list(model = model, models = models,
                 s = 0, stats = stats, info = info) )
  }

  # estimate order based on the Hankel singular values
  out.estorder = do.call(estorder, list(s.max = s.max, Hsv = Hsv,
                                        n.par = stats[, 'n.par'],
                                        m = m, n.obs = n.obs, Hsize=c(f,p), ...))

  s.hat = integer(0)
  if (!is.null(out.estorder)) {
    s.hat = out.estorder$s
    stats[, 'criterion'] = out.estorder$criterion
  }
  if (is.null(out.estorder) || (keep_models)) {
    s = (0:s.max)
  } else {
    s = s.hat
  }
  # cat('first round\n')
  # print(out.estorder)
  # cat(s, 's.hat', s.hat, '\n')

  # transform covariances for the model with the maximum order max(s)
  if (max(s) > 0) {
    # T = V' * inv(R_p')
    T = t(backsolve(Rp, svd.YX$v[, 1:max(s), drop = FALSE]))
    # check
    # print(T %*% XX %*% t(T)) # = identity
    yX = yX %*% t(T)
    yX1 = yX1 %*% t(T)
    X1X = T %*% X1X %*% t(T)
  }

  models = vector('list', length(s))
  D = diag(m)
  for (i in (1:length(s))) {
    # estimate model with order s = s[i]
    # state x[t] = X[1:si,t]
    si = s[i]
    j = which(stats[, 's'] == si)
    if (si == 0) {
      # there is nothing to do
      A = matrix(0, nrow = 0, ncol = 0)
      B = matrix(0, nrow = 0, ncol = m)
      C = matrix(0, nrow = m, ncol = 0)
      sigma = yy
    } else {
      # C is defined from the regression of y[t] on x[t]: C = Eyx (Exx^{-1})
      C = yX[, 1:si, drop = FALSE]

      # sigma = E u[t] u[t]' = gamma(0) - C E x[t] x[t]' C'
      sigma = yy - C %*% t(C)

      # A is defined from the regression of x[t+1] onto x[t]: A = Ex1x (Exx)^{-1}
      A = X1X[1:si, 1:si, drop = FALSE]

      # B is defined from the regression of x[t+1] onto u[t] = y[t] - Cx[t]
      # B = (E x[t+1] y[t]' - E x[t+1] x{t]' C' ) sigma^{-1}
      #   = (E x[t+1] y[t]' - A C' ) sigma^{-1}
      ux1 = yX1[ , 1:si, drop = FALSE] - C %*% t(A)
      B = t(solve(sigma, ux1))
    }
    sigma_L = t(chol(sigma))
    lndetSigma = 2*sum(log(diag(sigma_L)))

    models[[i]] = list(s = si, A=A, B=B, C=C, D = D, sigma = sigma, sigma_L = sigma_L)
    stats[j, 'lndetSigma'] = lndetSigma
  }

  if (is.null(out.estorder)) {
    # estimate order based on the lndetSigma values
    out.estorder = do.call(estorder, list(s.max = s.max, Hsv = Hsv,
                                          lndetSigma = stats[, 'lndetSigma'],
                                          n.par = stats[, 'n.par'],
                                          m = m, n.obs = n.obs, Hsize=c(f,p), ...))
    if (is.null(out.estorder)) stop('could not estimate the model order?!')
    # cat('second round\n')
    # print(out.estorder)
    # cat(s, 's.hat', s.hat, '#models', length(models), '\n')
    s.hat = out.estorder$s
    stats[, 'criterion'] = out.estorder$criterion
  }

  if (length(s) == 1) {
    model = models[[1]]
  } else {
    model = models[[s.hat+1]]
  }
  model = stspmod( sys = stsp(A = model$A, B = model$B, C = model$C, D = model$D),
                   sigma_L = model$sigma_L )
  if (!keep_models) models = NULL

  return( list(model = model, models = models,
               s = s.hat, stats = stats, info = info ) )
}


est_stsp_cca_sample = function(y, s.max, p, estorder = estorder_SVC, keep_models = FALSE,
                               mean_estimate = c('sample.mean', 'zero'), ...) {

  mean_estimate = match.arg(mean_estimate)

  y = try(as.matrix(y))
  if ( inherits(y, 'try-error') || (!is.numeric(y)) || (!is.matrix(y)) ) {
    stop('could not coerce the input "y" to a numeric matrix')
  }
  m = ncol(y)      # number of "outputs"
  n.obs = nrow(y)  # sample size
  if (m*n.obs == 0) stop('"y" contains no data')
  if (mean_estimate == 'sample.mean') {
    y.mean = colMeans(y)
    y = scale(y, center = y.mean, scale = FALSE)
  } else {
    y.mean = numeric(m)
  }

  # p <=> past
  p = as.integer(p)[1]
  if (p < 1) stop('the number of lags "p" must be positive')
  n.valid = n.obs - 2*p
  if ((n.valid-1) < m*(p+1)) {
    stop('the sample size is too small for the desired number of lags:',
         ' (n.obs-2*p-1)=', n.valid-1, ' < m*(p+1)=', m*(p+1))
  }
  # f <=> future
  f = p + 1

  s.max = as.integer(s.max)[1]
  if (s.max < 0) stop('maximum state dimension "s.max" must be non negative.')
  if (s.max > (p*m)) {
    stop('the number of lags "p" is too small for the desired maximum state dimension:',
         ' p*m=', p*m, ' < s.max=', s.max)
  }

  # "past" X[t] = (y[t-1]', y[t-2]',...,y[t-p]')'
  X = matrix(0, nrow = n.valid, ncol = m*p)
  for (i in (1:p)) X[, ((i-1)*m+1):(i*m)] = y[(p+1-i):(n.valid+p-i), ]
  # "future" y[t] = (y[t]', y[t+1]',...,y[t+f-1]')'
  Y = matrix(0, nrow = n.valid, ncol = m*f)
  for (i in (0:(f-1))) Y[, (i*m+1):((i+1)*m)] = y[(p+1+i):(n.valid+p+i), ]
  y = y[(p+1):(n.valid+p), ]

  # QR decomposition of Y,X
  X = qr.Q(qr(X)) # X is semiorthogonal: X'X = I
  Y = qr.Q(qr(Y)) # Y is semiorthogonal: Y'Y = I

  # compute SVD of weighted "Hankel matrix"
  svd.YX = svd((t(Y) %*% X), nu = 0)
  Hsv = svd.YX$d # Hankel singular values

  # construct a matrix to collect the selection criteria for the models with
  # state dimension s = 0:s.max
  stats = matrix(NA_real_, nrow = s.max+1, 5)
  colnames(stats) = c('s', 'n.par', 'Hsv', 'lndetSigma', 'criterion')
  stats[, 's'] = 0:s.max                # state dimension
  stats[, 'n.par'] = 2*stats[, 's']*m   # number of parameters
  # Hankel singular values
  # print(stats)
  # print(c(NA_real_, Hsv[iseq(1,s.max)]))
  stats[, 'Hsv'] = c(NA_real_, Hsv[iseq(1,s.max)])

  info = list(m = m, Hsize = c(f,p), n.obs = n.obs, s.max = s.max, Hsv = Hsv)

  if (s.max == 0) {
    # maximum order = 0 ==> there is nothing to compute
    sigma = ( t(y) %*% y ) / n.valid
    sigma_L = t(chol(sigma))
    model= list(s = 0, A = matrix(0, nrow = 0, ncol = 0), B = matrix(0, nrow = 0, ncol = m),
                C = matrix(0, nrow = m, ncol = 0), D = diag(m), sigma = sigma, sigma_L = sigma_L)
    stats[1, 'lndetSigma'] = 2*sum(log(diag(sigma_L)))
    if (keep_models) {
      models = list(model)
    } else {
      models = NULL
    }
    model = stspmod(sys = stsp(A = model$A, B = model$B, C = model$C, D = model$D),
                    sigma_L = model$sigma_L)
    return( list(model = model, models = models, y.mean = y.mean,
                 s = 0, stats = stats, info = info) )
  }

  # estimate order based on the Hankel singular values
  out.estorder = do.call(estorder, list(s.max = s.max, Hsv = Hsv,
                                        n.par = stats[, 'n.par'],
                                        m = m, n.obs = n.obs, Hsize=c(f,p), ...))

  s.hat = integer(0)
  if (!is.null(out.estorder)) {
    s.hat = out.estorder$s
    stats[, 'criterion'] = out.estorder$criterion
  }
  if (is.null(out.estorder) || (keep_models)) {
    s = (0:s.max)
  } else {
    s = s.hat
  }
  # cat('first round\n')
  # print(out.estorder)
  # cat(s, 's.hat', s.hat, '\n')

  # transform X -> X * V (for the model with the maximum order max(s))
  if (max(s) > 0) {
    X = X %*% svd.YX$v[, 1:max(s), drop = FALSE] #  X'X = identity!'
    # print(t(X) %*% X)
  }

  models = vector('list', length(s))
  D = diag(m)
  for (i in (1:length(s))) {
    # estimate model with order s = s[i]
    # state x[t] = X[1:si,t]
    si = s[i]
    j = which(stats[, 's'] == si)
    if (si == 0) {
      # there is nothing to do
      A = matrix(0, nrow = 0, ncol = 0)
      B = matrix(0, nrow = 0, ncol = m)
      C = matrix(0, nrow = m, ncol = 0)
      sigma = t(y) %*% y / n.valid
    } else {
      # C is defined from the regression of y[t] on x[t]: C = Eyx (Exx^{-1})
      # note X'X = I
      # cat(si, dim(X), s, '\n')
      C = t(y) %*% X[, 1:si, drop = FALSE]

      # sigma = E u[t] u[t]' = gamma(0) - C E x[t] x[t]' C'
      # u[t] = y[t] - C x[t]
      u = y - X[ , 1:si, drop = FALSE] %*% t(C)
      sigma = (t(u) %*% u ) / n.valid

      # AB is defined from the regression of x[t+1] onto (x[t]', u[t}')'
      AB = stats::lsfit(cbind(X[1:(n.valid-1), 1:si, drop = FALSE],
                              u[1:(n.valid-1), ,drop = FALSE]),
                        X[2:n.valid, 1:si, drop = FALSE],
                        intercept = FALSE) $coefficients
      AB = unname(AB)
      AB = matrix(AB, nrow = si+m, ncol = si)
      A = t(AB[1:si,,drop = FALSE])
      B = t(AB[(si+1):(si+m),,drop = FALSE])
    }
    sigma_L = t(chol(sigma))
    lndetSigma = 2*sum(log(diag(sigma_L)))

    models[[i]] = list(s = si, A=A, B=B, C=C, D = D, sigma = sigma, sigma_L = sigma_L)
    stats[j, 'lndetSigma'] = lndetSigma
  }

  if (is.null(out.estorder)) {
    # estimate order based on the lndetSigma values
    out.estorder = do.call(estorder, list(s.max = s.max, Hsv = Hsv,
                                          lndetSigma = stats[, 'lndetSigma'],
                                          n.par = stats[, 'n.par'],
                                          m = m, n.obs = n.obs, Hsize=c(f,p), ...))
    if (is.null(out.estorder)) stop('could not estimate the model order?!')
    # cat('second round\n')
    # print(out.estorder)
    # cat(s, 's.hat', s.hat, '#models', length(models), '\n')
    s.hat = out.estorder$s
    stats[, 'criterion'] = out.estorder$criterion
  }

  if (length(s) == 1) {
    model = models[[1]]
  } else {
    model = models[[s.hat+1]]
  }
  model = stspmod( sys = stsp(A = model$A, B = model$B, C = model$C, D = model$D),
                   sigma_L = model$sigma_L )
  if (!keep_models) models = NULL

  return( list(model = model, models = models, y.mean = y.mean,
               s = s.hat, stats = stats, info = info ) )
}

est_stsp_ss = function(obj, method = c('cca', 'aoki'),
                       s.max = NULL, p = NULL, p.ar.max = NULL, p.factor = 2,
                       extend_acf = FALSE, sample2acf = TRUE,
                       estorder = estorder_SVC, keep_models = FALSE,
                       mean_estimate = c('sample.mean', 'zero'), n.obs = NULL, ...) {
  method = match.arg(method)
  mean_estimate = match.arg(mean_estimate)

  if (!inherits(obj, 'autocov')) {
    y = try(as.matrix(obj))
    if ( inherits(y, 'try-error') || (!is.numeric(y)) || (!is.matrix(y)) ) {
      stop('could not coerce "obj" to a numeric matrix')
    }
    m = ncol(y)      # number of "outputs"
    n.obs = nrow(y)  # sample size (optional parameter n.obs is ignored)
    if (m == 0) stop('the "time series object" (obj) contains no data')
    if (n.obs < 3) stop('the sample size must be larger than 2')
    if (mean_estimate == 'sample.mean') {
      y.mean = colMeans(y)
      y = scale(y, center = y.mean, scale = FALSE)
    } else {
      y.mean = numeric(m)
    }
    lag.max = Inf
    gamma = NULL
  } else {
    gamma = obj$gamma
    d = dim(gamma)
    m = d[1]
    if (m == 0) {
      stop('the "autocov" object (obj) is empty.')
    }
    lag.max = d[3] - 1
    if (lag.max <= 2) {
      stop('the "autocov" object (obj) must contain at least 2 lags.',
           ' lag.max=', lag.max)
    }

    if (is.null(n.obs)) {
      n.obs = obj$n.obs
      if (is.null(n.obs)) n.obs = Inf
    }
    n.obs = as.numeric(n.obs)[1]
    if (is.finite(n.obs)) {
      n.obs = as.integer(n.obs)[1]
      if (n.obs < 3) stop('the sample size must be larger than 2')
    } else {
      n.obs = Inf
    }

    y.mean = rep(NA_real_, m)
    y = NULL
  }

  if ((is.null(gamma)) && (!sample2acf) && (method == 'aoki')) {
    stop('the "AOKI" method is only implemented for ACFs.')
  }

  # the "orders" p, p.factor, p.ar.max, ... must satisfy the following restrictions
  #
  # if extend_acf==TRUE (extrapolate ACF via AR(p) model): set ext = 2, else ext = 1
  #
  # 2 <= lag.max <= n.obs-1   (for a sample of size N = n.obs)  (=> n.obs >= 3)
  #
  # p.factor >= 1
  #
  # regression of y[t],...,y[t+p] on y[t-1],...,y[t-p]
  #   sample: n.obs - 2*p >= m*p                       => p <= n.obs/(m+2)
  #   ACF (Hankel matrix (p+1) x p): 2*p <= lag.max    => p <= lag.max/2 = ext*lag.max/2
  #      if extend_acf (extrapolate ACF via AR(p))     => p <= lag.max   = ext*lag.max/2
  #
  # p >= 1
  # p*m >= s.max                                       => p >= max(1, s.max/m)
  #
  # if p is undefined, then p is "estimated" via a long autoregression:
  #
  # AR model regression of y[t] on intercept and y[t-1],...,y[t-p]
  #   sample: n.obs - p >= p*m + 1                     => p.ar.max <= (n.obs-1)/(m+1)
  #   ACF:    p <= lag.max                             => p.ar.max <= lag.max
  #   default value for p.ar.max is                       p.ar.max = 10*log[10](n.obs)
  #
  #
  # p = p.factor*p.ar.hat <= p.factor*p.ar.max implies
  #                                         => p.factor*p.ar.max <= n.obs/(m+2)
  #                                         => p.factor*p.ar.max <= ext*lag.max/2
  #
  #                                         => p.factor*p.ar.max >= max(1, s.max/m)
  #
  # combining these restrictions (with p.factor >=1), we require:
  #                                         => p.ar.max <= (n.obs-1)/((m+2)*p.factor)
  #                                         => p.ar.max <= ext*lag.max/(2*p.factor)
  #
  #                                         => p.ar.max >= max(1, s.max/m)/p.factor

  ext = ifelse(extend_acf, 2, 1)
  if (is.null(s.max)) s.max = NA_integer_
  s.max = as.integer(s.max)[1]

  p.upper = floor(min(n.obs/(m+2), ext*lag.max/2))
  p.lower = ceiling(max(1, s.max/m, na.rm = TRUE))
  if (p.upper < p.lower) stop('the sample size or the number of lags is too small')

  if (is.null(p)) {
    # estimate the number of future/past lags (= size of the Hankel matrix)
    # by fitting an AR model
    p.ar.upper = floor(min((n.obs-1)/(m+2), ext*lag.max/2) / p.factor)
    p.ar.lower = ceiling(max(1, s.max/m, na.rm = TRUE) / p.factor)
    if (p.ar.upper < p.ar.lower) stop('the sample size or the number of lags is too small')
    if (is.null(p.ar.max)) {
      p.ar.max = max(p.ar.lower, min(p.ar.upper, round(10*log10(n.obs))))
    }
    # check p.ar.max
    if ((p.ar.lower > p.ar.max) || (p.ar.upper < p.ar.max)) {
      stop('the maximum AR order "p.ar.max" is too small or too large: ',
           'lower bound=',p.ar.lower, ', upper bound=', p.ar.upper, 'but p.ar.max=', p.ar.max)
    }

    ar_model = NULL
    p.ar.hat = -1
    # estimate the number of future/past lags (= size of the Hankel matrix)
    # by fitting an AR model
    if ( (is.null(gamma)) && (!sample2acf) ) {
      p.ar.hat = est_ar_ols(y, p.max = p.ar.max, penalty = 2/n.obs)$p
    } else {
      if (is.null(gamma)) {
        gamma = autocov(y, lag.max = 2*p.upper/ext, demean = FALSE)$gamma
      }
      ar_model = est_ar_yw(gamma, p.max = p.ar.max, penalty = 2/n.obs)
      p.ar.hat = ar_model$p
    }
    p = max(p.lower, as.integer(p.factor * p.ar.hat))
  }

  # check p
  p = as.integer(p)[1]
  if ( (p < p.lower) || (p > p.upper)) {
    stop('the number of (past) lags "p" is to small or too large for the given data.',
         ' p=', p)
  }


  # compute ACF if needed
  if ((is.null(gamma)) && (sample2acf)) {
    gamma = autocov(y, lag.max = 2*p/ext, demean = FALSE)$gamma
  }

  # exzend ACF if needed
  if ((!is.null(gamma)) && (extend_acf)) {
    if (is.null(ar_model) || (p.ar.hat != p)) {
      ar_model = est_ar_yw(gamma, p.max = p, penalty = -1)
    }
    gamma = rbind(bmatrix(gamma[,,1:(p+1),drop = FALSE], rows = c(1,3)),
                  matrix(0, nrow = m*p, ncol = m))
    a = bmatrix(ar_model$a[,,p:1,drop = FALSE], rows = 1)
    for (lag in ((p+1):(2*p))) {
      # print(c(lag,(m*lag+1),(m*(lag+1)),(m*(lag-p)+1),(m*lag)))
      gamma[(m*lag+1):(m*(lag+1)), ] = a %*% gamma[(m*(lag-p)+1):(m*lag), ,drop = FALSE]
    }
    dim(gamma) = c(m, 2*p+1, m)
    gamma = aperm(gamma,c(1,3,2))
  }

  # check s.max
  if (is.na(s.max)) {
    s.max = p*m
  }
  s.max = as.integer(s.max)[1]
  if (s.max < 0) stop('maximum state dimension "s.max" must be non negative.')
  if (s.max > (p*m)) {
    # this should not happen, due to the restrictions on p
    stop('the number of lags "p" is too small for the desired maximum state dimension:',
         ' p*m=', p*m, ' < s.max=', s.max)
  }

  n.valid = n.obs - 2*p

  if (is.null(gamma)) {
    out = est_stsp_cca_sample(y, s.max = s.max, p = p, estorder = estorder,
                              keep_models = keep_models,
                              mean_estimate = 'zero', ...)
  } else {
    if (method == 'aoki') {
      out = est_stsp_aoki(gamma, s.max = s.max, p = p, estorder = estorder,
                          keep_models = keep_models, n.obs = n.obs, ...)
    } else {
      out = est_stsp_cca(gamma, s.max = s.max, p = p, estorder = estorder,
                         keep_models = keep_models, n.obs = n.obs, ...)
    }
  }

  out$y.mean = y.mean
  # out$gamma = gamma
  return(out)
}


# templates (WS) #####################################
#

tmpl_sigma_L = function(sigma_L, structure = c('as_given', 'chol', 'symm', 'identity', 'full_normalized')) {
  if ( (!is.numeric(sigma_L)) || (!is.matrix(sigma_L)) || (ncol(sigma_L) != nrow(sigma_L)) ) {
    stop('"sigma_L" is not a square, numeric matrix')
  }
  structure = match.arg(structure)

  n = ncol(sigma_L)
  if (n == 0) {
    h = numeric(0)
    H = matrix(0, nrow =0, ncol = 0)
    return(list(h = h, H = H, n.par = 0))
  }
  u = upper.tri(sigma_L, diag = FALSE)

  if (structure == 'identity') {
    sigma_L = diag(n)
  }
  if (structure == 'full_normalized') {
    sigma_L = diag(n)
    sigma_L[sigma_L == 0] = NA
  }
  if (structure == 'chol') {
    sigma_L[u] = 0
  }
  if (structure == 'symm') {
    sigma_L = ( sigma_L + t(sigma_L) ) /2
    sigma_L[u] = 0
  }

  # print(sigma_L)
  h = as.vector(sigma_L)
  ix = which(is.na(h))
  h[ix] = 0
  n.par = length(ix)
  # cat(n.par, ix)
  H = matrix(0, nrow = n^2, ncol = n.par)
  # print(H)
  if (n.par > 0) H[ix, ] = diag(n.par)

  if ((structure == 'symm') && (n > 1)) {
    i = matrix(1:(n*n), ncol = n, nrow = n)
    iU = i[u]
    iL = t(i)[u]
    h[iU] = h[iL]
    H[iU, ] = H[iL, ]
  }

  return(list(h = h, H = H, n.par = n.par))
}

model2template = function(model, sigma_L = c("as_given", "chol", "symm", "identity", "full_normalized")) {

  # Check inputs and obtain integer-valued parameters
  if ( !( inherits(model, 'armamod') || inherits(model, 'stspmod') || inherits(model, 'rmfdmod') ) ) {
    stop('model is not an "armamod", "rmfdmod", or "stspmod" object.')
  }
  order = unname(dim(model$sys))
  n = order[2]

  # Affine parametrisation
  h = as.vector(unclass(model$sys))
  ix = which(!is.finite(h))
  h[ix] = 0
  n.par = length(ix)
  H = matrix(0, nrow = length(h), ncol = n.par)
  H[ix,] = diag(n.par)

  # Set noise parameters
  sigma_L_structure = match.arg(sigma_L)
  sigma_L = model$sigma_L

  junk = tmpl_sigma_L(sigma_L, structure = sigma_L_structure)
  h = c(h, junk$h)
  H = bdiag(H, junk$H)

  return(list(h = h, H = H, class = class(model)[1],
              order = order, n.par = ncol(H)))
}

tmpl_rmfd_echelon = function(nu, m = length(nu), sigma_L = c("chol", "symm", "identity", "full_normalized")) {

  # (m,n) transfer function k = d(z) c^(-1)(z) with degrees deg(c) = p, deg(d) = q
  sigma_L = match.arg(sigma_L)
  nu = as.integer(nu)
  n = length(nu)
  if ( (n < 1) || (m < 1) ) stop('illegal dimension, (m,n) must be positive')
  if (min(nu) < 0) stop('Kronecker indices must be non negative')

  p = max(nu)
  order = c(m, n, p, p)

  # code the position of the basis columns of the Hankel matrix
  basis = rationalmatrices::nu2basis(nu)

  # coefficients of c(z) in reverse order!!!
  # c = [c[p]',...,c[0]']'
  c = matrix(0, nrow = n*(p+1), ncol = n)
  # d = [d[0]',...,d[p]']'
  d = matrix(0, nrow = m*(p+1), ncol = n)

  # code free entries with NA's
  for (i in (1:n)) {
    shift = (p-nu[i])*n
    k = nu[i]*n + i
    basis_i = basis[basis < k]
    # cat(i, shift, k, basis_i, '\n')
    c[k + shift, i] = 1
    c[basis_i + shift, i] = NA_real_
    d[iseq(m + 1, (nu[i] + 1)*m), i] = NA_real_   # i-th column has degree nu[i]
    d[iseq(n + 1, m), i] = NA_real_               # the last (m-n) rows of b[0] are free
  }
  # print(c)
  # reshuffle c -> c = [c[0]',...,c[p]']'
  dim(c) = c(n, p+1, n)
  c = c[, (p+1):1, , drop = FALSE]
  dim(c) = c(n*(p+1), n)

  # print(rbind(c, d))

  sL = matrix(NA_real_, nrow = n, ncol = n)

  # create a helper model
  sys = structure(rbind(c, d), order = order, class = c('rmfd','ratm'))
  model = list(sys = sys, sigma_L = sL, names = NULL, label = NULL)
  model = structure(model, class = c('rmfdmod', 'rldm'))
  # print(model$sys)

  # create template
  tmpl = model2template(model, sigma_L = sigma_L)
  # add Kronecker indices
  tmpl$nu = nu

  # the first min(n,m) rows of d[0] and c[0] are equal!
  # matrix of linear indices
  i = matrix(1:((n+m)*(p+1)*n), nrow = (n+m)*(p+1), ncol = n)
  ic = as.vector(i[1:min(n,m), 1:n])
  id = as.vector(i[(n*(p+1)+1):(n*(p+1)+min(n,m)), 1:n])
  # print(rbind(ia,ib))
  tmpl$h[id] = tmpl$h[ic]
  tmpl$H[id, ] = tmpl$H[ic, ]

  return(tmpl)
}

tmpl_arma_pq = function(m, n, p, q, sigma_L = c("chol", "symm", "identity", "full_normalized")) {

  m = as.integer(m)[1]
  n = as.integer(n)[1]
  p = as.integer(p)[1]
  q = as.integer(q)[1]
  if (min(m-1, n-1, p, q) < 0) stop('illegal dimensions/orders: m,n >= 1 and p,q >= 0 must hold')

  sigma_L = match.arg(sigma_L)

  order = c(m, n, p, q)

  sys = matrix(NA_real_, nrow = m, ncol = m*(p+1) + n*(q+1))
  # set a[0] to the identity matrix
  sys[ , 1:m] = diag(m)
  # the first min(m,n) rows of b[0] and a[0] coincide
  sys[1:min(m,n), (m*(p+1) + 1):(m*(p+1) + n)] = diag(x = 1, nrow = min(m,n), ncol = n)

  sL = matrix(NA_real_, nrow = n, ncol = n)

  # create a helper model
  sys = structure(sys, order = order, class = c('lmfd','ratm'))
  model = list(sys = sys, sigma_L = sL, names = NULL, label = NULL)
  model = structure(model, class = c('armamod', 'rldm'))

  # create template
  tmpl = model2template(model, sigma_L = sigma_L)

  return(tmpl)
}

tmpl_stsp_full = function(m, n, s, sigma_L = c("chol", "symm", "identity", "full_normalized")) {

  m = as.integer(m)[1]
  n = as.integer(n)[1]
  s = as.integer(s)[1]
  if (min(m-1, n-1, s) < 0) stop('illegal dimensions/orders: m,n >= 1 and s >= 0 must hold')

  sigma_L = match.arg(sigma_L)

  order = c(m, n, s)

  sys = matrix(NA_real_, nrow = (m+s), ncol = (n+s))
  sys[(s+1):(s+min(m,n)),(s+1):(s+n)] = diag(x = 1, nrow = min(m,n), ncol = n)

  sL = matrix(NA_real_, nrow = n, ncol = n)

  # create a helper model
  sys = structure(sys, order = order,  class = c('stsp','ratm'))
  model = list(sys = sys, sigma_L = sL, names = NULL, label = NULL)
  model = structure(model, class = c('stspmod', 'rldm'))

  # create template
  tmpl = model2template(model, sigma_L = sigma_L)

  return(tmpl)
}

tmpl_stsp_ar = function(m, p, sigma_L = c("chol", "symm", "identity", "full_normalized")) {

  m = as.integer(m)[1]
  p = as.integer(p)[1]
  if (min(m-1, p) < 0) stop('illegal dimensions/orders: m >= 1 and p >= 0 must hold')

  sigma_L = match.arg(sigma_L)

  n = m
  s = m*p
  order = c(m, n, s)

  if (s == 0) {
    sys = diag(m)
  } else {
    C = matrix(NA_real_, nrow = m, ncol = s)
    A = rbind(matrix(0, nrow = m, ncol = s),
              diag(x = 1, nrow = m*(p-1), ncol = s))
    B = diag(x = 1, nrow = s, ncol = m)
    D = diag(m)
    sys = rbind( cbind(A,B), cbind(C, D))
  }

  sL = matrix(NA_real_, nrow = n, ncol = n)

  # create a helper model
  sys = structure(sys, order = order,  class = c('stsp', 'ratm'))
  model = list(sys = sys, sigma_L = sL, names = NULL, label = NULL)
  model = structure(model, class = c('stspmod', 'rldm'))

  # create template
  tmpl = model2template(model, sigma_L = sigma_L)

  # take care of the restriction A[1:m,] = C
  if (p > 0) {
    i = matrix(1:((m+s)*(m+s)), nrow = m+s, ncol = m+s)
    iA = i[1:m, 1:s]
    iC = i[(s+1):(s+m), 1:s]
    tmpl$H[iA, ] = tmpl$H[iC, ]
  }

  return(tmpl)
}

tmpl_stsp_echelon = function(nu, n = length(nu), sigma_L = c("chol", "symm", "identity", "full_normalized")) {

  # (m,n) transfer function C(I z^-1 -A)^-1 + D
  sigma_L = match.arg(sigma_L)
  nu = as.integer(nu)
  m = length(nu)
  if ( (n < 1) || (m < 1) ) stop('illegal dimension, (m,n) must be positive')
  if (min(nu) < 0) stop('Kronecker indices must be non negative')
  s = sum(nu) # state space dimension

  D = diag(x = 1, nrow = m, ncol = n)
  if (m > n) D[(n+1):m, ] = NA_real_

  if (s == 0) {
    sys = structure(D, order = c(m, n, s),  class = c('stsp','ratm'))
  } else {
    basis = nu2basis(nu)
    AC = matrix(0, nrow = s + m, ncol = s)
    dependent = c(basis + m, 1:m)
    for (i in (1:length(dependent))) {
      d = abs(basis-dependent[i])
      if (min(d) == 0) {
        # dependent[i]-th row is in basis
        j = which(d == 0)
        AC[i, j] = 1
      } else {
        j = which(basis < dependent[i])
        AC[i, j] = NA_real_
      }
    }
    B = matrix(NA_real_, nrow = s, ncol = n)
    sys = structure(cbind( AC, rbind(B,D)), order = c(m, n, s),  class = c('stsp','ratm'))
  }

  sL = matrix(NA_real_, nrow = n, ncol = n)

  # create a helper model
  model = list(sys = sys, sigma_L = sL, names = NULL, label = NULL)
  model = structure(model, class = c('stspmod', 'rldm'))
  # print(model$sys)

  # create template
  tmpl = model2template(model, sigma_L = sigma_L)
  # add Kronecker indices
  tmpl$nu = nu

  return(tmpl)
}


fill_template = function(th, template) {
  if (template$class == 'armamod') {
    order = template$order
    m = order[1]
    n = order[2]
    p = order[3]
    q = order[4]
    P = template$h + template$H %*% th
    sys = matrix(P[1:(length(P) - n*n)], nrow = m, ncol = m*(p+1)+n*(q+1))
    sys = structure(sys, order = order, class = c('lmfd','ratm'))

    sigma_L = matrix(P[((length(P) - n*n + 1)):length(P)],
                     nrow = n, ncol = n)
    model = list(sys = sys, sigma_L = sigma_L, names = NULL, label = NULL)
    model = structure(model, class = c('armamod', 'rldm'))
    return(model)
  }
  if (template$class == 'rmfdmod') {
    order = template$order
    m = order[1]
    n = order[2]
    p = order[3]
    q = order[4]
    P = template$h + template$H %*% th
    sys = matrix(P[1:(length(P) - n*n)], nrow = n*(p+1) + m*(q+1), ncol = n)
    sys = structure(sys, order = order, class = c('rmfd','ratm'))

    sigma_L = matrix(P[((length(P) - n*n + 1)):length(P)],
                     nrow = n, ncol = n)
    model = list(sys = sys, sigma_L = sigma_L, names = NULL, label = NULL)
    model = structure(model, class = c('rmfdmod', 'rldm'))
    return(model)
  }
  if (template$class == 'stspmod') {
    order = template$order
    m = order[1]
    n = order[2]
    s = order[3]
    P = template$h + template$H %*% th
    sys = matrix(P[1:(length(P) - n*n)], nrow = (m+s), ncol = (n+s))
    sys = structure(sys, order = order, class = c('stsp','ratm'))

    sigma_L = matrix(P[((length(P) - n*n + 1)):length(P)],
                     nrow = n, ncol = n)
    model = list(sys = sys, sigma_L = sigma_L, names = NULL, label = NULL)
    model = structure(model, class = c('stspmod', 'rldm'))
    return(model)
  }

  stop('illegal template')
}


extract_theta = function(model, template,
                         tol = sqrt(.Machine$double.eps), on_error = c('ignore','warn','stop'), ignore_sigma_L = FALSE) {

  on_error = match.arg(on_error)

  if ( !( (class(model)[1] == template$class) && all(dim(model$sys) == template$order) ) ) {
    stop('model and template are not compatible.')
  }

  P = c(as.vector(unclass(model$sys)), as.vector(model$sigma_L))
  if (ncol(template$H) == 0) return(double(0))

  out = stats::lsfit(template$H, P - template$h, intercept = FALSE)

  if (on_error != 'ignore') {
    res = out$residuals
    if (ignore_sigma_L) {
      n = template$order[2]
      res = res[iseq(1,length(res) - n^2)]   # ignore sigma_L
    }
    # is length(res)=0 possible?
    if (length(res) > 0) {
      if (max(abs(res)) > tol) {
        if (on_error == 'warn') warning(paste0('model does not match template. max(abs(res)) = ', max(abs(res))))
        if (on_error == 'stop') stop(paste0('model does not match template. max(abs(res)) = ', max(abs(res))))
      }
    }
  }

  th = unname(out$coef)
  return(th)
}

r_model = function(template, ntrials.max = 100, bpoles = NULL, bzeroes = NULL,
                       rand.gen = stats::rnorm, ...) {

  # Initialize variable
  n.par = template$n.par
  constraint_satisfied = FALSE
  ntrials = 0

  # take care of the non-square case
  if (template$order[1] != template$order[2]) bzeroes = NULL

  # Generate random models until the constraints are satisfied
  while ( (!constraint_satisfied) && (ntrials < ntrials.max) ) {
    ntrials = ntrials + 1
    theta = rand.gen(n.par, ...)
    model = fill_template(theta, template)

    # Check constraints on poles and zeros
    constraint_satisfied = TRUE
    if (!is.null(bpoles)) {
      poles = poles(model, print_message = FALSE)
      if (length(poles) > 0) {
        constraint_satisfied = (min(abs(poles)) > bpoles)
      }
    }
    if ( constraint_satisfied && (!is.null(bzeroes)) ) {
      zeroes = zeroes(model, print_message = FALSE)
      if (length(zeroes) > 0) {
        constraint_satisfied = (min(abs(zeroes)) > bzeroes)
      }
    }
  }

  # Throw error if constraint is not satisfied after a certain number of trials
  if (!constraint_satisfied){
    stop('Could not generate a suitable model with ', ntrials, ' trials!')
  }

  return(model)
}


# tools.R ###################################################
#

test_armamod = function(dim = c(1,1), degrees = c(1,1), b0 = NULL, sigma_L = NULL,
                        digits = NULL, bpoles = NULL, bzeroes = NULL, n.trials = 100) {
  # check input parameter "dim"
  dim = as.integer(dim) # note: as.integer converts to vector!
  if ((length(dim) != 2) || (dim[1] <= 0) || (dim[2] < 0)) {
    stop('argument "dim" must be a vector of integers with length 2, dim[1] > 0 and dim[2] >= 0!')
  }
  # check input parameter "degrees"
  degrees = as.integer(degrees) # note: as.integer converts to vector!
  if ((length(degrees) != 2) || (degrees[1] < 0) || (degrees[2] < (-1))) {
    stop('argument "degrees" must be a vector of integers with length 2, degrees[1] >= 0 and degrees[2] >= -1!')
  }

  m = dim[1]
  n = dim[2]
  p = degrees[1]
  q = degrees[2]

  if (p == 0) bpoles = NULL
  if ( (m != n) || (q <= 0) ) bzeroes = NULL

  # check input parameter "b0"
  if (!is.null(b0)) {
    if ( (!is.numeric(b0)) || (!is.matrix(b0)) || any(dim(b0) != dim) ) {
      stop('"b0" must be a compatible, numeric matrix')
      i = is.na(b0)
      theta = stats::rnorm(sum(i))
      if (!is.null(digits)) theta = round(theta, digits)
      b0[i] = theta
    }
  }
  else {
    b0 = diag(x = 1, nrow = m, ncol = n)
  }

  # check input parameter "sigma_L"
  if (!is.null(sigma_L)) {
    if ( (!is.numeric(sigma_L)) || (!is.matrix(sigma_L)) || any(dim(sigma_L) != n) ) {
      stop('"sigma_L" must be a compatible, numeric matrix')
    }
  }
  else {
    sigma_L = matrix(stats::rnorm(n*n), nrow = n, ncol = n)
    if (n >0) {
      sigma_L = t(chol(sigma_L %*% t(sigma_L)))
    }
    if (!is.null(digits)) sigma_L = round(sigma_L, digits)
  }

  i.trial = 0
  err = TRUE
  sd = 1
  while ( (i.trial < n.trials) && (err) ) {
    a = cbind(diag(m), matrix(stats::rnorm(m*m*p, sd = sd), nrow = m, ncol = m*p))
    dim(a) = c(m,m,p+1)
    if (q >= 0) {
      b = cbind(b0, matrix(stats::rnorm(m*n*q, sd = sd), nrow = m, ncol = n*q))
      dim(b) = c(m, n, q+1)
    } else {
      b = array(0, dim = c(m, n, q+1))
    }
    if (!is.null(digits)) {
      a = round(a, digits)
      b = round(b, digits)
    }
    sys = lmfd(a, b)

    err = FALSE
    if ( !is.null(bpoles) ) {
      err = try(min(abs(poles(sys, print_message = FALSE))) <= bpoles, silent = TRUE)
      if (inherits(err, 'try-error')) err = TRUE
    }
    if ( (!err) && (!is.null(bzeroes)) ) {
      err = try((min(abs(zeroes(sys, print_message = FALSE))) <= bzeroes), silent = TRUE)
      if (inherits(err, 'try-error')) err = TRUE
    }
    i.trial = i.trial + 1
    sd = sd/1.1
  }
  if (err) {
    stop('Could not generate a suitable ARMA model with ', n.trials, ' trials!')
  }

  model = armamod(sys = sys, sigma_L = sigma_L)
  return(model)
}

test_stspmod = function(dim = c(1,1), s = NULL, nu = NULL, D = NULL, sigma_L = NULL,
                     digits = NULL, bpoles = NULL, bzeroes = NULL, n.trials = 100) {


  # check input parameter "dim"
  dim = as.integer(dim) # note: as.integer converts to vector!
  if ((length(dim) != 2) || (min(dim) < 0)) {
    stop('argument "dim" must be a vector of non negative integers with length 2!')
  }
  m = dim[1]
  n = dim[2]

  # check input parameter "D"
  if (!is.null(D)) {
    if ( (!is.numeric(D)) || (!is.matrix(D)) || any(dim(D) != dim) ) {
      stop('"D" must be a compatible, numeric matrix')
    }
  }
  else {
    D = diag(x = 1, nrow = m, ncol = n)
  }

  # check input parameter "sigma_L"
  if (!is.null(sigma_L)) {
    if ( (!is.numeric(sigma_L)) || (!is.matrix(sigma_L)) || any(dim(sigma_L) != n) ) {
      stop('"sigma_L" must be a compatible, numeric matrix')
    }
  }
  else {
    sigma_L = matrix(stats::rnorm(n*n), nrow = n, ncol = n)
    if (n >0) {
      sigma_L = t(chol(sigma_L %*% t(sigma_L)))
    }
    if (!is.null(digits)) sigma_L = round(sigma_L, digits)
  }

  # check input parameter "nu"
  if (!is.null(nu)) {
    nu = as.integer(nu) # as.integer converts to vector
    if ((length(nu) != m) || (min(nu) < 0)) {
      stop('"nu" must be a vector of non negative integers with length equal to dim[1]!')
    }
  } else {
    if (is.null(s)) stop('either "s" or "nu" must be specified')
    s = as.integer(s)[1]
    if (s < 0) stop('parameter "s" must be a nonnegative integer')
  }

  sys = try(test_stsp(dim = dim, s = s, nu = nu, D = D, digits = digits,
                      bpoles = bpoles, bzeroes = bzeroes, n.trials = n.trials))
  if (inherits(sys, 'try-error')) {
    stop('Could not generate a suitable statespace model with ', n.trials, ' trials!')
  }

  model = stspmod(sys = sys, sigma_L = sigma_L)
  return(model)
}
# tools_riccati.R

riccati = function(A, M, C, G, only.X = TRUE) {
  m = nrow(M) # number of states
  n = ncol(M) # number of outputs
  if (is.null(C)) C = matrix(0, nrow = n, ncol = m)

  # setup generalized eigenvalue problem (of dimension (2m-by-2m))
  # AA*X = BB*X*Lambda
  AA = diag(1,2*m)
  BB = AA

  tA = A - M %*% solve(G, C)
  AA[1:m,1:m] = t(tA)
  AA[(m+1):(2*m),1:m] = - M %*% solve(G, t(M))

  BB[(m+1):(2*m),(m+1):(2*m)] = tA
  BB[1:m,(m+1):(2*m)] = -t(C) %*% solve(G, C)

  # compute QZ decomposition of (AA,BB):
  # AA = Q*S*Z^H, BB = Q*T*Z^H where S,Q are upper triangular and Q,Z are unitary
  out = QZ::qz.dgges(AA, BB, vsl = FALSE, vsr = TRUE, LWORK = NULL)

  # throw a warning if qz.dgges 'fails'
  if (out$INFO<0) {
    warning('qz.dgges was not successful!')
  }

  # eigenvalues
  lambda = out$ALPHA / out$BETA

  # select eigenvalues with modulus < 1
  select = abs(lambda)<1

  # the number of stable eigenvalues should be equal to m (number of states)
  if (sum(select)!=m) {
    warning('number of eigenvalues with modulus less than one is not equal to m!')
  }

  # reorder the QZ decomposition such that the "stable" eigenvalues are on the top!
  out = QZ::qz.dtgsen(S=out$S, T=out$T, Q = out$Z, Z=out$Z, select = select, ijob = 0L,
                      want.Q = FALSE, want.Z = TRUE, LWORK = NULL, LIWORK = NULL)

  # throw a warning if qz.dtgsen fails
  if (out$INFO<0) {
    warning('qz.dtgsen was not successful!')
  }
  # recompute (ordered) eigenvalues
  lambda = out$ALPHA / out$BETA

  # compute X from the eigenvectors Z[,1:m]
  X = solve(t(out$Z[1:m,1:m,drop=FALSE]), t(out$Z[(m+1):(2*m),1:m,drop=FALSE])) # note X should be symmetric!
  X = Re( X + t(X) )/2

  # if (only.X) {
  #   return(list(X = X, info = out$INFO, lambda = lambda, AA = AA, BB = BB, out))
  # } else {
  #   sigma = G - C %*% X %*% t(C)
  #   B = t(solve(t(sigma),t(M - A %*% X %*% t(C))))
  #   return(list(X = X, B = B, sigma = sigma, info = out$INFO, lambda = lambda, AA = AA, BB = BB, out))
  # }
  if (only.X) {
    return(X)
  } else {
    sigma = G - C %*% X %*% t(C)
    B = t(solve(t(sigma),t(M - A %*% X %*% t(C))))
    return(list(X = X, B = B, sigma = sigma, lambda = lambda))
  }

}
