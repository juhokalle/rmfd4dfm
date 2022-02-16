upgrade_objects = function(force = TRUE, ...) {
  args = list(...)
  n_args = length(args)
  if (n_args == 0) {
    # return empty list
    return(as.vector(integer(0), mode = 'list'))
  }

  # skip NULL's
  not_null = sapply(args, FUN = function(x) {!is.null(x)})
  args = args[not_null]
  n_args = length(args)
  if (n_args == 0) {
    # return empty list
    return(as.vector(integer(0), mode = 'list'))
  }
  if ((n_args == 1) && (!force)) {
    # just return arg(s)
    return(args)
  }

  classes = sapply(args, FUN = function(x) {class(x)[1]} )
  grades = match(classes, c('polm', 'lpolm', 'lmfd','rmfd','stsp','pseries','zvalues'), nomatch = 0)

  # coerce arguments to a common class = "highest" class
  k = which.max(grades)
  max_grade = grades[k]

  # coerce matrix to polm object
  if (max_grade == 0) {
    # this should not happen ?!
    obj = try(polm(args[[k]]))
    if (inherits(obj, 'try-error')) stop('could not coerce object to "polm" object')
    args[[k]] = obj
    max_grade = 1
    classes[k] = 'polm'
    grades[k] = 1
  }

  # If the argument with the highest grade is an mfd,
  # then coerce it to a stsp object and update max_grade
  if (max_grade %in% c(3,4)) {
    obj = as.stsp(args[[k]])
    args[[k]] = obj
    max_grade = 5
    classes[k] = 'stsp'
    grades[k] = 5
  }

  for (i in (1:n_args)) {
    if ((grades[i] < max_grade) && ( i != k )) {
      if (grades[i] == 0) {
        obj = try(polm(args[[i]]))
        if (inherits(obj, 'try-error')) stop('could not coerce object to "polm" object')
        args[[i]] = obj
      }

      # max_grade of all objects involved corresponds to lpolm()
      if (max_grade == 2) {
        obj = try(as.lpolm(args[[i]]))
        if (inherits(obj, 'try-error')) stop('could not coerce object to "stsp" object')
        args[[i]] = obj
      }

      # max_grade of all objects involved corresponds to state space model
      if (max_grade == 5) {
        if (grades[i] == 2) { stop('cannot coerce lpolm to state space model') }

        obj = try(as.stsp(args[[i]]))
        if (inherits(obj, 'try-error')) stop('could not coerce object to "stsp" object')
        args[[i]] = obj
      }

      # max_grade of all objects involved corresponds to power series model
      if (max_grade == 6) {
        if (grades[i] == 2) { stop('cannot coerce lpolm to pseries') }
        obj = pseries(args[[i]], lag.max = dim(args[[k]])[3])
        args[[i]] = obj
      }

      # max_grade of all objects involved corresponds to zvalues object
      if (max_grade == 7) {
        if (grades[i] == 2) { stop('cannot coerce lpolm to frequency response') }
        if (grades[i] == 6) { stop('cannot coerce power series to frequency response') }
        z = attr(args[[k]], 'z')
        obj = zvalues(args[[i]], z = z)
        args[[i]] = obj
      }
    }
  }

  return(args)
}

inflate_object = function(e, m, n) {
  cl = class(e)[1]

  if (cl == 'polm') {
    e = unclass(prune(e, tol = 0))
    e = as.vector(e)
    if ((length(e)*m*n) == 0) {
      return(polm(matrix(0, nrow = m, ncol = n)))
    }
    a = array(e, dim = c(length(e), m, n))
    a = aperm(a, c(2,3,1))
    return(polm(a))
  }

  if (cl == 'lpolm') {

    # Some hacking such that prune() can be applied
    attr_e = attributes(e)
    attributes(e) = NULL
    e = array(e, dim = attr_e$dim) %>% polm()

    # prune returns a polm object, so it is necessary to transform it back to lpolm()
    e = unclass(prune(e, tol = 0))
    e = as.vector(e)
    if ((length(e)*m*n) == 0) {
      return(lpolm(matrix(0, nrow = m, ncol = n), min_deg = 0))
    }
    a = array(e, dim = c(length(e), m, n))
    a = aperm(a, c(2,3,1))
    return(lpolm(a, min_deg = attr_e$min_deg))

  }

  if (cl == 'stsp') {
    if ((m*n) == 0) {
      return(stsp(A = matrix(0, nrow = 0, ncol = 0),
                  B = matrix(0, nrow = 0, ncol = n),
                  C = matrix(0, nrow = m, ncol = 0),
                  D = matrix(0, nrow = m, ncol = n)))
    }
    A = e$A
    s = ncol(A)
    B = e$B
    C = e$C
    D = e$D[1,1]
    if (n < m) {
      A = do.call(bdiag, args = rep(list(A), n))
      B = do.call(bdiag, args = rep(list(B), n))
      C = as.vector(C)
      C = do.call(cbind, args = rep(list(matrix(C, nrow = m, ncol = length(C), byrow = TRUE)), n))
    } else {
      A = do.call(bdiag, args = rep(list(A), m))
      C = do.call(bdiag, args = rep(list(C), m))
      B = as.vector(B)
      B = do.call(rbind, args = rep(list(matrix(B, nrow = length(B), ncol = n)), m))
    }
    D = matrix(D, nrow = m, ncol = n)
    return(stsp(A, B, C, D))
  }

  if (cl == 'pseries') {
    max.lag = unname(dim(e)[3])
    if (((max.lag + 1)*m*n) == 0) {
      a = array(0, dim = c(m,n,max.lag+1))
      a = structure(a, class = c('pseries', 'ratm'))
      return(a)
    }
    e = unclass(e)
    e = as.vector(e)
    a = array(e, dim = c(max.lag+1, m, n))
    a = aperm(a, c(2,3,1))
    a = structure(a, class = c('pseries', 'ratm'))
    return(a)
  }

  if (cl == 'zvalues') {
    z = attr(e, 'z')
    if ((length(z)*m*n) == 0) {
      a = array(0, dim = c(m,n,length(z)))
      a = structure(a, z = z, class = c('zvalues', 'ratm'))
      return(a)
    }
    e = unclass(e)
    e = as.vector(e)
    a = array(e, dim = c(length(z), m, n))
    a = aperm(a, c(2,3,1))
    a = structure(a, z = z, class = c('zvalues', 'ratm'))
    return(a)
  }

  stop('unsupported class "', cl, '" for inflation of scalars')
}

# Internal functions acting on ARRAYS for matrix multiplication ####

convolve_3D = function(a, b, truncate = FALSE) {

  # a,b must be two compatible arrays
  # a is an (m,n) dimensional polynomial of degree p
  d = dim(a)
  m = d[1]
  n = d[2]
  p = d[3] - 1

  # d is an (n,o) dimensional polynomial of degree q
  d = dim(b)
  if (d[1] != n) stop('arguments are not compatible')
  o = d[2]
  q = d[3] - 1

  # output c = a*b is an (m,o) dimensional polynomial of degree r

  if (truncate) {
    # multiplication of two "pseries" objects,
    # in this case we only need the powers of c(z) up to degree r = min(p,q)
    r = min(p, q)

    # if any of the arguments is an empty pseries, or a pseries with lag.max -1
    # ensure that lag.max = r >= 0
    if (min(c(m, n, o, r+1)) == 0) return(array(0, dim = c(m, o, max(r, 0)+1)))

    if (p > r) a = a[,,1:(r+1), drop = FALSE] # chop a to degree min(p,q)
    if (q > r) b = b[,,1:(r+1), drop = FALSE] # chop b to degree min(p,q)
    p = r
    q = r
  } else {
    # multiplication of two polynomials

    # if any of the arguments is an empty polynomial, or a polynomial of degree -1
    if (min(c(m, n, o, p+1, q+1)) == 0) return(array(0, dim = c(m, o, 0)))

    # degree = sum of the two degrees
    r = p + q
  }

  # cat('truncate', truncate,'\n')
  if (p <= q) {
    c = matrix(0, nrow = m, ncol = o*(r+1))
    b = matrix(b, nrow = n, ncol = o*(q+1))
    for (i in (0:p)) {
      # compute a[i] * (b[0], ..., b[q]) and add to (c[i], ..., c[i+q])
      # however, if i+q > r (truncate = TRUE) then
      # compute a[i] * (b[0], ..., b[r-i]) and add to (c[i], ..., c[r])
      j1 = i
      j2 = min(j1+q+1, r+1)
      # cat('a*b', i, j1, j2, '|', dim(c), ':', (j1*o + 1), (j2*o),
      #      '|', dim(b), ':', 1, ((j2-j1)*o), '\n')
      c[ , (j1*o + 1):(j2*o)] =
        c[ , (j1*o + 1):(j2*o)] +
        matrix(a[,,i+1], nrow = m, ncol = n) %*% b[,1:((j2-j1)*o) , drop = FALSE]
    }
    c = array(c, dim = c(m,o,r+1))
    return(c)
  } else {
    # if the degree of b is smaller than the degree of a (q < p)
    # then we first compute b' * a'
    c = matrix(0, nrow = o, ncol = m*(r+1))
    a = matrix(aperm(a, c(2,1,3)), nrow = n, ncol = m*(p+1))
    b = aperm(b, c(2,1,3))
    for (i in (0:q)) {
      j1 = i
      j2 = min(j1+p+1, r+1)
      # cat('b*a', i, j1, j2, '|', dim(c), ':', (j1*m + 1), (j2*m),
      # '|', dim(a), ':', 1, ((j2-j1)*m), '\n')
      c[ , (j1*m + 1):(j2*m)] =
        c[ , (j1*m + 1):(j2*m)] +
        matrix(b[,,i+1], nrow = o, ncol = n) %*% a[,1:((j2-j1)*m) , drop = FALSE]
    }
    c = array(c, dim = c(o,m,r+1))
    c = aperm(c, c(2,1,3))
  }
  return(c)
}

# # multiplication of matrix polynomials
# # this function performs only basic checks on the inputs!
# mmult_poly = function(a, b) {
#   # a,b must be two compatible arrays
#   da = dim(a)
#   db = dim(b)
#   if (da[2] != db[1]) stop('arguments are not compatible')
#   # if any of the arguments is an empty polynomial, or a polynomial of degree -1
#   if (min(c(da,db)) == 0) return(array(0, dim = c(da[1], db[2], 0)))
#
#   # skip zero leading coefficients
#   if (da[3] > 0) {
#     a = a[ , , rev(cumprod(rev(apply(a == 0, MARGIN = 3, FUN = all)))) == 0, drop = FALSE]
#     da = dim(a)
#   }
#   # skip zero leading coefficients
#   if (db[3] > 0) {
#     b = b[ , , rev(cumprod(rev(apply(b == 0, MARGIN = 3, FUN = all)))) == 0, drop = FALSE]
#     db = dim(b)
#   }
#   # if any of the arguments is an empty polynomial, or a polynomial of degree -1
#   if (min(c(da,db)) == 0) return(array(0, dim = c(da[1], db[2], 0)))
#
#   pa = da[3] - 1
#   pb = db[3] - 1
#
#   x = array(0, dim = c(da[1], db[2], pa + pb + 1))
#   # the 'convolution' of the coefficients is computed via a double loop
#   # of course this could be implemented more efficiently!
#   for (i in (0:(pa+pb))) {
#     for (k in iseq(max(0, i - pb), min(pa, i))) {
#       x[,,i+1] = x[,,i+1] +
#         matrix(a[,,k+1], nrow = da[1], ncol = da[2]) %*% matrix(b[,,i-k+1], nrow = db[1], ncol = db[2])
#     }
#   }
#   return(x)
# }
#
#
# # internal function
# # multiplication of two impulse response functions
# # this function performs only basic checks on the inputs!
# # almost equal to mmult_poly
# mmult_pseries = function(a, b) {
#   # a,b must be two compatible arrays
#   da = dim(a)
#   db = dim(b)
#   if (da[2] != db[1]) stop('arguments are not compatible')
#   # if any of the arguments is an empty pseries, or a pseries with lag.max = -1
#   if (min(c(da,db)) == 0) return(array(0, dim = c(da[1], db[2], min(da[3], db[3]))))
#
#   lag.max = min(da[3], db[3]) - 1
#   # truncate to the minimum lag.max
#   # a = a[ , , 1:(lag.max+1), drop = FALSE]
#   # b = b[ , , 1:(lag.max+1), drop = FALSE]
#
#   x = array(0, dim = c(da[1], db[2], lag.max + 1))
#   # the 'convolution' of the impulse response coefficients is computed via a double loop
#   # of course this could be implemented more efficiently!
#   for (i in (0:lag.max)) {
#     for (k in (0:i)) {
#       x[,,i+1] = x[,,i+1] +
#         matrix(a[,,k+1], nrow = da[1], ncol = da[2]) %*% matrix(b[,,i-k+1], nrow = db[1], ncol = db[2])
#     }
#   }
#   return(x)
# }


# Internal functions acting on ARRAYS for elementwise matrix multiplication ####

# internal function
# univariate polynomial division c = a / b
poly_div = function(a, b) {
  # a,b are vectors

  # take care of the case that the leading coefficients of b are zero ( e.g. b = c(1,2,0,0))
  if (length(b) > 0) {
    b = b[rev(cumprod(rev(b == 0))) == 0]
  }
  lb = length(b)
  if (lb == 0) {
    stop('illegal polynomial division (b is zero)')
  }

  # take care of the case that the leading coefficients of a are zero ( e.g. a = c(1,2,0,0))
  if (length(a) > 0) {
    a = a[rev(cumprod(rev(a == 0))) == 0]
  }
  la = length(a)

  if (la < lb) return(0)   # deg(a) < deg(b)

  if (lb == 1) return(a/b) # deg(b) = 0

  a = rev(a)
  b = rev(b)
  c = rep.int(0, la - lb + 1)
  i = la - lb + 1
  while (i > 0) {
    d = a[1]/b[1]
    c[i] = d
    a[1:lb] = a[1:lb] - d*b
    a = a[-1]
    i = i - 1
  }
  return(c)
}

# internal function
# univariate polynomial remainder r: a = b * c + r
poly_rem = function(a, b) {
  # a,b are vectors

  # take care of the case that the leading coefficients of b are zero ( e.g. b = c(1,2,0,0))
  if (length(b) > 0) {
    b = b[rev(cumprod(rev(b == 0))) == 0]
    lb = length(b)
  } else {
    lb = 0
  }
  if (lb == 0) {
    stop('illegal polynomial division (b is zero)')
  }

  # take care of the case that the leading coefficients of a are zero ( e.g. a = c(1,2,0,0))
  if (length(a) > 0) {
    a = a[rev(cumprod(rev(a == 0))) == 0]
    la = length(a)
  } else {
    la = 0
  }

  if (la < lb) return(a)   # deg(a) < deg(b)

  if (lb == 1) {
    return(0)
  }

  a = rev(a)
  b = rev(b)
  while (length(a) >= lb) {
    d = a[1]/b[1]
    a[1:lb] = a[1:lb] - d*b
    a = a[-1]
  }
  return( rev(a) )
}

# internal function
# elementwise multiplication of matrix polynomials
# this function performs only basic checks on the inputs!
emult_poly = function(a, b) {
  # a,b must be two compatible arrays
  da = dim(a)
  db = dim(b)
  if ( (da[1] != db[1]) || (da[2] != db[2]) ) stop('arguments are not compatible')
  # if any of the arguments is an empty polynomial, or a polynomial of degree -1
  if (min(c(da,db)) == 0) return(array(0, dim = c(da[1], da[2], 0)))

  # skip zero leading coefficients
  if (da[3] > 0) {
    a = a[ , , rev(cumprod(rev(apply(a == 0, MARGIN = 3, FUN = all)))) == 0, drop = FALSE]
    da = dim(a)
  }
  # skip zero leading coefficients
  if (db[3] > 0) {
    b = b[ , , rev(cumprod(rev(apply(b == 0, MARGIN = 3, FUN = all)))) == 0, drop = FALSE]
    db = dim(b)
  }
  # if any of the arguments is an empty polynomial, or a polynomial of degree -1
  if (min(c(da,db)) == 0) return(array(0, dim = c(da[1], da[2], 0)))

  pa = da[3] - 1
  pb = db[3] - 1

  x = array(0, dim = c(da[1], db[2], pa + pb + 1))
  # the 'convolution' of the coefficients is computed via a double loop
  # of course this could be implemented more efficiently!
  for (i in (0:(pa+pb))) {
    for (k in iseq(max(0, i - pb), min(pa, i))) {
      x[,,i+1] = x[,,i+1] +
        matrix(a[,,k+1], nrow = da[1], ncol = da[2]) * matrix(b[,,i-k+1], nrow = db[1], ncol = db[2])
    }
  }
  return(x)
}


# internal function
# elementwise multiplication of two impulse response functions
# this function performs only basic checks on the inputs!
# almost equal to emult_poly
emult_pseries = function(a, b) {
  # a,b must be two compatible arrays
  da = dim(a)
  db = dim(b)
  if ( (da[1] != db[1]) || (da[2] != db[2]) ) stop('arguments are not compatible')
  # if any of the arguments is an empty pseries, or a pseries with lag.max = -1
  if (min(c(da,db)) == 0) return(array(0, dim = c(da[1], da[2], min(da[3], db[3]))))

  lag.max = min(da[3], db[3]) - 1
  # truncate to the minimum lag.max
  # a = a[ , , 1:(lag.max+1), drop = FALSE]
  # b = b[ , , 1:(lag.max+1), drop = FALSE]

  x = array(0, dim = c(da[1], db[2], lag.max + 1))
  # the 'convolution' of the impulse response coefficients is computed via a double loop
  # of course this could be implemented more efficiently!
  for (i in (0:lag.max)) {
    for (k in (0:i)) {
      x[,,i+1] = x[,,i+1] +
        matrix(a[,,k+1], nrow = da[1], ncol = da[2]) * matrix(b[,,i-k+1], nrow = db[1], ncol = db[2])
    }
  }
  return(x)
}


# internal function
# elementwise multiplication of two rational vectors (in stsp form)
emult_stsp_vek = function(a,b) {

  dim_a = dim(a)
  m = dim_a[1]
  s_a = dim_a[3]
  dim_b = dim(b)
  s_b = dim_b[3]

  # convert a to a diagonal matrix
  A = a$A
  B = a$B
  C = a$C
  D = a$D
  A = do.call(bdiag, args = rep(list(A), m))
  B = do.call(bdiag, args = rep(list(B), m))
  C = do.call(bdiag, args = lapply(as.vector(1:m, mode = 'list'),
                                   FUN = function(x) C[x,,drop = FALSE]))
  D = diag(x = as.vector(D), nrow = m, ncol = m)
  da = stsp(A, B, C, D)
  # print(da)

  # multiply with b
  ab = da %r% b
  # print(ab)
  # print(pseries(ab) - pseries(a)*pseries(b))

  # controllability matrix
  Cm = ctr_matrix(ab)
  svd_Cm = svd(Cm, nv = 0)
  # print(svd_Cm$d)

  # Cm has rank <= s = s_a+s_b (should we check this?)
  # state transformation
  A = t(svd_Cm$u) %*% ab$A %*% svd_Cm$u
  B = t(svd_Cm$u) %*% ab$B
  C = ab$C %*% svd_Cm$u
  D = ab$D
  s = s_a + s_b

  # print(ctr_matrix(stsp(A,B,C,D)))

  # skip the "non-controllable" states
  ab = stsp(A[1:s,1:s, drop = FALSE], B[1:s,,drop = FALSE],
            C[,1:s,drop = FALSE], D)

  # print(ctr_matrix(ab))
  # print(svd(ctr_matrix(ab))$d)

  return(ab)
}


# internal function
# elementwise multiplication of two rational vectors (in stsp form)
emult_stsp = function(a,b) {
  dim_a = dim(a)
  m = dim_a[1]
  n = dim_a[2]
  if ((m*n) == 0) {
    return(stsp(A = matrix(0, nrow = 0, ncol = 0), B = matrix(0, nrow = 0, ncol = n),
                C = matrix(0, nrow = m, ncol = 0), D = matrix(0, nrow = m, ncol = n)))
  }

  s_a = dim_a[3]
  s_b = dim(b)[3]

  if ((s_a+s_b)==0) {
    return(stsp(a$A, a$B, a$C, a$D * b$D))
  }

  if (n <= m) {
    cols = vector(n, mode = 'list')
    # compute elementwise multiplication of the columns of a and b
    for (i in (1:n)) {
      cols[[i]] = emult_stsp_vek(a[,i], b[,i])
    }
    # print(cols)
    # bind the columns
    ab = do.call(cbind, cols)
    return(ab)
  }

  # consider the transposed marices
  ab = t(emult_stsp(t(a), t(b)))
  return(ab)
}


# Bind methods ####

# bind methods
#

rbind.ratm = function(...) {

  args = upgrade_objects(force = FALSE, ...)
  n_args = length(args)
  if (n_args == 0) return(NULL)
  if (n_args == 1) return(args[[1]])

  # check number of columns
  n_cols = vapply(args, FUN = function(x) {dim(x)[2]}, 0)
  if (min(n_cols) != max(n_cols)) stop('the number of columns of the matrices does not coincide')

  # combine the first two arguments ##############
  cl1 = class(args[[1]])[1]

  # polynomial matrices ##########################
  if (cl1 == 'polm') {
    x1 = unclass(args[[1]])
    d1 = dim(x1)
    x2 = unclass(args[[2]])
    d2 = dim(x2)
    d = c(d1[1]+d2[1], d1[2], max(d1[3],d2[3]))
    if (min(d) == 0) {
      # empty polynomial
      args[[2]] = polm(array(0, dim = d))
    } else {
      if (d1[3] < d2[3]) x1 = dbind(d = 3, x1, array(0, dim = c(d1[1], d1[2], d2[3]-d1[3])))
      if (d2[3] < d1[3]) x2 = dbind(d = 3, x2, array(0, dim = c(d2[1], d2[2], d1[3]-d2[3])))
      args[[2]] = polm(dbind(d = 1, x1, x2))
    }
    if (n_args == 2) return(args[[2]])
    return(do.call(rbind, args[-1]))
  }

  # Laurent polynomial matrices ##########################
  if (cl1 == 'lpolm') {

      attr_e1 = attributes(args[[1]])
      attr_e2 = attributes(args[[2]])

      dim1 = attr_e1$dim
      dim2 = attr_e2$dim

      q1 = attr_e1$min_deg
      p1 = dim1[3]-1+q1

      q2 = attr_e2$min_deg
      p2 = dim2[3]-1+q2

      min_q = min(q1, q2)
      max_p = max(p1, p2)

      e1 = unclass(args[[1]])
      e2 = unclass(args[[2]])

      e1 = dbind(d = 3,
                 array(0, dim = c(dim1[1], dim1[2], -min_q+q1)),
                 e1,
                 array(0, dim = c(dim1[1], dim1[2], max_p-p1)))
      e2 = dbind(d = 3,
                 array(0, dim = c(dim2[1], dim2[2], -min_q+q2)),
                 e2,
                 array(0, dim = c(dim2[1], dim1[2], max_p-p2)))
      dim1 = dim(e1)
      dim2 = dim(e2)

      d = c(dim1[1]+dim2[1], dim1[2], max(dim1[3],dim2[3]))
    if (min(d) == 0) {
      # empty polynomial
      args[[2]] = lpolm(array(0, dim = d), min_deg = min_q)
    } else {
      args[[2]] = lpolm(dbind(d = 1, e1, e2), min_deg = min_q)
    }
    if (n_args == 2) return(args[[2]])
    return(do.call(rbind, args[-1]))
  }

  # statespace realizations ##########################
  if (cl1 == 'stsp') {
    x1 = args[[1]]
    x2 = args[[2]]
    A = bdiag(x1$A, x2$A)
    B = rbind(x1$B, x2$B)
    C = bdiag(x1$C, x2$C)
    D = rbind(x1$D, x2$D)
    args[[2]] = stsp(A,B,C,D)
    if (n_args == 2) return(args[[2]])
    return(do.call(rbind, args[-1]))
  }

  # pseries functions ##############################
  if (cl1 == 'pseries') {
    x1 = unclass(args[[1]])
    d1 = dim(x1)
    x2 = unclass(args[[2]])
    d2 = dim(x2)
    d = c(d1[1]+d2[1], d1[2], min(d1[3],d2[3]))
    if (min(d) == 0) {
      # empty pseries
      x = array(0, dim = d)
    } else {
      x = dbind(d = 1, x1[,,1:d[3],drop = FALSE], x2[,,1:d[3],drop = FALSE])
    }
    x = structure(x, class = c('pseries','ratm'))
    args[[2]] = x
    if (n_args == 2) return(args[[2]])
    return(do.call(rbind, args[-1]))
  }

  # zvalues functions ##############################
  if (cl1 == 'zvalues') {
    # print(args[[1]])
    # print(args[[2]])
    z1 = attr(args[[1]],'z')
    z2 = attr(args[[2]],'z')
    if (!isTRUE(all.equal(z1,z2))) {
      stop('the complex points z of the two frequency responses do not coincide')
    }
    x1 = unclass(args[[1]])
    x2 = unclass(args[[2]])
    x = dbind(d = 1, x1, x2)
    x = structure(x, z = z1, class = c('zvalues','ratm'))
    args[[2]] = x
    if (n_args == 2) return(args[[2]])
    return(do.call(rbind, args[-1]))
  }

  stop('(rbind) this should not happen')
}


cbind.ratm = function(...) {
  # args = list(...)
  args = upgrade_objects(force = FALSE, ...)
  n_args = length(args)
  if (n_args == 0) return(NULL)
  if (n_args == 1) return(args[[1]])

  # check number of rows
  n_rows = sapply(args, FUN = function(x) {dim(x)[1]})
  if (min(n_rows) != max(n_rows)) stop('the number of rows of the matrices does not coincide')

  # combine the first two arguments ##############
  cl1 = class(args[[1]])[1]

  # polynomial matrices ##########################
  if (cl1 == 'polm') {
    x1 = unclass(args[[1]])
    d1 = dim(x1)
    x2 = unclass(args[[2]])
    d2 = dim(x2)
    d = c(d1[1], d1[2] + d2[2], max(d1[3],d2[3]))
    if (min(d) == 0) {
      # empty polynomial
      args[[2]] = polm(array(0, dim = d))
    } else {
      if (d1[3] < d2[3]) x1 = dbind(d = 3, x1, array(0, dim = c(d1[1], d1[2], d2[3]-d1[3])))
      if (d2[3] < d1[3]) x2 = dbind(d = 3, x2, array(0, dim = c(d2[1], d2[2], d1[3]-d2[3])))
      args[[2]] = polm(dbind(d = 2, x1, x2))
    }
    if (n_args == 2) return(args[[2]])
    return(do.call(cbind, args[-1]))
  }

  if (cl1 == 'lpolm') {

    attr_e1 = attributes(args[[1]])
    attr_e2 = attributes(args[[2]])

    dim1 = attr_e1$dim
    dim2 = attr_e2$dim

    q1 = attr_e1$min_deg
    p1 = dim1[3]-1+q1

    q2 = attr_e2$min_deg
    p2 = dim2[3]-1+q2

    min_q = min(q1, q2)
    max_p = max(p1, p2)

    e1 = unclass(args[[1]])
    e2 = unclass(args[[2]])

    min_deg_e = min(attr_e1$min_deg, attr_e2$min_deg)

    e1 = dbind(d = 3,
               array(0, dim = c(dim1[1], dim1[2], -min_q+q1)),
               e1,
               array(0, dim = c(dim1[1], dim1[2], max_p-p1)))
    e2 = dbind(d = 3,
               array(0, dim = c(dim2[1], dim2[2], -min_q+q2)),
               e2,
               array(0, dim = c(dim2[1], dim1[2], max_p-p2)))

    d = c(dim1[1], dim1[2]+dim2[2], max(dim1[3],dim2[3]))
    if (min(d) == 0) {
      # empty polynomial
      args[[2]] = lpolm(array(0, dim = d), min_deg = min_deg_e)
    } else {
      args[[2]] = lpolm(dbind(d = 2, e1, e2), min_deg = min_deg_e)
    }
    if (n_args == 2) return(args[[2]])
    return(do.call(cbind, args[-1]))
  }

  # statespace realizations ##########################
  if (cl1 == 'stsp') {
    x1 = args[[1]]
    x2 = args[[2]]
    A = bdiag(x1$A, x2$A)
    B = bdiag(x1$B, x2$B)
    C = cbind(x1$C, x2$C)
    D = cbind(x1$D, x2$D)
    args[[2]] = stsp(A,B,C,D)
    if (n_args == 2) return(args[[2]])
    return(do.call(cbind, args[-1]))
  }

  # pseries functions ##############################
  if (cl1 == 'pseries') {
    x1 = unclass(args[[1]])
    d1 = dim(x1)
    x2 = unclass(args[[2]])
    d2 = dim(x2)
    d = c(d1[1], d1[2] + d2[2], min(d1[3],d2[3]))
    if (min(d) == 0) {
      # empty pseries
      x = array(0, dim = d)
    } else {
      x = dbind(d = 2, x1[,,1:d[3],drop = FALSE], x2[,,1:d[3],drop = FALSE])
    }
    x = structure(x, class = c('pseries','ratm'))
    args[[2]] = x
    if (n_args == 2) return(args[[2]])
    return(do.call(cbind, args[-1]))
  }

  # zvalues functions ##############################
  if (cl1 == 'zvalues') {
    # print(args[[1]])
    # print(args[[2]])
    z1 = attr(args[[1]],'z')
    z2 = attr(args[[2]],'z')
    if (!isTRUE(all.equal(z1,z2))) {
      stop('the complex points z of the two frequency responses do not coincide')
    }
    x1 = unclass(args[[1]])
    x2 = unclass(args[[2]])
    x = dbind(d = 2, x1, x2)
    x = structure(x, z = z1, class = c('zvalues','ratm'))
    args[[2]] = x
    if (n_args == 2) return(args[[2]])
    return(do.call(cbind, args[-1]))
  }

  stop('(cbind) this should not happen')
}



# Important functions: %r% and group arithmetic ####

'%r%' = function(e1, e2) {

  if ( !(inherits(e1, 'ratm') || inherits(e2, 'ratm') )) {
    stop('one of the arguments must be a rational matrix object (ratm)')
  }

  # print(class(e1))
  # print(class(e2))

  out = upgrade_objects(force = TRUE, e1, e2)
  e1 = out[[1]]
  e2 = out[[2]]

  # print(class(e1))
  # print(class(e2))

  d1 = unname(dim(e1))
  cl1 = class(e1)[1]
  # Not needed according to RStudio: gr1 = match(cl1, c('polm','lmfd','rmfd','stsp','pseries','zvalues'), nomatch = 0)
  d2 = unname(dim(e2))
  # Not needed according to RStudio: cl2 = class(e2)[1]
  # Not needed according to RStudio: gr2 = match(cl1, c('polm','lmfd','rmfd','stsp','pseries','zvalues'), nomatch = 0)

  if (d1[2] != d2[1]) stop('the rational matrices e1, e2 are not compatible (ncol(e1) != nrow(e2))')

  # finally do the computations
  if (cl1 == 'polm') {
    # e = polm(mmult_poly(unclass(e1), unclass(e2)))
    e = polm(convolve_3D(unclass(e1), unclass(e2)))
    e = prune(e, tol = 0)
    return(e)
  }

  if (cl1 == 'lpolm') {
    attr_e1 = attributes(e1)
    attr_e2 = attributes(e2)
    # e = lpolm(mmult_poly(unclass(e1), unclass(e2)),
    #           min_deg = attr_e1$min_deg + attr_e2$min_deg)
    e = lpolm(convolve_3D(unclass(e1), unclass(e2)),
              min_deg = attr_e1$min_deg + attr_e2$min_deg)

    attr_e = attributes(e)
    e = prune(polm(unclass(e)), tol = 0)

    e = lpolm(array(c(e), dim = attr_e$dim), min_deg = attr_e$min_deg)
    return(e)
  }


  if (cl1 == 'stsp') {
    A1 = e1$A
    B1 = e1$B
    C1 = e1$C
    D1 = e1$D
    A2 = e2$A
    B2 = e2$B
    C2 = e2$C
    D2 = e2$D
    A = rbind(cbind(A1,                                         B1 %*% C2),
              cbind(matrix(0,nrow = nrow(A2), ncol = ncol(A1)), A2))
    B = rbind(B1 %*% D2, B2)
    C = cbind(C1, D1 %*% C2)
    D = D1 %*% D2
    # print(A)
    e = stsp(A, B, C, D)
    return(e)
  }

  if (cl1 == 'pseries') {
    # e = mmult_pseries(unclass(e1), unclass(e2))
    e = convolve_3D(unclass(e1), unclass(e2), TRUE)
    class(e) = c('pseries','ratm')
    return(e)
  }

  if (cl1 == 'zvalues') {
    z1 = attr(e1,'z')
    z2 = attr(e2,'z')
    if (!isTRUE(all.equal(z1, z2))) {
      stop('the complex points z of the two frequency responses do not coincide')
    }
    e1 = unclass(e1)
    e2 = unclass(e2)
    e = array(0, dim = c(d1[1], d2[2], length(z1)))
    for (i in (1:length(z1))) {
      e[,,i] = matrix(e1[,,i], nrow = d1[1], ncol = d1[2]) %*% matrix(e2[,,i], nrow = d2[1], ncol = d2[2])
    }

    e = structure(e, z = z1, class = c('zvalues', 'ratm'))
    return(e)
  }

  stop('this should not happen')
}





Ops.ratm = function(e1, e2) {

# unary operator +/- ###############################################################
  if (missing(e2)) {

    if (.Generic == '+') {
      return(e1)
    }

    if (.Generic == '-') {
      cl1 = class(e1)[1]
      if (cl1 == 'polm') {
        return(polm(-unclass(e1)))
      }
      if (cl1 == 'lpolm') {
        attr_e1 = attributes(e1)
        return(lpolm(-unclass(e1), min_deg = attr_e1$min_deg))
      }
      if (cl1 == 'lmfd') {
        b = polm(-unclass(e1$b))
        return(lmfd(a = e1$a, b = b))
      }
      if (cl1 == 'rmfd') {
        d = polm(-unclass(e1$d))
        return(rmfd(c = e1$c, d = d))
      }
      if (cl1 == 'stsp') {
        return(stsp(A = e1$A, B = e1$B, C = -e1$C, D = -e1$D))
      }
      if (cl1 == 'pseries') {
        e1 = -unclass(e1)
        e1 = structure(e1, class = c('pseries', 'ratm'))
        return(e1)
      }
      if (cl1 == 'zvalues') {
        z = attr(e1,'z')
        e1 = -unclass(e1)
        e1 = structure(e1, z = z, class = c('zvalues', 'ratm'))
        return(e1)
      }
      stop('unsupported class: "',cl1,'"')
    }
    stop('unsupported unary operator: "',.Generic,'"')
  }

# power operator e1^n ###########################################
  if (.Generic == '^') {
    d1 = unname(dim(e1))
    cl1 = class(e1)[1]
    a2 = unclass(e2)

    if ( ( length(a2) != 1 ) || (a2 %% 1 != 0 ) ) {
      stop('unsupported power!')
    }

    if ( d1[1] != d1[2] ) {
      stop('power operation is only defined for non empty, square rational matrices!')
    }

    # __a2 == 1 ######################
    if (a2 == 1) {
      return(e1)
    } # a2 = 1

    # __a2 == 0 ######################
    if (a2 == 0) {
      if (cl1 == 'polm') {
        return(polm(diag(d1[1])))
      }
      if (cl1 == 'lpolm') {
        return(lpolm(diag(d1[1]), min_deg = 0))
      }
      if (cl1 == 'lmfd') {
        return(lmfd(a = diag(d1[1]), b = diag(d1[1])))
      }
      if (cl1 == 'rmfd') {
        return(rmfd(c = diag(d1[1]), d = diag(d1[1])))
      }
      if (cl1 == 'stsp') {
        return(stsp(A = matrix(0, nrow = 0, ncol = 0),
                    B = matrix(0, nrow = 0, ncol = d1[1]),
                    C = matrix(0, nrow = d1[1], ncol = 0),
                    D = diag(d1[1])))
      }
      if (cl1 == 'pseries') {
        if (d1[1] > 0) {
          e1 = unclass(e1)
          e1[,,-1] = 0
          e1[,,1] = diag(d1[1])
          e1 = structure(e1, class = c('pseries', 'ratm'))
        }
        return(e1)
      }
      if (cl1 == 'zvalues') {
        if (d1[1] > 0) {
          z = attr(e1,'z')
          e1 = array(diag(d1[1]), dim = d1)
          e1 = structure(e1, z = z, class = c('zvalues', 'ratm'))
        }
        return(e1)
      }
      stop('unsupported class: "',cl1,'"')
    } # a2 = 0

    # upgrade "lmfd" to "stsp" objects
    if (cl1 == 'lmfd') {
      e1 = as.stsp(e1)
      cl1 = 'stsp'
      d1 = unname(dim(e1))
    }

    # upgrade "rmfd" to "stsp" objects
    if (cl1 == 'rmfd') {
      e1 = as.stsp(e1)
      cl1 = 'stsp'
      d1 = unname(dim(e1))
    }

    # __a2 > 1 ######################
    if (a2 > 1) {
      e = e1
      for (i in (2:a2)) {
        e = e %r% e1
      }
      return(e)
    } # a2 > 1

    # __a2 < 0 ######################

    if (d1[1] <= 0) {
      stop('power operation with negative power is only defined for non empty, square rational matrices!')
    }

    # convert "polm" to "stsp" objects
    if (cl1 == 'polm') {
      e1 = as.stsp(e1)
      cl1 = 'stsp'
      d1 = unname(dim(e1))
    }

    if (cl1 == 'lpolm') {
      stop("Negative powers of Laurent polynom object cannot be taken.")
    }

    if (cl1 == 'stsp')  {
      # compute inverse
      A = e1$A
      B = e1$B
      C = e1$C
      D = e1$D
      D = try(solve(D), silent = TRUE)
      if (inherits(D, 'try-error')) {
        stop('could not compute state space representation of inverse (D is singular)')
      }
      B = B %*% D
      e1 = stsp(A = A - B %*% C, B = B, C = -D %*% C, D = D)

      e = e1
      for (i in iseq(2,abs(a2))) {
        e = e %r% e1
      }
      return(e)
    }

    if (cl1 == 'pseries')  {
      # compute inverse
      a = unclass(e1)
      m = dim(a)[1]           # we could also use d1!
      lag.max = dim(a)[3] - 1
      if (lag.max < 0) {
        # this should not happen?!
        stop('impulse response contains no lags!')
      }

      # b => inverse impulse response
      b = array(0, dim = c(m,m,lag.max+1))

      # a[0] * b[0] = I
      b0 = try(solve(matrix(a[,,1], nrow = m, ncol = m)))
      if (inherits(b0, 'try-error')) {
        stop('impulse response is not invertible (lag zero coefficient is singular)')
      }
      b[,,1] = b0
      for (i in iseq(1,lag.max)) {
        # a[i] * b[0] + ... + a[0] b[i] = 0
        for (j in (1:i)) {
          b[,,i+1] = b[,,i+1] + matrix(a[,,j + 1], nrow = m, ncol = m) %*%
            matrix(b[,,i - j + 1], nrow = m, ncol = m)
        }
        b[,,i+1] = -b0 %*% matrix(b[,,i+1], nrow = m, ncol = m)
      }

      # convert to pseries object
      class(b) = c('pseries','ratm')

      e = b
      for (i in iseq(2,abs(a2))) {
        e = e %r% b
      }
      return(e)
    }

    if (cl1 == 'zvalues')  {
      z = attr(e1, 'z')
      e1 = unclass(e1)
      # compute inverse
      for (i in (1:length(z))) {
        ifr = try(solve(matrix(e1[,,i], nrow = d1[1], ncol = d1[1])), silent = TRUE)
        if (inherits(ifr, 'try-error')) {
          ifr = matrix(NA_real_, nrow = d1[1], ncol=d1[1])
        }
        e1[,,i] = ifr
      }
      e1 = structure(e1, z = z, class = c('zvalues', 'ratm'))

      e = e1
      for (i in iseq(2,abs(a2))) {
        e = e %r% e1
      }
      return(e)
    }

    # this should not happen!
    stop('unsupported class: "',cl1,'"')
  }


  # elementwise operations '*', '+', '-', '%/%', '%%' ################################

  # make sure that both arguments have the same class!
  out = upgrade_objects(force = TRUE, e1, e2)
  e1 = out[[1]]
  e2 = out[[2]]

  d1 = unname(dim(e1))
  cl1 = class(e1)[1]
  # not needed according to RStudio: gr1 = match(cl1, c('polm','lmfd','rmfd','stsp','pseries','zvalues'), nomatch = 0)
  d2 = unname(dim(e2))
  # not needed according to RStudio: cl2 = class(e2)[1]
  # not needed according to RStudio: gr2 = match(cl1, c('polm','lmfd','rmfd','stsp','pseries','zvalues'), nomatch = 0)

  if (d1[1]*d1[2] == 1) {
    # e1 is a scalar
    e1 = inflate_object(e1, m = d2[1], n = d2[2])
    d1 = unname(dim(e1))
  }
  if (d2[1]*d2[2] == 1) {
    # e2 is a scalar
    e2 = inflate_object(e2, m = d1[1], n = d1[2])
    d2 = unname(dim(e2))
  }

  if (any(d1[1:2] != d2[1:2]))  {
    stop('the rational matrices e1, e2 are not compatible (nrow(e1) != nrow(e2) or ncol(e1) != ncol(e2))')
  }

  # __elementwise multiplication '*' ################################
  if (.Generic == '*') {
    # elementwise addition/substraction

    if (cl1 == 'polm') {
      e = polm(emult_poly(unclass(e1), unclass(e2)))
      e = prune(e, tol = 0)
      return(e)
    }

    if (cl1 == 'lpolm') {
      attr_e1 = attributes(e1)
      attr_e2 = attributes(e2)
      e = lpolm(emult_poly(unclass(e1), unclass(e2)), min_deg = attr_e1$min_deg + attr_e2$min_deg)

      attr_e = attributes(e)
      e = prune(polm(unclass(e)), tol = 0)

      e = lpolm(array(c(e), dim = attr_e$dim), min_deg = attr_e$min_deg)
      return(e)
    }


    if (cl1 == 'stsp') {
      e = emult_stsp(e1, e2)
      return(e)
    }

    if (cl1 == 'pseries') {
      e = emult_pseries(unclass(e1), unclass(e2))
      class(e) = c('pseries','ratm')
      return(e)
    }

    if (cl1 == 'zvalues') {
      z1 = attr(e1,'z')
      z2 = attr(e2,'z')
      if (!isTRUE(all.equal(z1, z2))) {
        stop('the complex points z of the two frequency responses do not coincide')
      }
      e1 = unclass(e1)
      e2 = unclass(e2)
      e = array(0, dim = c(d1[1], d2[2], length(z1)))
      for (i in (1:length(z1))) {
        e[,,i] = matrix(e1[,,i], nrow = d1[1], ncol = d1[2]) * matrix(e2[,,i], nrow = d2[1], ncol = d2[2])
      }

      e = structure(e, z = z1, class = c('zvalues', 'ratm'))
      return(e)
    }

  } # elementwise multiplication '*'

  # __elementwise addition/substraction '+', '-' ################################
  if ((.Generic == '+') || (.Generic == '-')) {
    # elementwise addition/substraction

    if (.Generic == '-') e2 = -e2

    if (cl1 == 'polm') {
      # polynomial matrices
      e1 = unclass(e1)
      e2 = unclass(e2)
      e = array(0, dim = c(d1[1], d1[2], max(d1[3], d2[3])+1))
      if (d1[3]>=0) e[,,1:(d1[3]+1)] = e1
      if (d2[3]>=0) e[,,1:(d2[3]+1)] = e[,,1:(d2[3]+1),drop = FALSE] + e2
      return(polm(e))
    }

    if (cl1 == 'lpolm') {
      # Laurent polynomial matrices
      attr_e1 = attributes(e1)
      attr_e2 = attributes(e2)

      dim1 = attr_e1$dim
      dim2 = attr_e2$dim

      q1 = attr_e1$min_deg
      p1 = dim1[3]-1+q1

      q2 = attr_e2$min_deg
      p2 = dim2[3]-1+q2

      min_q = min(q1, q2)
      max_p = max(p1, p2)

      e1 = unclass(e1)
      e2 = unclass(e2)

      min_deg_e = min(attr_e1$min_deg, attr_e2$min_deg)

      e = dbind(d = 3,
                array(0, dim = c(dim1[1], dim1[2], -min_q+q1)),
                e1,
                array(0, dim = c(dim1[1], dim1[2], max_p-p1))) +
        dbind(d = 3,
              array(0, dim = c(dim2[1], dim2[2], -min_q+q2)),
              e2,
              array(0, dim = c(dim2[1], dim1[2], max_p-p2)))

      return(lpolm(e, min_deg = min_deg_e))
    }

    if (cl1 == 'stsp') {
      # statespace realization
      e1 = unclass(e1)
      e2 = unclass(e2)

      A = bdiag(e1[iseq(1,       d1[3]),       iseq(1,       d1[3])        , drop = FALSE],
                e2[iseq(1,       d2[3]),       iseq(1,       d2[3])        , drop = FALSE])
      B = rbind(e1[iseq(1,       d1[3]),       iseq(d1[3]+1, d1[3]+d1[2])  , drop = FALSE],
                e2[iseq(1,       d2[3]),       iseq(d2[3]+1, d2[3]+d2[2])  , drop = FALSE])
      C = cbind(e1[iseq(d1[3]+1, d1[3]+d1[1]), iseq(1,       d1[3])        , drop = FALSE],
                e2[iseq(d2[3]+1, d2[3]+d2[1]), iseq(1,       d2[3])        , drop = FALSE])
      D =       e1[iseq(d1[3]+1, d1[3]+d1[1]), iseq(d1[3]+1, d1[3]+d1[2])  , drop = FALSE]     +
        e2[iseq(d2[3]+1, d2[3]+d2[1]), iseq(d2[3]+1, d2[3]+d2[2])  , drop = FALSE]
      return(stsp(A = A, B = B, C = C, D = D))
    }

    if (cl1 == 'pseries') {
      # print(attributes(e1))
      # print(attributes(e2))
      # impulse response
      e1 = unclass(e1)
      e2 = unclass(e2)
      if (d1[3] <= d2[3]) {
        e1 = e1 + e2[,,iseq(1, d1[3]+1), drop = FALSE]
        e1 = structure(e1, class = c('pseries', 'ratm'))
        return(e1)
      }
      e2 = e2 + e1[,,iseq(1, d2[3]+1), drop = FALSE]
      e2 = structure(e2, class = c('pseries', 'ratm'))
      return(e2)
    }

    if (cl1 == 'zvalues') {
      # frequency response
      z1 = attr(e1,'z')
      z2 = attr(e2,'z')
      if (!isTRUE(all.equal(z1, z2))) {
        stop('the complex points z of the two frequency responses do not coincide')
      }
      e1 = unclass(e1)
      e2 = unclass(e2)
      e1 = e1 + e2
      e1 = structure(e1, z = z1, class = c('zvalues', 'ratm'))
      return(e1)
    }

  }  # elementwise addition/substraction '+', '-'


  # __elementwise (polynomial) division #########################################
  if (.Generic == '%/%') {
    if ( cl1 != 'polm' ) {
      stop('elementwise (polynomial) divsision "%/%" is only implemented for "polm" objects')
    }

    e = polm(array(0, dim = c(d1[1], d1[2], 1)))
    if (d1[1]*d1[2] == 0) {
      return( e )
    }

    a1 = unclass(e1)
    a2 = unclass(e2)
    for ( i in (1:d1[1]) ) {
      for (j in (1:d1[2])) {
        #        print(str(a))
        #        print(polm_div(a1[i,j,], a2[i,j,]))
        e[i,j] = poly_div(a1[i,j,], a2[i,j,])
      }
    }
    # skip leading zero coefficient matrices
    e = prune(e, tol = 0)
    return( e )
  }

  # __elementwise (polynomial) remainder #########################################
  if (.Generic == '%%') {
    if ( cl1 != 'polm' ) {
      stop('elementwise (polynomial) remainder "%%" is only implemented for "polm" objects')
    }

    e = polm(array(0, dim = c(d1[1], d1[2], 1)))

    if (d1[1]*d1[2] == 0) {
      return( e )
    }

    a1 = unclass(e1)
    a2 = unclass(e2)
    for ( i in (1:d1[1]) ) {
      for (j in (1:d1[2])) {
        #        print(str(a))
        #        print(polm_div(a1[i,j,], a2[i,j,]))
        e[i,j] = poly_rem(a1[i,j,], a2[i,j,])
      }
    }
    # skip leading zero coefficient matrices
    e = prune(e, tol = 0)
    return( e )
  }

  stop('unsupported operator: "',.Generic,'"')
}



# DELETE: Not used anywhere ####

# internal function
# create a diagonal rational matrix from a scalar
# brauche ich das noch ????
diag_object = function(e, d) {
  cl = class(e)[1]

  if (cl == 'polm') {
    e = unclass(prune(e, tol = 0))
    e = as.vector(e)
    if ((length(e) == 0) || (d == 0)) {
      return(polm(matrix(0, nrow = d, ncol = d)))
    }
    a = matrix(0, nrow = length(e), ncol = d*d)
    a[, diag(matrix(1:(d*d), nrow = d, ncol = d))] = e
    a = t(a)
    dim(a) = c(d,d,length(e))
    return(polm(a))
  }

  if (cl == 'stsp') {
    if (d == 0) {
      return(stsp(A = matrix(0, nrow = 0, ncol = 0),
                  B = matrix(0, nrow = 0, ncol = 0),
                  C = matrix(0, nrow = 0, ncol = 0),
                  D = matrix(0, nrow = 0, ncol = 0)))
    }
    A = e$A
    B = e$B
    C = e$C
    D = e$D
    A = do.call(bdiag, args = rep(list(A), d))
    B = do.call(bdiag, args = rep(list(B), d))
    C = do.call(bdiag, args = rep(list(C), d))
    D = diag(x = D[1,1], nrow = d, ncol = d)
    return(stsp(A, B, C, D))
  }

  if (cl == 'zvalues') {
    z = attr(e, 'z')
    if ((length(z) == 0) || (d == 0)) {
      a = array(0, dim = c(d,d,length(z)))
      a = structure(a, z = z, class = c('zvalues', 'ratm'))
      return(a)
    }
    e = unclass(e)
    e = as.vector(e)
    a = matrix(0, nrow = length(z), ncol = d*d)
    a[, diag(matrix(1:(d*d), nrow = d, ncol = d))] = e
    a = t(a)
    dim(a) = c(d,d,length(z))
    a = structure(a, z = z, class = c('zvalues', 'ratm'))
    return(a)
  }

  stop('unsupported class "', cl, '" for diagonalization of scalars')
}

# Conversion between polm() and lpolm() ####

as.polm = function(obj, ...){
  UseMethod("as.polm", obj)
}

as.polm.lpolm = function(obj, ...){
  obj = prune(obj)
  min_deg = attr(obj, "min_deg")
  stopifnot("The *min_deg* attribute needs to be non-negative for coercion to polm-obj! Use get_bwd() for discarding negative powers." = min_deg >= 0)
  polm_offset = array(0, dim = c(dim(obj)[1:2], min_deg))
  return(polm(dbind(d = 3, polm_offset, unclass(obj))))
}

as.lpolm = function(obj, ...){
  UseMethod("as.lpolm", obj)
}

as.lpolm.polm = function(obj, ...){
  attr(obj, "min_deg") = 0
  class(obj)[1] = "lpolm"
  return(obj)
}



# as.stsp.____ methods ##################################################

as.stsp = function(obj, ...){
  UseMethod("as.stsp", obj)
}

as.stsp.lpolm = function(obj, ...){
  stop("A lpolm object cannot be coerced to a state space model.")
}

as.stsp.polm = function(obj, ...){
  obj = unclass(obj)
  d = dim(obj)
  m = d[1]
  n = d[2]
  p = d[3] - 1

  if (p < 0) {
    x = stsp(A = matrix(0, nrow = 0, ncol = 0), B = matrix(0, nrow = 0, ncol = n),
             C = matrix(0, nrow = m, ncol = 0), D = matrix(0, nrow = m, ncol = n))
    return(x)
  }
  if (p == 0) {
    x = stsp(A = matrix(0, nrow = 0, ncol = 0), B = matrix(0, nrow = 0, ncol = n),
             C = matrix(0, nrow = m, ncol = 0), D = matrix(obj, nrow = m, ncol = n))
    return(x)
  }
  if (m >= n) {
    x = stsp(A = rbind(matrix(0, nrow = n, ncol = p*n), diag(x = 1, nrow = (p-1)*n, ncol = p*n)),
             B = diag(x = 1, nrow = p*n, ncol = n),
             C = matrix(obj[,,-1], nrow = m, ncol = p*n), D = matrix(obj[,,1], nrow = m, ncol = n))
    return(x)
  }
  B = obj[,,-1,drop = FALSE]
  B = aperm(B, c(1,3,2))
  dim(B) = c(p*m, n)
  x = stsp(A = cbind(matrix(0, nrow = p*m, ncol = m), diag(x = 1, nrow = p*m, ncol = (p-1)*m)),
           B = B, C = diag(x = 1, nrow = m, ncol = p*m), D = matrix(obj[,,1], nrow = m, ncol = n))
  return(x)
}

as.stsp.lmfd = function(obj, ...){
  d = attr(obj, 'order')
  m = d[1]
  n = d[2]
  p = d[3]
  q = d[4]

  # note for a valid lmfd object, m > 0, p >= 0 must hold!
  if ((m*(p+1)) == 0) stop('input object is not a valid "lmfd" object')

  if ( (n*(q+1)) == 0 ) {
    x = stsp(A = matrix(0, nrow = 0, ncol = 0), B = matrix(0, nrow = 0, ncol = n),
             C = matrix(0, nrow = m, ncol = 0), D = matrix(0, nrow = m, ncol = n))
    return(x)
  }

  ab = unclass(obj)
  a = ab[,1:(m*(p+1)), drop = FALSE]
  b = ab[,(m*(p+1)+1):(m*(p+1)+n*(q+1)), drop = FALSE]

  # check a(z)
  a0 = matrix(a[, 1:m, drop = FALSE], nrow = m, ncol = m)
  junk = try(solve(a0))
  if (inherits(junk, 'try-error')) stop('left factor "a(0)" is not invertible')

  if ((p == 0) && (q == 0)) {
    # static system
    x = stsp(A = matrix(0, nrow = 0, ncol = 0), B = matrix(0, nrow = 0, ncol = n),
             C = matrix(0, nrow = m, ncol = 0), D = solve(a0, matrix(b, nrow = m, ncol = n)))
    return(x)
  }

  # a[,,i] -> a0^{-1} a[,,i], convert to matrix and append zeroes if p < q
  if (p > 0) {
    a = solve(a0, a[, (m+1):(m*(p+1)), drop = FALSE])
  } else {
    a = matrix(0, nrow = m, ncol = 0)
  }
  if (p < q) a = cbind(a, matrix(0, nrow = m, ncol = (q-p)*m))
  p = max(p,q)

  # compute impulse response
  # this is not very efficient,
  # e.g. the scaling of the coefficients a0^(-1)a[,,i] and a0^(-1)b[,,i] is done twice
  k = unclass(pseries(obj, lag.max = p))

  A = rbind(-a, diag(x = 1, nrow = m*(p-1), ncol = m*p))
  D = matrix(k[,,1], nrow = m, ncol = n)
  C = cbind(matrix(0, nrow = m, ncol = (p-1)*m), diag(m))
  B = aperm(k[,,(p+1):2,drop = FALSE], c(1,3,2))
  dim(B) = c(p*m, n)
  x = stsp(A = A, B = B, C = C, D = D)
  return(x)
}


as.stsp.rmfd = function(obj, ...){
  d = attr(obj, 'order')
  m = d[1] # output dimension
  n = d[2] # input dimension
  p = d[3] # degree of c(z)
  q = d[4] # degree of d(z)

  # otherwise c(z) is not invertible
  stopifnot("Input object is not a valid rmfd object" = (m*(p+1)) > 0)

  # if d(z) is identically zero or has no column, return early
  if ( (n*(q+1)) == 0 ) {
    x = stsp(A = matrix(0, nrow = 0, ncol = 0), B = matrix(0, nrow = 0, ncol = n),
             C = matrix(0, nrow = m, ncol = 0), D = matrix(0, nrow = m, ncol = n))
    return(x)
  }

  # c(0) must be the identity matrix. Note that this trafo only changes the covariance matrix of the inputs
  c0 = matrix(c(unclass(obj$c))[1:(n^2)], nrow = n, ncol = n)
  c0inv = tryCatch(solve(c0),
                   error = function(cnd) stop(' "c(0)" is not invertible'))
  cd = unclass(obj) %*% c0inv

  # Separate c(z) and d(z)
  c = cd[1:(n*(p+1)),, drop = FALSE]
  d = cd[(n*(p+1)+1):(n*(p+1)+m*(q+1)),, drop = FALSE]

  # Case 1: Static System
  if ((p == 0) && (q == 0)) {
    x = stsp(A = matrix(0, nrow = 0, ncol = 0), B = matrix(0, nrow = 0, ncol = n),
             C = matrix(0, nrow = m, ncol = 0), D = matrix(d, nrow = m, ncol = n))
    return(x)
  }

  # Case 2: (p>0 or q>0): Construct state space system
  if (p < q+1){
    c = rbind(c,
              matrix(0, nrow = (q+1-p)*n, ncol = n))
  }
  p = max(p,q+1)

  # Transpose and reshuffle c(z) such that the coefficients are in a wide matrix, then create stsp matrices (A,B)
  c = t(c)[,-(1:n), drop = FALSE] %>% array(dim = c(n,n,p)) %>% aperm(perm = c(2,1,3)) %>% matrix(nrow = n)

  A = rbind(-c,
            diag(x = 1, nrow = n*(p-1), ncol = n*p))
  B = rbind(c0inv,
            matrix(0, nrow = n*(p-1), ncol = n))

  # Same for d(z), and add zeros to d(z) if p>q+1, then create stsp matrices (C,D)
  d = t(d) %>% array(dim = c(n,m,q+1)) %>% aperm(perm = c(2,1,3)) %>% matrix(nrow = m)
  C = cbind(d, matrix(0, nrow = m, ncol = n*(p-q-1)))
  D = C %*% B
  C = C %*% A

  # Create output
  x = stsp(A = A, B = B,
           C = C, D = D)
  return(x)
}




pseries2stsp = function(obj, method = c('balanced', 'echelon'),
                        Hsize = NULL, s = NULL, nu = NULL,
                        tol = sqrt(.Machine$double.eps), Wrow = NULL, Wcol = NULL) {
  # no input checks

  method = match.arg(method)

  # construct Hankel matrix with the helper function pseries2hankel
  H = try(pseries2hankel(obj, Hsize = Hsize))
  if (inherits(H, 'try-error')) {
    stop('computation of Hankel matrix failed')
  }

  k0 = attr(H, 'k0')
  d = attr(H, 'order')
  m = d[1]
  n = d[2]
  f = d[3]
  p = d[4]

  # take care of the case of an empty impulse response (m*n=0)
  if ((m*n) == 0) {
    s = 0
    Xs = stsp(A = matrix(0, nrow = s, ncol = s), B = matrix(0, nrow = s, ncol = n),
              C = matrix(0, nrow = m, ncol = s), D = matrix(0, nrow = m, ncol = n))
    return(list(Xs = Xs, Hsv = numeric(0), nu = integer(0)))
  }

  # compute statespace model in "balanced" form ###############
  if (method == 'balanced') {
    if (!is.null(Wrow)) {
      H = Wrow %*% H
    }
    if (!is.null(Wcol)) {
      H = H %*% t(Wcol)
    }

    svd.H = svd(H)
    Hsv = svd.H$d    # singular values of (weighted) Hankel matrix

    if (is.null(s)) {
      # determine state dimension from singular values
      s = ifelse(svd.H$d[1]>.Machine$double.eps, sum(svd.H$d >= (tol*svd.H$d[1])),0)
    }

    if (s>0) {
      if (s > m*(f-1)) {
        stop('number of block rows of "H" is too small for the (desired) state dimension "s"!')
      }
      sv2 = sqrt(Hsv[1:s])
      # (approximately) factorize H as H = U V
      U = svd.H$u[,1:s,drop = FALSE] * matrix(sv2, nrow = nrow(H), ncol = s, byrow = TRUE)
      V = t( svd.H$v[,1:s,drop = FALSE] * matrix(sv2, nrow = ncol(H), ncol = s, byrow = TRUE) )
      if (!is.null(Wrow)) {
        U = solve(Wrow, U)
      }
      if (!is.null(Wcol)) {
        V = t( solve(Wcol, t(V)) )
      }
      C = U[1:m,,drop = FALSE]
      B = V[,1:n,drop = FALSE]
      A = stats::lsfit(U[1:(m*(f-1)), ,drop = FALSE],
                       U[(m+1):(m*f),,drop = FALSE], intercept = FALSE)$coef
      A = unname(A)
    } else {
      A = matrix(0, nrow=s,  ncol=s)
      B = matrix(0, nrow=s,  ncol=n)
      C = matrix(0, nrow=m, ncol=s)
    }
    D = k0

    Xs = stsp(A = A, B = B, C = C, D = D)

    return( list(Xs = Xs, Hsv = Hsv, nu = NULL) )
  } # method == balanced

  # compute statespace model in "echelon" form ###############

  # if (n < m) {
  #   stop('method "echelon" for the case (n < m) is not yet implemented!')
  # }

  # compute Kronecker indices
  if (is.null(nu)) {
    nu = try(hankel2nu(H, tol = tol))
    if (inherits(nu, 'try-error')) stop('computation of Kronecker indices failed')
  }
  # check nu
  nu = as.vector(as.integer(nu))
  if ( (length(nu) != m) || (min(nu) < 0) || (max(nu) > (f-1)) ) {
    stop('Kronecker indices are not compatible with the impulse response')
  }
  basis = nu2basis(nu)
  s = length(basis) # state space dimension

  # consider the transposed Hankel matrix!
  H = t(H)

  if (s > 0) {
    AC = matrix(0, nrow = s+m, ncol = s)
    dependent = c(basis + m, 1:m)
    # cat('basis', basis,'\n')
    for (j in (1:length(dependent))) {
      i = dependent[j]
      if (min(abs(i - basis)) == 0) {
        # row i is in the basis!
        k = which(i == basis)
        AC[j, k] = 1
      } else {
        # explain row i by the previous basis rows
        k = which(basis < i)
        AC[j, k] = stats::lsfit(H[ , k, drop=FALSE], H[ , i], intercept=FALSE)$coef
      }
      # cat('j=', j, 'i=', i, 'k=', k, 'basis[k]=', basis[k], '\n')
    }
    B = t(H[(1:n), basis, drop=FALSE])
    A = AC[1:s,,drop = FALSE]
    C = AC[(s+1):(s+m), , drop = FALSE]
  } else {
    A = matrix(0, nrow=s,  ncol=s)
    B = matrix(0, nrow=s,  ncol=n)
    C = matrix(0, nrow=m, ncol=s)
  }
  D = k0

  Xs = stsp(A = A, B = B, C = C, D = D)

  return( list(Xs = Xs, Hsv = NULL, nu = nu) )
}


as.stsp.pseries = function(obj, method = c('balanced','echelon'), ...){
  out = pseries2stsp(obj, method = method)
  return(out$Xs)
}

# as.lmfd.____ methods ##################################################


pseries2lmfd = function(obj, Hsize = NULL, nu = NULL, tol = sqrt(.Machine$double.eps)) {


  # construct Hankel matrix with the helper function pseries2hankel
  H = try(pseries2hankel(obj, Hsize = Hsize))
  if (inherits(H, 'try-error')) {
    stop('computation of Hankel matrix failed')
  }

  k0 = attr(H, 'k0')
  d = attr(H, 'order')
  m = d[1]
  n = d[2]
  f = d[3]
  # p = d[4]

  if (m == 0) stop('the number of rows (m) must be positive!')

  # compute Kronecker indices
  if (is.null(nu)) {
    nu = try(hankel2nu(H, tol = tol))
    if (inherits(nu, 'try-error')) stop('computation of Kronecker indices failed')
  }
  # check nu
  nu = as.vector(as.integer(nu))
  if ( (length(nu) != m) || (min(nu) < 0) || (max(nu) > (f-1)) ) {
    stop('Kronecker indices are not compatible with the impulse response')
  }

  k = unclass(obj)
  H = t(H) # use the transposed of H in the following!

  p = max(nu)

  # a(z),b(z) have degree zero => c(z) is constant
  if (p == 0) {
    Xl = lmfd(b = k0)
    return(list(Xl = Xl, nu = nu))
  }

  # Note that the block diagonal is zero because
  # the element R[,,1] in the argument R of btoeplitz()
  # overwrites the block diagonal element from argument C
  # (corresponding to the zero lag coefficient of the transfer function)
  Tk = btoeplitz(R = array(0, dim = c(m, n, p)), C = k[, , (1:(p+1)), drop=FALSE])
  # print(Tk)

  basis = nu2basis(nu)
  # index of the "dependent" rows of H
  dependent = m*nu + (1:m)
  # matrices with coefficients in reversed order
  a = matrix(0, nrow = m, ncol = m*(p+1))  # [a[p], a[p-1], ..., a[0]]
  # note b relates to b(z) - b(0) = b(z) - a(0)*k(0)
  b = matrix(0, nrow = m, ncol = n*(p+1))  # [b[p], b[p-1], ..., b[0]]
  for (i in (1:m)) {
    j = basis[basis < dependent[i]]
    # cat('i=', i,', dependent=', dependent[i], ', j=', j,'\n',sep = ' ')
    if (length(j) > 0) {
      ai = stats::lsfit(H[ , j,drop = FALSE], H[ , dependent[i]], intercept = FALSE)$coef
      a[i, j + m*(p-nu[i])] = -ai
    }
    a[i, dependent[i] + m*(p-nu[i])] = 1
    # print(a)
    # note b(0) = 0
    # cat(i,':',nu[i],':',p,':',iseq(1 + n*(p-nu[i]), n*p),'\n')
    b[i, iseq(1 + n*(p-nu[i]), n*p)] =
      a[i,] %*% Tk[, iseq(1 + n*(p-nu[i]), n*p), drop = FALSE]
    # print(b)
  }
  # print(t(H))
  # print(a)
  # print(a %*% t(H)[1:ncol(a),,drop = FALSE])
  # print(b)
  dim(a) = c(m, m, p+1)
  dim(b) = c(m, n, p+1)

  # ak0 <=> a(z)*k(0)
  ak0 = aperm(a, c(1,3, 2))
  dim(ak0) = c(m*(p+1), m)
  ak0 = ak0 %*% k0
  dim(ak0) = c(m, p+1, n)
  ak0 = aperm(ak0, c(1, 3, 2))

  # b -> b + a(z)*k(0)
  b = b + ak0

  # reshuffle
  a = a[,,(p+1):1]
  b = b[,,(p+1):1]

  Xl = lmfd(polm(a), polm(b))
  return(list(Xl = Xl, nu = nu))
}


as.lmfd = function(obj, method, ...){
  UseMethod("as.lmfd", obj)
}

as.lmfd.pseries = function(obj, method, ...){
  return(pseries2lmfd(obj)$Xl)
}

# as.rmfd.____ methods ####

pseries2rmfd = function(obj, Hsize = NULL, mu = NULL, tol = sqrt(.Machine$double.eps)) {

  # Construct Hankel matrix with the helper function pseries2hankel ####
  H = try(pseries2hankel(obj, Hsize = Hsize))
  if (inherits(H, 'try-error')) {
    stop('computation of Hankel matrix failed')
  }

  # __ Dimensions of Hankel ####
  k0 = attr(H, 'k0')
  H_order = attr(H, 'order')
  dim_out = H_order[1]
  dim_in = H_order[2]
  block_rows = H_order[3]
  block_cols = H_order[4] # this was commented out for pseries2lmfd, here we need it though

  if (dim_out == 0) stop('The number of rows (m) must be positive!')

  # Compute right-Kronecker indices mu ####
  if (is.null(mu)) {
    mu = try(hankel2mu(H, tol = tol))
    if (inherits(mu, 'try-error')){
      stop('Computation of Kronecker indices failed')
    }
  }
  # check mu
  mu = as.vector(as.integer(mu))
  if ( (length(mu) != dim_in) || (min(mu) < 0) || (max(mu) > (block_cols-1)) ) {
    stop('Kronecker indices are not compatible with the impulse response')
  }

  # Transfer function as array (i.e. power series coefficients) ####
  k = unclass(obj)
  max_mu = max(mu) # difficulty in notation: p is used as deg(c(z)) or deg(a(z)), but also as number of block columns of the Hankel matrix

  # Special case: Zero poly (c(z), d(z) have degree zero => c(z) is constant ) ####
  if (max_mu == 0) {
    Xr = rmfd(d = k0)
    return(list(Xr = Xr, mu = mu))
  }

  # Needed to determine the coefficients in d(z) ####
  # Note how the coefficients are stacked when using this Toeplitz matrix! ([\tilde{d}[0]', \tilde{d}[1]', ..., \tilde{d}[p]']')
  Tk = btoeplitz(R = array(c(rep(0,dim_out*dim_in),
                           k[,,2:(max_mu+1)]), dim = c(dim_out, dim_in, max_mu+1)),
                 C = array(rep(0, dim_out*dim_in*(max_mu+1)), dim =c(dim_out, dim_in, max_mu+1)))

  # Which columns of Hankel span its column space? ####
  basis = nu2basis(mu) # nu2basis is okay also for right-Kronecker indices. No need to create new function mu2basis

  # First dependent column of Hankel for each variable (index of the "dependent" columns of H) ####
  dependent = dim_in*mu + (1:dim_in)

  # Tilde coefficient matrices ####
  # \tilde{c}(z) = \tilde{c_0} + \tilde{c_1} z + ... for describing \tilde{k}(z) = k_1 z^(-1) + k_2 z^(-2) + ... (see Hannan/Deistler Chapter 2.4 and 2.5)

  # [\tilde{c}[0]', \tilde{c}[1]', ..., \tilde{c}[p]']' in echelon form, but the columns get shifted immediately with z^(p-mu_i) (in the for-loop)!
  c = matrix(0, nrow = dim_in*(max_mu+1), ncol = dim_in)
  # [\tilde{d}[0]', \tilde{d}[1]', ..., \tilde{d}[p]']': also shifted in the for-loop below
  d = matrix(0, nrow = dim_out*(max_mu+1), ncol = dim_in)

  # Obtain the coefficients for each (dependent) column of Hankel ####
  for (ix_var in (1:dim_in)) {

    # indices of basis columns of Hankel to the left of the first dependent column pertaining to variable ix_var
    ix_basis = basis[basis < dependent[ix_var]]

    # regression coefficients
    if (length(ix_basis) > 0) {
      ci = stats::lsfit(H[ , ix_basis,drop = FALSE], H[ , dependent[ix_var]], intercept = FALSE)$coef #lsfit(x,y)
      c[ix_basis, ix_var] = -ci
    }
    # Coefficient of dependent column corresponding to variable ix_var is set to one
    c[dependent[ix_var], ix_var] = 1

    # Obtain \tilde{d}(z) from the Toeplitz matrix (corresponding to non-negative coefficients in comparison)
    # and the corresponding column in the matrix representation of \tilde{c}(z)
    d[, ix_var] =  Tk %*% c[, ix_var, drop = FALSE]

    # Multiply z^(max_mu - mu[ix_var]) on column "ix_var" of matrix representation of \tilde{c}(z) and \tilde{d}(z)
    # such that \tilde{c}_p has ones on the diagonal and is upper triangular
    c[, ix_var] = c(rep(0, dim_in*(max_mu - mu[ix_var])),
                    c[1:(dim_in*(mu[ix_var]+1)), ix_var])
    d[, ix_var] = c(rep(0, dim_out*(max_mu - mu[ix_var])),
                    d[1:(dim_out*(mu[ix_var]+1)), ix_var])

  }

  # Transform \tilde{c}(z) to array ####
  c = t(c)
  dim(c) = c(dim_in, dim_in, max_mu+1)
  c = aperm(c, perm = c(2,1,3))

  # k0c <=> k(0)*c(z) (to be added to \tilde{d}(z)) ####
  k0c = purrr::map(1:(max_mu+1), ~ k0 %*% c[,,.x]) %>% unlist() %>% array(dim = c(dim_out, dim_in, max_mu+1))

  # Transform \tilde{d}(z) to array ####
  d = t(d)
  dim(d) = c(dim_in, dim_out, max_mu+1) # (dim_in, dim_out) because we work here with tranposed!
  d = aperm(d, perm = c(2,1,3))

  # Add k_0 * \tilde(c)(z) to \tilde{d}(z) ####
  d = d + k0c

  # Obtain representation in backward shift (instead of forward shift) ####
  c = c[,,(max_mu+1):1, drop = FALSE]
  d = d[,,(max_mu+1):1, drop = FALSE]

  Xr = rmfd(polm(c), polm(d))
  return(list(Xr = Xr, mu = mu))
}


as.rmfd = function(obj, method, ...){
  UseMethod("as.rmfd", obj)
}

as.rmfd.pseries = function(obj, method, ...){
  return(pseries2rmfd(obj)$Xr)
}
lpolm = function(a, min_deg = 0) {

  stopifnot("Only array, matrix, or vector objects can be supplied to this constructor." = any(class(a) %in% c("matrix", "array")) || is.vector(a))

  stopifnot('Input "a" must be a numeric or complex vector/matrix/array!' = (is.numeric(a) || is.complex(a)))

  if (is.vector(a)) {
    dim(a) = c(1,1,length(a))
  }
  if (is.matrix(a)) {
    dim(a) = c(dim(a),1)
  }

  stopifnot('could not coerce input parameter "a" to a valid 3-D array!' = is.array(a) && (length(dim(a)) == 3),
            "Minimal degree 'min_deg' must be set, numeric, and of length 1!" = !is.null(min_deg) && is.numeric(min_deg) && length(min_deg) == 1)

  # achtung
  min_deg = as.integer(min_deg)[1]

  # Initially, lpolm object were implicitly coerced to polm object when deg_min >= 0.
  # This was not optimal in terms of type consistency.
  # It is now necessary to coerce such lpolm objects to polm objects explicitly with as.polm.lpolm()
  # The code below is used in as.polm.lpolm()
  # if (min_deg >= 0) {
  #   polm_offset = array(0, dim = c(dim(a)[1:2], min_deg))
  #   return(polm(dbind(d = 3, polm_offset, a)))
  # }

  return(structure(a, class = c("lpolm", 'ratm'), min_deg = min_deg))
}


polm = function(a) {

  stopifnot("Only array, matrix, or vector objects can be supplied to this constructor." = any(class(a) %in% c("matrix", "array")) || is.vector(a))

  if (!(is.numeric(a) || is.complex(a))) {
    stop('input "a" must be a numeric or complex vector/matrix/array!')
  }
  if (is.vector(a)) {
    dim(a) = c(1,1,length(a))
  }
  if (is.matrix(a)) {
    dim(a) = c(dim(a),1)
  }

  if ( (!is.array(a)) || (length(dim(a)) !=3) ) {
    stop('could not coerce input parameter "a" to a valid 3-D array!')
  }
  class(a) = c('polm','ratm')
  return(a)
}

lmfd = function(a, b) {
  if (missing(a) && missing(b)) {
    stop('no arguments have been provided')
  }
  if (!missing(a)) {
    if (!inherits(a,'polm')) {
      a = try(polm(a))
      if (inherits(a, 'try-error')) {
        stop('could not coerce "a" to a "polm" object!')
      }
    }
    dim_a = dim(unclass(a))
    if ((dim_a[1] == 0) || (dim_a[1] != dim_a[2]) || (dim_a[3] == 0)) {
      stop('"a" must represent a square polynomial matrix with degree >= 0')
    }
  }
  if (!missing(b)) {
    if (!inherits(b,'polm')) {
      b = try(polm(b))
      if (inherits(b, 'try-error')) {
        stop('could not coerce "b" to a "polm" object!')
      }
    }
    dim_b = dim(unclass(b))
  }
  if (missing(b)) {
    b = polm(diag(dim_a[2]))
    dim_b = dim(unclass(b))
  }
  if (missing(a)) {
    a = polm(diag(dim_b[1]))
    dim_a = dim(unclass(a))
  }
  if (dim_a[2] != dim_b[1]) {
    stop('the arguments "a", "b" are not compatible')
  }

  c = matrix(0, nrow = dim_b[1], ncol = dim_a[2]*dim_a[3] + dim_b[2]*dim_b[3])
  if ((dim_a[2]*dim_a[3]) > 0) c[, 1:(dim_a[2]*dim_a[3])] = a
  if ((dim_b[2]*dim_b[3]) > 0) c[, (dim_a[2]*dim_a[3]+1):(dim_a[2]*dim_a[3] + dim_b[2]*dim_b[3])] = b

  c = structure(c, order = as.integer(c(dim_b[1], dim_b[2], dim_a[3]-1, dim_b[3] - 1)),
                class = c('lmfd','ratm'))
  return(c)
}


rmfd = function(c = NULL, d = NULL) {

  # Check inputs: At least one polm object must be non-empty ####

  if (is.null(c) && is.null(d)) {
    stop('At least one of c(z) or d(z) needs to be provided.')
  }

  # Convert array to polm objects ####
  if (!is.null(c)) {
    if (!inherits(c,'polm')) {
      c = try(polm(c))
      if (inherits(c, 'try-error')) {
        stop('Could not coerce "c" to a "polm" object!')
      }
    }
    dim_c = dim(unclass(c))
    if ((dim_c[1] == 0) || (dim_c[1] != dim_c[2]) || (dim_c[3] == 0)) {
      stop('"c" must represent a square polynomial matrix with degree >= 0')
    }
  }

  if (!is.null(d)) {
    if (!inherits(d,'polm')) {
      d = try(polm(d))
      if (inherits(d, 'try-error')) {
        stop('Could not coerce "d" to a "polm" object!')
      }
    }
    dim_d = dim(unclass(d))
  }

  # If one polm object is NULL, set to identity matrix ####

  if (is.null(d)) {
    d = polm(diag(dim_c[2]))
    dim_d = dim(unclass(d))
  }

  if (is.null(c)) {
    c = polm(diag(dim_d[2]))
    dim_c = dim(unclass(c))
  }

  # Check input: Dimensions ####

  if (dim_c[1] != dim_d[2]) {
    stop('The dimensions of "c" and "d" are not compatible.')
  }

  h = matrix(0,
             nrow = dim_c[2], # number of cols!
             ncol = dim_c[1]*dim_c[3] + dim_d[1]*dim_d[3]) # (cols of c(z))* (lag.max + 1) + (cols of d(z))* (lag.max + 1)
  if ((dim_c[1]*dim_c[3]) > 0) h[, 1:(dim_c[1]*dim_c[3])] = aperm(c, c(2,1,3))
  if ((dim_d[1]*dim_d[3]) > 0) h[, (dim_c[1]*dim_c[3]+1):(dim_c[1]*dim_c[3] + dim_d[1]*dim_d[3])] = aperm(d, c(2,1,3))
  h = t(h)

  h = structure(h,
                order = as.integer(c(dim_d[1], dim_d[2], dim_c[3]-1, dim_d[3] - 1)),
                class = c('rmfd','ratm'))
  return(h)
}

stsp = function(A, B, C, D) {
  # only D has been given => statespace dimension s = 0
  if (missing(A) && missing(B) && missing(C)) {
    if (missing(D)) stop('no parameters')
    if (!( is.numeric(D) || is.complex(D) )) stop('parameter D is not numeric or complex')
    if ( (is.vector(D)) && (length(D) == 1) ) {
      D = matrix(D, nrow = 1, ncol = 1)
    }
    if (!is.matrix(D)) {
      stop('parameter D must be a numeric or complex matrix')
    }
    m = nrow(D)
    n = ncol(D)
    s = 0
    ABCD = structure(D, order = as.integer(c(m,n,s)),  class = c('stsp','ratm'))
    return(ABCD)
  }

  if (missing(A)) stop('parameter A is missing')
  if ( !( is.numeric(A) || is.complex(A) ) ) stop('parameter A is not numeric or complex')
  if ( is.vector(A) && (length(A) == (as.integer(sqrt(length(A))))^2) )  {
    A = matrix(A, nrow = sqrt(length(A)), ncol = sqrt(length(A)))
  }
  if ( (!is.matrix(A)) || (nrow(A) != ncol(A)) ) stop('parameter A is not a square matrix')
  s = nrow(A)

  if (missing(B)) stop('parameter B is missing')
  if (!( is.numeric(B) || is.complex(B) )) stop('parameter B is not numeric or complex')
  if ( is.vector(B) &&  (s > 0) && ((length(B) %% s) == 0) ) {
    B = matrix(B, nrow = s)
  }
  if ( (!is.matrix(B)) || (nrow(B) != s) ) stop('parameter B is not compatible to A')
  n = ncol(B)

  if (missing(C)) stop('parameter C is missing')
  if (!( is.numeric(C) || is.complex(C) )) stop('parameter C is not numeric or complex')
  if ( is.vector(C) &&  (s > 0) && ((length(C) %% s) == 0) ) {
    C = matrix(C, ncol = s)
  }
  if ( (!is.matrix(C)) || (ncol(C) != s)) stop('parameter C is not compatible to A')
  m = nrow(C)

  if (missing(D)) D = diag(x = 1, nrow = m, ncol = n)
  if (!( is.numeric(D) || is.complex(D) )) stop('parameter D is not numeric or complex')
  if ( is.vector(D) && (length(D) == (m*n)) ) {
    D = matrix(D, nrow = m, ncol = n)
  }
  if ( (!is.matrix(D)) || (nrow(D) != m)  || (ncol(D) != n) ) {
    stop('parameter D is not compatible to B,C')
  }

  ABCD = rbind( cbind(A,B), cbind(C, D) )
  ABCD = structure(ABCD, order = as.integer(c(m,n,s)),  class = c('stsp','ratm'))
  return(ABCD)
}
# derivative.R #######################################
# compute derivatives of rational matrices



derivative = function(obj, ...){
  UseMethod("derivative", obj)
}

derivative.lpolm = function(obj, ...) {
  stop('computation of the derivative is not implemented for "lpolm" objects')
}


derivative.polm = function(obj, ...) {
  obj = unclass(obj)
  dim = dim(obj)
  m = dim[1]
  n = dim[2]
  p = dim[3] - 1

  if (p <= 0) {
    obj = array(0, dim = c(m,n,0))
    class(obj) = c('polm','ratm')
    return(obj)
  }

  obj = obj[,,-1,drop = FALSE]
  k = array(1:p, dim = c(p,m,n))
  k = aperm(k, c(2,3,1))
  obj = obj*k
  class(obj) = c('polm','ratm')
  return(obj)
}


derivative.lmfd = function(obj, ...) {
  stop('computation of the derivative is not implemented for "lmfd" objects')
}

derivative.rmfd = function(obj, ...) {
  stop('computation of the derivative is not implemented for "rmfd" objects')
}

derivative.stsp = function(obj, ...) {
  dim = dim(obj)
  m = dim[1]
  n = dim[2]
  s = dim[3]
  if (min(dim) <= 0) {
    obj = stsp(matrix(0, nrow = 0, ncol = 0), matrix(0, nrow = 0, ncol = n),
               matrix(0, nrow = m, ncol = 0), matrix(0, nrow = m, ncol = n))
    return(obj)
  }
  A = obj$A
  B = obj$B
  C = obj$C

  obj = stsp(rbind( cbind(A, diag(s)), cbind(matrix(0, nrow = s, ncol = s), A) ),
             rbind( B, A %*% B), cbind(C %*%A, C), C %*% B)
  return(obj)
}


derivative.pseries = function(obj, ...) {
  obj = unclass(obj)
  dim = dim(obj)
  m = dim[1]
  n = dim[2]
  p = dim[3] - 1

  if (p <= 0) {
    obj = array(0, dim = c(m,n,1))
    class(obj) = c('pseries','ratm')
    return(obj)
  }

  obj = obj[,,-1,drop = FALSE]
  k = array(1:p, dim = c(p,m,n))
  k = aperm(k, c(2,3,1))
  obj = obj*k
  class(obj) = c('pseries','ratm')
  return(obj)
}


derivative.zvalues = function(obj, ...) {
  stop('computation of the derivative is not implemented for "zvalues" objects')
}
# dim.____ methods ##############################################################

dim.polm = function(x) {
  d = dim(unclass(x))
  d[3] = d[3] - 1
  names(d) = c('m','n','p')
  return(d)
}

dim.lpolm = function(x) {
  d = dim(unclass(x))
  min_deg = attr(x, which = "min_deg")
  d[3] = d[3] - 1 + min_deg
  d = c(d, min_deg)
  names(d) = c('m','n','p', 'min_deg')
  return(d)
}

dim.lmfd = function(x) {
  d = attr(x, 'order')
  names(d) = c('m','n','p','q')
  return(d)
}


dim.rmfd = function(x) {
  d = attr(x, 'order')
  names(d) = c('m','n','p','q')
  return(d)
}



dim.stsp = function(x) {
  d = attr(x, 'order')
  names(d) = c('m','n','s')
  return(d)
}

dim.pseries = function(x) {
  d = dim(unclass(x))
  d[3] = d[3] - 1
  names(d) = c('m','n','lag.max')
  return(d)
}

dim.zvalues = function(x) {
  d = dim(unclass(x))
  names(d) = c('m','n','n.f')
  return(d)
}
# subsetting for rational matrices ###############################################

# helper function:
# check the result of subsetting x[]
# based on the result of subsetting an "ordinary" matrix
extract_matrix_ = function(m, n, n_args, is_missing, i, j) {
  # Write linear indices
  idx = matrix(iseq(1, m*n), nrow = m, ncol = n)

  # Case 1: Everything missing: x[] and x[,]
  if ( ((n_args==1) && is_missing[1]) || all(is_missing) ) return(idx)

  # Case 2: One argument: x[i] (input j is ignored if available (which will not be the case))
  if (n_args == 1) {
    # x[i]
    idx = idx[i]
    if (any(is.na(idx))) stop('index out of bounds')
    idx = matrix(idx, nrow = length(idx), ncol = 1)
    return(idx)
  }

  # Case 3a: Two arguments, first missing: x[,j]
  if (is_missing[1]) {
    # x[,j]
    return(idx[, j, drop = FALSE])
  }

  # Case 3b: Two arguments, second missing: x[i,]
  if (is_missing[2]) {
    # x[i,]
    return(idx[i, , drop = FALSE])
  }

  # Case 3c: Two arguments, none missing: x[i,j]
  return(idx[i, j, drop = FALSE])
}


"[<-.polm" = function(x,i,j,value) {
  names_args = names(sys.call())
  # print(sys.call())
  # print(names_args)
  if (!all(names_args[-length(names_args)] == '')) {
    stop('named dimensions are not supported')
  }
  n_args = nargs() - 2
  # print(n_args)
  is_missing = c(missing(i), missing(j))
  # print(is_missing)

  x = unclass(x)
  dx = dim(x)
  m = dx[1]
  n = dx[2]
  p = dx[3] - 1

  idx = try(extract_matrix_(m, n, n_args, is_missing, i, j), silent = TRUE)
  if (inherits(idx, 'try-error')) stop('index/subscripts out of bounds')
  # print(idx)
  idx = as.vector(idx)
  # print(length(idx))

  if (!inherits(value,'polm')) {
    # coerce right hand side 'value' to polm object
    value = try( polm(value) )
    if (inherits(value, 'try-error')) {
      stop('Could not coerce the right hand side to a polm object!')
    }
  }
  value = unclass(value)
  dv = dim(value)

  # no items to replace: return original object
  if (length(idx) == 0) return(polm(x))

  if ((dv[1]*dv[2]) == 0) stop('replacement has length 0')

  # bring degrees of 'x' and of 'value' in line
  if (dv[3] > dx[3]) {
    x = dbind(d = 3, x, array(0, dim = c(dx[1], dx[2], dv[3] - dx[3])) )
    p = dv[3] - 1
  }
  if (dv[3] < dx[3]) {
    value = dbind(d = 3, value, array(0, dim = c(dv[1], dv[2], dx[3] - dv[3])) )
  }

  # coerce 'x' and 'value' to "vectors"
  dim(x) = c(m*n, p+1)
  dim(value) = c(dv[1]*dv[2], p+1)
  # print(value)

  # extend 'value' if needed
  if ( (length(idx) %% nrow(value)) != 0 ) {
    warning('number of items to replace is not a multiple of replacement length')
  }
  value = value[rep_len(1:nrow(value), length(idx)),,drop = FALSE]
  # print(value)

  # plug in new values
  x[idx,] = value
  # reshape
  dim(x) = c(m, n, p+1)
  # re-coerce to polm
  x = polm(x)

  return(x)
}

'[<-.lpolm' = function(x, i, j, value){

  # Check inputs
  names_args = names(sys.call())
  if (!all(names_args[-length(names_args)] == '')) {
    stop('named dimensions are not supported')
  }

  # Info about inputs
  n_args = nargs() - 2
  is_missing = c(missing(i), missing(j))

  # Extract attributes and transform to array
  attr_x = attributes(x)
  attributes(x) = NULL
  dim_x = attr_x$dim
  dim(x) = dim_x
  min_deg_x = attr_x$min_deg

  # x = unclass(x)
  # dx = dim(x)
  # m = dx[1]
  # n = dx[2]
  # p = dx[3] - 1

  # Extract linear indices
  idx = try(extract_matrix_(dim_x[1], dim_x[2], n_args, is_missing, i, j), silent = TRUE)
  if (inherits(idx, 'try-error')) stop('index/subscripts out of bounds')
  idx = as.vector(idx)

  if (inherits(value, "polm")){
    value = as.lpolm(value)
  }

  if (!inherits(value,'lpolm')) {
    # coerce right hand side 'value' to lpolm object
    value = try(lpolm(value, min_deg = 0) )
    if (inherits(value, 'try-error')) {
      stop('Could not coerce the right hand side to a lpolm object with min_deg = 0!')
    }
  }
  attr_v = attributes(value)
  attributes(value) = NULL
  dim_v = attr_v$dim
  dim(value) = dim_v
  min_deg_v = attr_v$min_deg

  # value = unclass(value)
  # dim_v = dim(value)

  # no items to replace: return original object
  if (length(idx) == 0) return(lpolm(x, min_deg = min_deg_x))

  if ((dim_v[1]*dim_v[2]) == 0) stop('replacement has length 0')

  p1 = dim_x[3]-1+min_deg_x
  p2 = dim_v[3]-1+min_deg_v

  min_q = min(min_deg_x, min_deg_v)
  max_p = max(p1, p2)

  x = dbind(d = 3,
             array(0, dim = c(dim_x[1], dim_x[2], -min_q+min_deg_x)),
             x,
             array(0, dim = c(dim_x[1], dim_x[2], max_p-p1)))
  value = dbind(d = 3,
             array(0, dim = c(dim_v[1], dim_v[2], -min_q+min_deg_v)),
             value,
             array(0, dim = c(dim_v[1], dim_v[2], max_p-p2)))
  dim_x = dim(x)
  dim_v = dim(value)
#
#
#   # bring degrees of 'x' and of 'value' in line
#   if (dim_v[3] > dim_x[3]) {
#     x = dbind(d = 3,
#               x,
#               array(0, dim = c(dim_x[1], dim_x[2], dim_v[3] - dim_x[3])))
#     dim_x[3] = dim_v[3]
#   }
#   if (dim_v[3] < dim_x[3]) {
#     value = dbind(d = 3,
#                   value,
#                   array(0, dim = c(dim_v[1], dim_v[2], dim_x[3] - dim_v[3])))
#     dim_v[3] = dim_x[3]
#   }

  # coerce 'x' and 'value' to "vectors"
  dim(x) = c(dim_x[1]*dim_x[2], max(dim_x[3], dim_v[3]))
  #cat(dim_v)
  #cat(dim(value))
  dim(value) = c(dim_v[1]*dim_v[2], dim_v[3])

  # extend 'value' if needed
  if ( (length(idx) %% nrow(value)) != 0 ) {
    warning('number of items to replace is not a multiple of replacement length')
  }
  value = value[rep_len(1:nrow(value), length(idx)),,drop = FALSE]
  # print(value)

  # plug in new values
  x[idx,] = value
  # reshape
  dim(x) = c(dim_x[1], dim_x[2], dim_x[3])
  # re-coerces to polm
  x = lpolm(x, min_deg = min_q)

  return(x)

}



'[.polm' = function(x,i,j) {
  # print(sys.call())
  if (!is.null(names(sys.call()))) {
    stop('named dimensions are not supported')
  }
  n_args = nargs() - 1
  is_missing = c(missing(i), missing(j))
  # x[] or x[,]
  if ( ((n_args==1) && is_missing[1]) || all(is_missing) ) return(x)

  x = unclass(x)
  d = dim(x)
  m = d[1]
  n = d[2]
  p = d[3] - 1

  idx = try(extract_matrix_(m, n, n_args, is_missing, i, j), silent = TRUE)
  if (inherits(idx, 'try-error')) stop('index/subscripts out of bounds')
  # print(idx)

  if (length(idx) == 0) {
    # result is an "empty" rational matrix
    x = polm(array(0, dim = c(nrow(idx), ncol(idx), 1)))
    return(x)
  }

  if (n_args == 1) {
    dim(x) = c(m*n,1,p+1)
    x = polm(x[i,1,,drop = FALSE])
    return(x)
  }

  if (is_missing[1]) {
    # x[,j]
    x = polm(x[,j,,drop = FALSE])
    return(x)
  }

  if (is_missing[2]) {
    # x[i,]
    x = polm(x[i,,,drop=FALSE])
    return(x)
  }

  # x[i,j]
  x = polm(x[i,j,,drop=FALSE])
  return(x)
}

'[.lpolm' = function(x,i,j){

  # print(sys.call())
  if (!is.null(names(sys.call()))) {
    stop('named dimensions are not supported')
  }
  n_args = nargs() - 1
  is_missing = c(missing(i), missing(j))
  # x[] or x[,]
  if ( ((n_args==1) && is_missing[1]) || all(is_missing) ) return(x)

  attr_x = attributes(x)
  min_deg = attr_x$min_deg
  attributes(x) = NULL
  d = attr_x$dim
  dim(x) = d

  # x = unclass(x)
  # d = dim(x)
  # m = d[1]
  # n = d[2]
  # p = d[3] - 1

  idx = try(extract_matrix_(d[1], d[2], n_args, is_missing, i, j), silent = TRUE)
  if (inherits(idx, 'try-error')) stop('index/subscripts out of bounds')
  # print(idx)

  if (length(idx) == 0) {
    # result is an "empty" rational matrix
    x = polm(array(0, dim = c(nrow(idx), ncol(idx), 1)))
    return(x)
  }

  if (n_args == 1) {
    dim(x) = c(d[1]*d[2],1,d[3])
    x = lpolm(x[i,1,,drop = FALSE], min_deg = attr_x$min_deg)
    return(x)
  }

  if (is_missing[1]) {
    # x[,j]
    x = lpolm(x[,j,,drop = FALSE], min_deg = attr_x$min_deg)
    return(x)
  }

  if (is_missing[2]) {
    # x[i,]
    x = lpolm(x[i,,,drop=FALSE], min_deg = attr_x$min_deg)
    return(x)
  }

  # x[i,j]
  x = lpolm(x[i,j,,drop=FALSE], min_deg = attr_x$min_deg)
  return(x)


}

'[.lmfd' = function(x,i,j) {
  n_args = nargs() - 1
  is_missing = c(missing(i), missing(j))
  if(!all(is_missing == c(TRUE, FALSE)) && n_args != 2){
    stop("Only columns of lmfd()s can be subset.")
  }

  lmfd_obj = x
  x = lmfd_obj$b
  x = unclass(x)
  d = dim(x)
  m = d[1]
  n = d[2]

  idx = try(extract_matrix_(m, n, n_args, is_missing, i, j), silent = TRUE)
  if (inherits(idx, 'try-error')) stop('index/subscripts out of bounds')

  if (length(idx) == 0) {
    # result is an "empty" rational matrix
    x = polm(array(0, dim = c(nrow(idx), ncol(idx), 1)))
    return(lmfd(a = lmfd_obj$a, b = x))
  }

  if (is_missing[1]) {
    # x[,j]
    x = polm(x[,j,,drop = FALSE])
    return(lmfd(a = lmfd_obj$a, b = x))
  }

  stop("This should not happen.")
}

'[.rmfd' = function(x,i,j) {

  n_args = nargs() - 1
  is_missing = c(missing(i), missing(j))
  if(!all(is_missing == c(FALSE, TRUE)) && n_args != 2){
    stop("Only rows of rmfd()s can be subset.")
  }

  rmfd_obj = x
  x = rmfd_obj$d
  x = unclass(x)
  d = dim(x)
  m = d[1]
  n = d[2]

  idx = try(extract_matrix_(m, n, n_args, is_missing, i, j), silent = TRUE)
  if (inherits(idx, 'try-error')) stop('index/subscripts out of bounds')

  if (length(idx) == 0) {
    # result is an "empty" rational matrix
    x = polm(array(0, dim = c(nrow(idx), ncol(idx), 1)))
    return(rmfd(c = rmfd_obj$c, d = x))
  }

  if (is_missing[2]) {
    # x[i,]
    x = polm(x[i,,,drop = FALSE])
    return(rmfd(c = rmfd_obj$c, d = x))
  }

  stop("This should not happen.")
}

'[.stsp' = function(x,i,j) {
  # print(sys.call())
  if (!is.null(names(sys.call()))) {
    stop('named dimensions are not supported')
  }
  n_args = nargs() - 1
  is_missing = c(missing(i), missing(j))
  # print(is_missing)

  d = attr(x, 'order')
  m = d[1]
  n = d[2]
  s = d[3]

  idx = try(extract_matrix_(m, n, n_args, is_missing, i, j), silent = TRUE)
  if (inherits(idx, 'try-error')) stop('index/subscripts out of bounds')
  # print(idx)

  if (length(idx) == 0) {
    # result is an "empty" rational matrix
    x = stsp(A = matrix(0, nrow = 0, ncol = 0), B = matrix(0, nrow = 0, ncol = ncol(idx)),
             C = matrix(0, nrow = nrow(idx), ncol = 0), D = matrix(0, nrow = nrow(idx), ncol = ncol(idx)))
    return(x)
  }

  # x[] or x[,]
  if ( ((n_args==1) && is_missing[1]) || all(is_missing) ) return(x)

  x = unclass(x)
  A = x[iseq(1,s), iseq(1,s), drop = FALSE]
  B = x[iseq(1,s), iseq(s+1,s+n), drop = FALSE]
  C = x[iseq(s+1,s+m), iseq(1,s), drop = FALSE]
  D = x[iseq(s+1,s+m), iseq(s+1,s+n), drop = FALSE]

  if (n_args == 1) {
    transposed = FALSE
    if (m < n) {
      # in order to get a statespace model with a "small" statespace dimension
      # transpose matrix
      transposed = TRUE
      A = t(A)
      D = t(D)
      junk = B
      B = t(C)
      C = t(junk)
      i = matrix(1:(m*n), nrow = m, ncol = n)
      i = t(i)[idx]
      m = nrow(D)
      n = ncol(D)
    }
    A = do.call('bdiag', rep(list(A), n))
    C = do.call('bdiag', rep(list(C), n))
    dim(B) = c(s*n,1)
    dim(D) = c(m*n,1)
    if (transposed) {
      A = t(A)
      D = t(D)
      junk = B
      B = t(C)
      C = t(junk)
    }
    C = C[i,,drop = FALSE]
    D = D[i,,drop = FALSE]
    x = stsp(A = A, B = B, C = C, D = D)
    return(x)
  }

  if (is_missing[1]) {
    # x[,j]
    B = B[,j,drop = FALSE]
    D = D[,j,drop = FALSE]
    x = stsp(A = A, B = B, C = C, D = D)
    return(x)
  }
  if (is_missing[2]) {
    # x[i,]
    C = C[i,,drop = FALSE]
    D = D[i,,drop = FALSE]
    x = stsp(A = A, B = B, C = C, D = D)
    return(x)
  }
  # x[i,j]
  B = B[,j,drop = FALSE]
  C = C[i,,drop = FALSE]
  D = D[i,j,drop = FALSE]
  x = stsp(A = A, B = B, C = C, D = D)
  return(x)
}

'[.pseries' = function(x,i,j) {
  # print(sys.call())
  if (!is.null(names(sys.call()))) {
    stop('named dimensions are not supported')
  }
  n_args = nargs() - 1
  is_missing = c(missing(i), missing(j))
  # x[] or x[,]
  if ( ((n_args==1) && is_missing[1]) || all(is_missing) ) return(x)

  x = unclass(x)
  d = dim(x)
  m = d[1]
  n = d[2]
  lag.max = d[3] - 1

  idx = try(extract_matrix_(m, n, n_args, is_missing, i, j), silent = TRUE)
  if (inherits(idx, 'try-error')) stop('index/subscripts out of bounds')
  # print(idx)

  if (length(idx) == 0) {
    # result is an "empty" impulse response
    x = array(0, dim = c(nrow(idx), ncol(idx), lag.max+1))
    x = structure(x, class = c('pseries', 'ratm'))
    return(x)
  }

  if (n_args == 1) {
    dim(x) = c(m*n, 1, lag.max + 1)
    x = x[i,1,,drop = FALSE]
    x = structure(x, class = c('pseries', 'ratm'))
    return(x)
  }

  if (is_missing[1]) {
    # x[,j]
    x = x[,j,,drop = FALSE]
    x = structure(x, class = c('pseries', 'ratm'))
    return(x)
  }

  if (is_missing[2]) {
    # x[i,]
    x = x[i,,,drop=FALSE]
    x = structure(x, class = c('pseries', 'ratm'))
    return(x)
  }

  # x[i,j]
  x = x[i,j,,drop=FALSE]
  x = structure(x, class = c('pseries', 'ratm'))
  return(x)
}

'[.zvalues' = function(x,i,j) {
  # print(sys.call())
  if (!is.null(names(sys.call()))) {
    stop('named dimensions are not supported')
  }
  n_args = nargs() - 1
  is_missing = c(missing(i), missing(j))
  # x[] or x[,]
  if ( ((n_args==1) && is_missing[1]) || all(is_missing) ) return(x)

  x = unclass(x)
  z = attr(x,'z')
  d = dim(x)
  m = d[1]
  n = d[2]
  n.z = d[3]

  idx = try(extract_matrix_(m, n, n_args, is_missing, i, j), silent = TRUE)
  if (inherits(idx, 'try-error')) stop('index/subscripts out of bounds')
  # print(idx)

  if (length(idx) == 0) {
    # result is an "empty" frequency response
    x = array(0, dim = c(nrow(idx), ncol(idx), n.z))
    x = structure(x, z = z, class = c('zvalues','ratm'))
    return(x)
  }

  if (n_args == 1) {
    dim(x) = c(m*n, 1, n.z)
    x = x[i,1,,drop = FALSE]
    x = structure(x, z = z, class = c('zvalues','ratm'))
    return(x)
  }

  if (is_missing[1]) {
    # x[,j]
    x = x[,j,,drop = FALSE]
    x = structure(x, z = z, class = c('zvalues','ratm'))
    return(x)
  }

  if (is_missing[2]) {
    # x[i,]
    x = x[i,,,drop=FALSE]
    x = structure(x, z = z, class = c('zvalues','ratm'))
    return(x)
  }

  # x[i,j]
  x = x[i,j,,drop=FALSE]
  x = structure(x, z = z, class = c('zvalues','ratm'))
  return(x)
}


'$.lmfd' = function(x, name) {
  i = match(name, c('a','b'))
  if (is.na(i)) stop('reference to "',name, '" is not defined')
  d = attr(x,'order')
  m = d[1]
  n = d[2]
  p = d[3]
  q = d[4]
  x = unclass(x)
  if (i == 1) {
    return(polm(array(x[,iseq(1,m*(p+1))], dim = c(m,m,p+1))))
  }
  if (i == 2) {
    return(polm(array(x[,iseq(m*(p+1)+1,m*(p+1)+n*(q+1))], dim = c(m,n,q+1))))
  }
  # this should not happen
  stop('unknown reference')
}


'$.rmfd' = function(x, name) {
  i = match(name, c('c','d'))
  if (is.na(i)) stop('reference to "',name, '" is not defined')
  d = attr(x,'order')
  m = d[1]
  n = d[2]
  p = d[3]
  q = d[4]
  x = t(unclass(x)) # transposing is necessary because RMFDs are stacked one above the other
  if (i == 1) {
    w = array(x[,iseq(1,n*(p+1))], dim = c(n,n,p+1))
    return(polm(aperm(w, c(2,1,3))))
  }
  if (i == 2) {
    w = array(x[,iseq(n*(p+1)+1,n*(p+1)+m*(q+1))], dim = c(n,m,q+1))
    return(polm(aperm(w, c(2,1,3))))
  }
  # this should not happen
  stop('unknown reference')
}

'$.stsp' = function(x, name) {
  i = match(name, c('A','B','C','D'))
  if (is.na(i)) stop('reference to "',name, '" is not defined')
  d = attr(x,'order')
  m = d[1]
  n = d[2]
  s = d[3]
  x = unclass(x)
  if (i == 1) {
    return(x[iseq(1,s),iseq(1,s),drop = FALSE])
  }
  if (i == 2) {
    return(x[iseq(1,s),iseq(s+1,s+n),drop = FALSE])
  }
  if (i == 3) {
    return(x[iseq(s+1,s+m),iseq(1,s),drop = FALSE])
  }
  if (i == 4) {
    return(x[iseq(s+1,s+m),iseq(s+1,s+n),drop = FALSE])
  }
  # this should not happen
  stop('unknown reference')
}

'$.zvalues' = function(x, name) {
    i = match(name, c('z','f'))
    if (is.na(i)) stop('reference to "',name, '" is not defined')
    z = attr(x, 'z')
    if (i == 1) {
      return(z)
    }
    if (i == 2) {
      return( -Arg(z)/(2*pi) )
    }
    # this should not happen
    stop('unknown reference')
  }
  # is.____ methods ##############################################################

is.polm = function(x) {
  if (!inherits(x,'polm')) return(FALSE)
  if (!inherits(x,'ratm')) return(FALSE)

  x = unclass(x)
  if (! (is.numeric(x) || is.complex(x)) ) {
    return(FALSE)
  }

  if ( (!is.array(x)) || (!(length(dim(x))==3)) ) {
    return(FALSE)
  }

  return(TRUE)
}

is.lpolm = function(x) {
  if (!inherits(x,'lpolm')) return(FALSE)
  if (!inherits(x,'ratm')) return(FALSE)

  x = unclass(x)
  if (! (is.numeric(x) || is.complex(x)) ) {
    return(FALSE)
  }

  if ( (!is.array(x)) || (!(length(dim(x))==3)) ) {
    return(FALSE)
  }

  min_deg = attr(x,'min_deg')
  if ( (is.null(min_deg)) || (length(min_deg) != 1) ) {
    return(FALSE)
  }

  return(TRUE)
}


is.lmfd = function(x) {
  if (!inherits(x,'lmfd')) return(FALSE)
  if (!inherits(x,'ratm')) return(FALSE)

  x = unclass(x)
  if (! (is.numeric(x) || is.complex(x)) ) {
    return(FALSE)
  }

  if (!is.matrix(x)) {
    return(FALSE)
  }

  order = attr(x,'order')
  if ( (is.null(order)) || (!is.integer(order)) || (length(order) != 4) ||
       (nrow(x) == 0) || (order[1] != nrow(x)) || (order[3] < 0) ||
       (order[2] < 0) || (order[4] < -1) ||
       ((order[1]*(order[3]+1) + order[2]*(order[4]+1)) != ncol(x)) ) {
    return(FALSE)
  }

  return(TRUE)
}


is.rmfd = function(x) {
  if (!inherits(x,'rmfd')) return(FALSE)
  if (!inherits(x,'ratm')) return(FALSE)

  x = unclass(x)
  if (! (is.numeric(x) || is.complex(x)) ) {
    return(FALSE)
  }

  if (!is.matrix(x)) {
    return(FALSE)
  }

  order = attr(x,'order')
  if ( (is.null(order)) || (!is.integer(order)) || (length(order) != 4) ||
       (ncol(x) == 0) || (order[2] != ncol(x)) || (order[3] < 0) ||
       (order[1] < 0) || (order[4] < -1) ||
       ((order[1]*(order[4]+1) + order[2]*(order[3]+1)) != nrow(x)) ) {
    return(FALSE)
  }

  return(TRUE)
}


is.stsp = function(x) {
  if (!inherits(x,'stsp')) return(FALSE)
  if (!inherits(x,'ratm')) return(FALSE)

  x = unclass(x)
  if (! (is.numeric(x) || is.complex(x)) ) {
    return(FALSE)
  }

  if (!is.matrix(x)) {
    return(FALSE)
  }

  order = attr(x, 'order')
  if ( (is.null(order)) || (!is.integer(order)) || (length(order) != 3) ||
       (min(order) < 0) || ((order[1]+order[3]) != nrow(x)) ||
       ((order[2]+order[3]) != ncol(x)) ) {
    return(FALSE)
  }

  return(TRUE)
}

is.pseries = function(x) {
  if (!inherits(x,'pseries')) return(FALSE)
  if (!inherits(x,'ratm')) return(FALSE)

  x = unclass(x)
  if ( !( is.numeric(x) || is.complex(x) ) ) {
    return(FALSE)
  }

  if ( (!is.array(x)) || (!(length(dim(x))==3)) ) {
    return(FALSE)
  }

  return(TRUE)
}


is.zvalues = function(x) {
  if (!inherits(x,'zvalues')) return(FALSE)
  if (!inherits(x,'ratm')) return(FALSE)

  x = unclass(x)
  if ( !( is.numeric(x) || is.complex(x) ) ) {
    return(FALSE)
  }

  if ( (!is.array(x)) || (!(length(dim(x))==3)) ) {
    return(FALSE)
  }

  z = attr(x, 'z')
  if ( (is.null(z)) || ( !( is.numeric(z) || is.complex(z) ) ) ||
       (!is.vector(z)) || (length(z) != dim(x)[3]) ) {
    return(FALSE)
  }

  return(TRUE)
}


is.stable = function(x, ...) {
    UseMethod("is.stable", x)
}

is.stable.ratm = function(x, ...) {
  if (class(x)[1] %in% c('polm','lmfd','rmfd','stsp')) {
    z = poles(x)
    return(all(abs(z) > 1))
  }
  return(NA)
}

is.miniphase = function(x, ...) {
  UseMethod("is.miniphase", x)
}

is.miniphase.ratm = function(x, ...) {
  d = dim(x)
  if ((d[1] != d[2]) || (d[1]==0)) stop('the miniphase property is only defined for square non-empty matrices')

  if (class(x)[1] %in% c('polm','lmfd','rmfd','stsp')) {
    z = zeroes(x)
    return(all(abs(z) > 1))
  }
  return(NA)
}

# lyapunov.R
#
# These are essentially R wrapper functions for the Rcpp routines!

lyapunov = function(A, Q,
                    non_stable = c("ignore", "warn", "stop"),
                    attach_lambda = FALSE) {
  # check inputs
  if ( (!is.numeric(A)) || (!is.matrix(A)) || (nrow(A) != ncol(A)) ) {
    stop('"A" must be a square, numeric  matrix ')
  }

  if ( (!is.numeric(Q)) || (!is.matrix(Q)) || any(dim(Q) != dim(A)) ) {
    stop('"Q" must be a numeric matrix with the same dimension as "A"')
  }

  non_stable = match.arg(non_stable)
  attach_lambda = as.logical(attach_lambda)[1]

  m = ncol(A)
  if (m == 0) {
    stop('A,Q are "empty"')
  }

  P = matrix(0, nrow = m, ncol = m)
  lambda_r = numeric(m)
  lambda_i = numeric(m)
  stop_if_non_stable = (non_stable == 'stop')

  # call RcppArmadillo routine
  is_stable = lyapunov_cpp(A, Q, P, lambda_r, lambda_i, stop_if_non_stable)

  if ((stop_if_non_stable) && (!is_stable)) stop('"A" matrix is not stable')
  if ((non_stable == 'warn') && (!is_stable)) warning('"A" matrix is not stable')
  if (attach_lambda) {
      attr(P,'lambda') = complex(real = lambda_r, imaginary = lambda_i)
    }

  return(P)
}

lyapunov_Jacobian = function(A, Q, dA, dQ,
                             non_stable = c("ignore", "warn", "stop")) {
  # check inputs
  if ( (!is.numeric(A)) || (!is.matrix(A)) || (nrow(A) != ncol(A)) ) {
    stop('"A" must be a square, numeric  matrix ')
  }

  if ( (!is.numeric(Q)) || (!is.matrix(Q)) || any(dim(Q) != dim(A)) ) {
    stop('"Q" must be a numeric matrix with the same dimension as "A"')
  }

  m = nrow(A)

  if ( (!is.numeric(dA)) || (!is.matrix(dA)) || (nrow(dA) != (m^2) ) ) {
    stop('"dA" must be a (m^2, k) dimensional numeric  matrix ')
  }

  if ( (!is.numeric(dQ)) || (!is.matrix(dQ)) || any(dim(dA) != dim(dQ)) ) {
    stop('"dQ" is not compatible with "A" or "dA"')
  }

  n = ncol(dA)

  non_stable = match.arg(non_stable)

  if (m*n == 0) {
    stop('A, Q, dA or dQ are "empty"')
  }

  P = matrix(0, nrow = m, ncol = m)
  J = matrix(0, nrow = m^2, ncol = n)
  lambda_r = numeric(m)
  lambda_i = numeric(m)
  stop_if_non_stable = (non_stable == 'stop')

  # call RcppArmadillo routine
  is_stable = lyapunov_Jacobian_cpp(A, Q, P,
                    dA, dQ, J, lambda_r, lambda_i, stop_if_non_stable)

  if ((stop_if_non_stable) && (!is_stable)) stop('"A" matrix is not stable')
  if ((non_stable == 'warn') && (!is_stable)) warning('"A" matrix is not stable')

  return(list(P = P, J = J, lambda = complex(real = lambda_r, imaginary = lambda_i), is_stable = is_stable))
}# munkres.R

# use an extra file for this helper algorithm.

munkres = function(C) {
  # define some helper functions.
  # these helper functions are defined "inside" of munkres()
  # such that they have access to the environment of munkres().

  # this subroutine was just used for debugging
  # # print the current state and eventually what to next.
  # print_state = function(next_step = NULL) {
  #   if (!silent) {
  #     # browser()
  #     rownames(state$C) = state$row_is_covered
  #     colnames(state$C) = state$col_is_covered
  #     state$C[state$stars] = NA_real_
  #     state$C[state$primes] = Inf
  #     print(state$C)
  #     if (!is.null(next_step)) {
  #       cat('next step:', next_step, '\n')
  #     }
  #   }
  # }

  # mark zeroes as "stars"
  # each row, each column may contain at most one starred zeroe.
  # if m (<= n) zeroes are starred, then we have found the optimal assignment
  star_zeroes = function() {
    # browser()
    # C is square or wide: (m-by-n) and m <= n.
    m = nrow(state$C)
    state$stars = matrix(0, nrow = 0, ncol = 2)
    cc = logical(ncol(state$C)) # cc[j] is set to TRUE if the column j contains a star
    for (i in (1:m)) {
      j = which((state$C[i, ] == 0) & (!cc))
      if ( length(j) > 0) {
        cc[j[1]] = TRUE
        state$stars = rbind(state$stars, c(i,j[1]))
      }
    }

    state ->> state
    next_step = 'cover cols with stars'
    # print_state(next_step)

    return(next_step)
  }

  # cover/mark the columns which contain a starred zeroe,
  # if there are m (<= n) starred zeroes, then these
  # starred zeroes describe the optimal assignment
  cover_cols_with_stars = function() {
    # browser()
    m = nrow(state$C)
    state$col_is_covered[] = FALSE
    state$col_is_covered[state$stars[, 2]] = TRUE

    state ->> state
    if (sum(state$col_is_covered) == m) {
      next_step = 'done'
    } else {
      next_step = 'prime zeroes'
    }
    # print_state(next_step)

    return(next_step)
  }

  # find a zeroe within the non-covered entries of C
  find_zeroe = function() {
    state$C[state$row_is_covered, ] = 1 # in order to exclude the covered rows from the search
    state$C[, state$col_is_covered] = 1 # in order to exclude the covered cols from the search
    # browser()
    if (any(state$C == 0)) {
      i = which.max(apply(state$C == 0, MARGIN = 1, FUN = any))
      j = which.max(state$C[i, ] == 0)
      return(c(i,j))
    } else {
      return(NULL)
    }
  }


  # find a zeroe within the non-covered entries of C and mark it as "primed".
  # primed zeroes are "candidates" to replace starred zeroes
  prime_zeroes = function() {
    # iter = 1
    # while ((TRUE) && (iter<=100)) {
    while (TRUE) {
      # the function "returns", if no zeroe can be found, or if there is no
      # star in the respective row. otherwise the row is "covered".
      # hence the while loop will stop after at most m iterations.
      ij = find_zeroe()
      if (is.null(ij)) {
        state ->> state
        next_step = 'augment path'
        # print_state(next_step)

        return(next_step)
      } else {
        # if (!silent) cat('prime zeroe', ij, '\n')
        state$primes = rbind(state$primes, ij)
        i = ij[1]
        k = which(state$stars[, 1] == i) # find column with star in the i-th row
        if ( (length(k) == 0)) {
          state$zigzag_path = matrix(ij, nrow = 1, ncol = 2)
          state ->> state
          next_step = 'make zigzag path'
          # print_state(next_step)

          return(next_step)
        } else {
          j = state$stars[k[1], 2]
          state$row_is_covered[i] = TRUE  # cover the i-th row
          state$col_is_covered[j] = FALSE # uncover the j-th column
          state ->> state
        }
      }
      # iter = iter + 1
    }
    # if (iter >= 100)  stop('did not converge')
  }

  # create a path, which connects alternating primed and starred zeroes.
  #
  # the path starts and ends with a prime zeroe.
  # at the end, the starred zeroes of the path are "unstarred" and the primed zeroes
  # get starred zeroes. So we construct an additional starred zeroe!
  #
  make_zigzag_path = function() {
    pos = 1 # current length of zigzag path
    done = FALSE
    # iter = 1
    # while ((!done) && (iter <= 100)) {
    while (!done) {
      # the while loop stops after at most m iterations, since ...

      # the current zeroe is a primed zero
      # find a starred zero in this column
      k = which(state$stars[, 2] == state$zigzag_path[pos, 2])
      if (length(k) == 0) {
        done = TRUE   # if we cannot find a starred zeroe, the zigzag path is finished
      } else {
        pos = pos + 1 # add this starred zeroe to the path
        state$zigzag_path = rbind(state$zigzag_path, c(state$stars[k,]))
      }

      if (!done) {
        # find a primed zero in this row
        k = which(state$primes[, 1] == state$zigzag_path[pos, 1])
        if (length(k) == 0) {
          stop('this should not happen (1)')
        } else {
          pos = pos + 1 # add this primed zeroe to the zigzag path
          state$zigzag_path = rbind(state$zigzag_path, c(state$primes[k,]))
        }
      }
      # iter = iter + 1
    }
    # if (iter >= 100) stop('did not converge')
    # if (!silent) print(state$zigzag_path)

    # modify stars/primes along the path
    for (i in (1:pos)) {
      if ((i %% 2) == 1) {
        # primed zeroe -> starred zeroe
        state$stars = rbind(state$stars, state$zigzag_path[i,])
      } else {
        # starred zeroes get unstarred
        k = which( (state$stars[, 1] == state$zigzag_path[i, 1]) & (state$stars[, 2] == state$zigzag_path[i, 2]) )
        if (length(k == 1)) {
          state$stars = state$stars[-k, , drop = FALSE]
        } else {
          cat('i=', i, '\n')
          print(state$zigzag_path)
          print(state$stars)
          print(state$primes)
          print(k)
          stop('this should not happen (2)')
        }
      }
    }
    state$row_is_covered[] = FALSE   # uncover all rows
    state$col_is_covered[] = FALSE   # uncover all columns
    state$primes = matrix(0, nrow = 0, ncol = 2) # delete all primed zeroes

    state ->> state
    next_step = 'cover cols with stars'
    # print_state(next_step)

    return(next_step)
  }

  # augment path
  # this step is executed, if the matrix C does not contain "uncovered" zeroes.
  # let m be the minimum of the non-covered elements.
  # substract m from the uncovered elements and add m to the elements which
  # are double covered (column and ro is covered).
  augment_path = function() {
    m = min(state$C[!state$row_is_covered, !state$col_is_covered])
    state$C[!state$row_is_covered, !state$col_is_covered] =
      state$C[!state$row_is_covered, !state$col_is_covered] - m
    state$C[state$row_is_covered, state$col_is_covered] =
      state$C[state$row_is_covered, state$col_is_covered] + m

    state ->> state
    next_step = 'prime zeroes'
    # print_state(next_step)

    return(next_step)
  }

  C0 = C
  m = nrow(C)
  n = ncol(C)

  # convert to "wide" matrix
  if (m > n) {
    transposed = TRUE
    C = t(C)
    m = nrow(C)
    n = ncol(C)
  } else {
    transposed = FALSE
  }

  # adding/substracting a scalar to a row or column of C does not
  # change the optimal assignment(s).
  # (However, the minimal total cost is changed.)
  #
  # substract the row minima (from the respective rows) and the
  # column minima (from the respective columns).
  # this gives a matrix with nonnegative entries and at least
  # one zero in each column and each row.
  C = C - matrix(apply(C, MARGIN = 1, FUN = min), m, n)
  # C = C - matrix(apply(C, MARGIN = 2, FUN = min), m, n, byrow = TRUE)

  # construct a "state" variable, which stores the current state.
  state = list(C = C,
               row_is_covered = logical(m),
               col_is_covered = logical(n),
               stars = matrix(0, nrow = 0, ncol = 2),
               primes = matrix(0, nrow = 0, ncol = 2),
               zigzag_path = matrix(0, nrow = 0, ncol = 2))

  next_step = star_zeroes()

  # iter = 1
  done = FALSE
  while (!done) {
    # print(iter)

    if (next_step == 'cover cols with stars') {
      next_step = cover_cols_with_stars()
      if (next_step == 'done') break   # we are done!!!
    }

    if (next_step == 'prime zeroes') {
      next_step = prime_zeroes()
    }

    if (next_step == 'make zigzag path') {
      next_step = make_zigzag_path()
    }
    if (next_step == 'augment path') {
      next_step = augment_path()
    }
    # iter = iter + 1

  }

  # browser()
  state$C = C0

  # starred zeroes represent the optimal matching
  state$a = state$stars
  if (transposed) state$a = state$a[, c(2,1), drop = FALSE]
  state$a = state$a[order(state$a[, 1]), , drop = FALSE]
  state$c = sum(state$C[state$a])
  return(state[c('a','c','C')])
}
# plot.____ methods ##############################################################


plot.pseries = function(x, x_list = NULL,
                        xlim = c('global','column','subfig'), ylim = c('row','subfig','global'),
                        main = 'impulse response', xlab = 'lag (k)', ylab = NULL,
                        subfigure_main = NA, parse_subfigure_main = FALSE,
                        style = c('gray', 'bw', 'bw2', 'colored'),
                        col = NA, type = 'l', lty = 'solid', lwd = 1,
                        pch = 16, cex.points = 1, bg.points = 'black',
                        legend = NULL, legend_args = NA, ...) {
  style = match.arg(style)
  xlim = match.arg(xlim)
  ylim = match.arg(ylim)

  ir = unclass(x)
  m0 = dim(ir)[1]
  n0 = dim(ir)[2]
  lag.max  = dim(ir)[3] - 1
  if ((n0*m0*(lag.max+1)) == 0) stop('empty impulse response')

  y = list(ir)
  x = list(0:lag.max)

  if ((!is.null(x_list)) && (is.list(x_list))) {
    for (i in (1:length(x_list))) {
      if (!is.pseries(x_list[[i]])) stop('argument "x_list" must be a list of "pseries" objects')
      ir = unclass(x_list[[i]])
      m = dim(ir)[1]
      n = dim(ir)[2]
      lag.max  = dim(ir)[3] - 1
      if ( (m != m0) || (n != n0)) stop('impulse response objects are not compatible')
      if ((lag.max+1) == 0) stop('empty impulse response')

      y = c(y, list(ir))
      x = c(x, list(0:lag.max))
    }
  }
  k = length(x)

  if ((!is.null(subfigure_main)) && (!is.expression(subfigure_main)) && is.na(subfigure_main)) {
    if ((m0*n0) > 1) {
      subfigure_main = '(i_,j_)-th entry'
    } else {
      subfigure_main = NULL
    }
  }

  subfigure = plot_3D(x, y,
                      xlim = xlim, ylim = ylim,
                      main = main, ylab = ylab, xlab = xlab,
                      subfigure_main = subfigure_main, parse_subfigure_main = parse_subfigure_main,
                      style = style,
                      col = col, type = type, lty = lty, lwd = lwd, pch = pch,
                      legend = legend, legend_args = legend_args)
  return(invisible(subfigure))
}




plot.zvalues = function(x, x_list = NULL,  style = c('gray', 'bw', 'colored'),
                        which = c('modulus','phase','nyquist','real'),
                        subfigure_main = NA,
                        xlim = NA, ylim = NA,
                        main = NA, ylab = NA, xlab = NA,
                        legend = NULL, legend_args = NA,
                        col = NA, type = 'l', lty = 'solid', lwd = 1,
                        pch = 16, cex.points = 1, bg.points = 'black', ...) {
  style = match.arg(style)
  which = match.arg(which)

  z = attr(x, 'z')
  n.z = length(z)

  fr = unclass(x)
  attr(fr, 'z') = NULL
  m0 = dim(fr)[1]
  n0 = dim(fr)[2]
  if ((n0*m0*n.z) == 0) stop('empty frequency response')

  x = list(z)
  y = list(fr)

  if ((!is.null(x_list)) && (is.list(x_list))) {
    for (i in (1:length(x_list))) {
      if (!is.zvalues(x_list[[i]])) stop('argument "x_list" must be a list of "zvalues" objects')

      z = attr(x_list[[i]],'z')
      n.z = length(z)

      fr = unclass(x_list[[i]])
      attr(fr, 'z') = NULL
      m = dim(fr)[1]
      n = dim(fr)[2]
      if ( (m != m0) || (n != n0)) stop('zvalues objects are not compatible')
      if (n.z == 0) stop('empty frequency response')

      x = c(x, list(z))
      y = c(y, list(fr))
    }
  }
  k = length(x)

  if (which == 'modulus') {
    for (i in (1:k)) {
      x[[i]] = -Arg(x[[i]])/(2*pi)
      y[[i]] = abs(y[[i]])
      o = order(x[[i]])
      x[[i]] = x[[i]][o]
      y[[i]] = y[[i]][ , , o, drop = FALSE]
    }
    if (is.na(main)) main = 'frequency response'
    if (is.na(xlab)) xlab = 'frequency (f)'
    if (is.na(ylab)) ylab = 'modulus - frequency response'
    if (is.na(xlim)) {
      xlim = 'column'
    } else {
      xlim = match.arg(xlim, c('global','column','subfig'))
    }
    if (is.na(ylim)) {
      ylim = 'row'
    } else {
      ylim = match.arg(ylim, c('global','row','subfig'))
    }
  }
  if (which == 'phase') {
    for (i in (1:k)) {
      x[[i]] = -Arg(x[[i]])/(2*pi)
      y[[i]] = Arg(y[[i]])/(2*pi)
      o = order(x[[i]])
      x[[i]] = x[[i]][o]
      y[[i]] = y[[i]][ , , o, drop = FALSE]
    }
    if (is.na(main)) main = 'frequency response'
    if (is.na(xlab)) xlab = 'frequency (f)'
    if (is.na(ylab)) ylab = 'argument - frequency response'
    if (is.na(xlim)) {
      xlim = 'column'
    } else {
      xlim = match.arg(xlim, c('global','column','subfig'))
    }
    if (is.na(ylim)) {
      ylim = 'row'
    } else {
      ylim = match.arg(ylim, c('global','row','subfig'))
    }
  }
  if (which == 'nyquist') {
    for (i in (1:k)) {
      d = dim(y[[i]])
      # "close" the curves
      x[[i]] = Re(y[[i]])[,,c(1:d[3],1), drop = FALSE]
      y[[i]] = Im(y[[i]])[,,c(1:d[3],1), drop = FALSE]
    }
    if (is.na(main)) main = 'Nyquist plot'
    if (is.na(xlab)) xlab = 'real part - frequency response'
    if (is.na(ylab)) ylab = 'imaginary part - frequency response'
    if (is.na(xlim)) {
      xlim = 'subfig'
    } else {
      xlim = match.arg(xlim, c('global','column','subfig'))
    }
    if (is.na(ylim)) {
      ylim = 'subfig'
    } else {
      ylim = match.arg(ylim, c('global','row','subfig'))
    }
  }
  if (which == 'real') {
    for (i in (1:k)) {
      x[[i]] = Re(x[[i]])
      y[[i]] = Re(y[[i]])
      o = order(x[[i]])
      x[[i]] = x[[i]][o]
      y[[i]] = y[[i]][ , , o, drop = FALSE]
    }
    if (is.na(main)) main = 'plot - function values'
    if (is.na(xlab)) xlab = 'x'
    if (is.na(ylab)) ylab = 'f(x)'
    if (is.na(xlim)) {
      xlim = 'column'
    } else {
      xlim = match.arg(xlim, c('global','column','subfig'))
    }
    if (is.na(ylim)) {
      ylim = 'row'
    } else {
      ylim = match.arg(ylim, c('global','row','subfig'))
    }
  }

  if ((!is.null(subfigure_main)) && (!is.expression(subfigure_main)) && is.na(subfigure_main)) {
    if ((m0*n0) > 1) {
      subfigure_main = '(i_,j_)-th entry'
    } else {
      subfigure_main = NULL
    }
  }
  subfigure = plot_3D(x, y, style = style, subfigure_main = subfigure_main,
                      xlim = xlim, ylim = ylim,
                      main = main, ylab = ylab, xlab = xlab, legend = legend, legend_args = legend_args,
                      col = col, type = type, lty = lty, lwd = lwd, pch = pch)
  return(invisible(subfigure))
}

# plot_tools.R ##################################################

plot_3D = function(x, y,
                   xlim = c('subfig','column','global'), ylim = c('subfig','row','global'), log = '',
                   main = NULL, xlab = NULL, ylab = NULL,
                   subfigure_main = '(i_,j_)-th entry', parse_subfigure_main = FALSE,
                   style = c('gray', 'bw', 'bw2', 'colored'),
                   col = NA, type = 'l', lty = 'solid', lwd = 1,
                   pch = 16, cex.points = 1, bg.points = 'black',
                   legend = NULL, legend_args = NA) {
  style = match.arg(style)

  # __Call to style_parameters() ####
  style = style_parameters(style)

  if (!is.list(x)) x = list(x)
  k_x = length(x)
  if (!is.list(y)) y = list(y)
  k_y = length(y)

  if ((k_x == 0) || (k_y == 0)) stop('no x or no y data supplied')
  if (!( (k_x == 1) || (k_y == 1) || (k_x == k_y) )) stop('the lists "x","y" are not compatible')
  k = max(k_x, k_y)
  if ((k_x == 1) && (k > 1)) x = x[rep(1, k)]
  if ((k_y == 1) && (k > 1)) y = y[rep(1, k)]

  # only minimal input checks!
  chk = function(x) {
    (is.numeric(x) || is.complex(x)) && (length(x) > 0) && ( (is.vector(x)) || (is.array(x) && (length(dim(x)) == 3)) )
  }
  if (!all(sapply(x, FUN = chk))) {
    stop('argument "x" must be a list of (non-empty, numeric or complex) vectors or 3D arrays')
  }
  if (!all(sapply(y, FUN = chk))) {
    stop('argument "y" must be a list of (non-empty, numeric or complex) vectors or 3D arrays')
  }

  # get dimensions (m,n)
  dim2 = function(x) {
    if (is.vector(x)) return(c(1,1))
    return(dim(x)[1:2])
  }
  dim_x = sapply(x, FUN = dim2)
  dim_y = sapply(y, FUN = dim2)
  junk = apply(cbind(dim_x, dim_y), FUN = max, MARGIN = 1)
  m = junk[1]
  n = junk[2]

  # convert vectors to 3-D arrays
  conv3D = function(x, m, n) {
    if (is.vector(x)) {
      x = array(x, dim = c(length(x), m, n))
      x = aperm(x, c(2,3,1))
    }
    d = dim(x)[1:2]
    if (!all(d == c(m,n))) stop('array/vector is not compatible')
    return(x)
  }
  x = try(lapply(x, FUN = conv3D, m, n), silent = TRUE)
  if (inherits(x, 'try-error')) stop('x arrays are not compatible')
  y = try(lapply(y, FUN = conv3D, m, n), silent = TRUE)
  if (inherits(y, 'try-error')) stop('y arrays are not compatible')

  # check dim[3]
  for (i in (1:k)) {
    l_x = dim(x[[min(i, k_x)]])[3]
    l_y = dim(y[[min(i, k_y)]])[3]
    if ( l_x != l_y ) stop('x, y arrays are not compatible')
  }

  # data range
  #range_3D = function(x, MARGIN = integer(0)) {
  X = do.call(dbind, c(list(d=3), x))
  if (grepl('x', log)) X[X<=0] = NA_real_
  Y = do.call(dbind, c(list(d=3), y))
  if (grepl('y', log)) Y[Y<=0] = NA_real_
  if (is.character(xlim)) {
    xlim = match.arg(xlim, c('global','column','subfig'))

    # __Call to range_3D() ####
    x_range = switch(xlim,
                     global = range_3D(X),
                     column = range_3D(X, MARGIN = 2),
                     subfig = range_3D(X, MARGIN = 1:2))
  } else {
    xlim = as.vector(xlim)
    if ( (!is.numeric(xlim)) || (length(xlim) != 2) || any(!is.finite(xlim)) ) {
      stop('illegal "xlim" argument!')
    }
    if ( grepl('x', log) && any(xlim <= 0) ) stop('"xlim" has non negative entries and logarithmix x-axis')
    xlim = sort(xlim)
    # print(Y[1,1,])
    Y[(X < xlim[1]) | (X > xlim[2])] = NA_real_
    # print(Y[1,1,])
    x_range = aperm(array(xlim, c(2,m,n)), c(2,3,1))
    xlim = 'global'
  }
  ylim = match.arg(ylim, c('global','row','subfig'))
  y_range = switch(ylim,
                   global = range_3D(Y),
                   row = range_3D(Y, MARGIN = 1),
                   subfig = range_3D(Y, MARGIN = 1:2))

  # __Call to default_colmap() ####
  # ()||() and ()&&() conditions must have length 1
  if (is.null(col) || all(is.na(col))) col = default_colmap(k)
  col = as.vector(col)
  col = col[(0:(k-1)) %% length(col) + 1]

  type = as.vector(type)
  type = type[(0:(k-1)) %% length(type) + 1]
  lty = as.vector(lty)
  lty = lty[(0:(k-1)) %% length(lty) + 1]
  lwd = as.vector(lwd)
  lwd = lwd[(0:(k-1)) %% length(lwd) + 1]
  pch = as.vector(pch)
  pch = pch[(0:(k-1)) %% length(pch) + 1]
  cex.points = as.vector(cex.points)
  cex.points = cex.points[(0:(k-1)) %% length(cex.points) + 1]
  bg.points = as.vector(bg.points)
  bg.points = bg.points[(0:(k-1)) %% length(bg.points) + 1]

  # deal with subfigure_main
  if (!is.null(subfigure_main)) {

    stopifnot('Input argument *subfigure_main* must be of type *character*' = is.character(subfigure_main))

    if (length(subfigure_main) == 1) {
      subfigure_main = matrix(subfigure_main, nrow = m, ncol = n)

      # replace 'place holder' i_, j_ by i,j
      for (i in (1:m)) {
        for (j in (1:n)) {
          subfigure_main[i,j] = gsub('i_',paste(i), subfigure_main[i,j])
          subfigure_main[i,j] = gsub('j_',paste(j), subfigure_main[i,j])
        }
      }
    }

    stopifnot("Input argument *subfigure_main* must be a matrix (or scalar of type character) with appropriate dimensions." =
                is.matrix(subfigure_main) && all(dim(subfigure_main) == c(m, n)))
  }

  if (!is.null(legend)) {
    if (!is.list(legend_args)) {
      legend_args = list(legend = legend, fill = col, border = col, bty = 'n')
    }
    legend_args$legend = legend
    # ()||() and ()&&() conditions must have length 1
    if ( (!is.null(legend_args$fill)) && all(is.na(legend_args$fill)) ) legend_args$fill = col
    if ( (!is.null(legend_args$border)) && all(is.na(legend_args$border)) ) legend_args$border = col
    if ( (!is.null(legend_args$col)) && all(is.na(legend_args$col)) ) legend_args$col = col
    if ( (!is.null(legend_args$lty)) && all(is.na(legend_args$lty)) ) legend_args$lty = lty
    if ( (!is.null(legend_args$lwd)) && all(is.na(legend_args$lwd)) ) legend_args$lwd = lwd
    if ( (!is.null(legend_args$pch)) && all(is.na(legend_args$pch)) ) legend_args$pch = pch
    if ( (!is.null(legend_args$pt.cex)) && all(is.na(legend_args$pt.cex)) ) legend_args$pt.cex = cex.points
    if ( (!is.null(legend_args$pt.bg)) && all(is.na(legend_args$pt.bg)) ) legend_args$pt.bg = bg.points
  } else {
    legend_args = NULL
  }

  tick_labels = array(FALSE, dim = c(m,n,4))
  if (xlim == 'subfig') {
    tick_labels[,,1] = TRUE
  } else {
    tick_labels[m,,1] = TRUE
  }
  if (ylim == 'subfig') {
    tick_labels[,,2] = TRUE
  } else {
    tick_labels[,1,2] = TRUE
  }
  axes = tick_labels
  axes[,,1] = TRUE
  axes[,,2] = TRUE
  titles = array(NA_character_, dim = c(m,n,4))
  if (!is.null(subfigure_main)) titles[,,3] = subfigure_main

  # cat('plot_3D:\n', xlim, ylim, '\n')
  # print(tick_labels[,,1])
  # print(tick_labels[,,2])

  # compute margins
  margins = list(bottom = 0.2, left = 0.2, top = 0.2, right = 0.2)
  for (i in (1:4)) {
    margins[[i]] = margins[[i]] +
      pmax(1.2 * apply(matrix(tick_labels[,,i], nrow = m, ncol = n), MARGIN = (i-1) %% 2 + 1, FUN = any),
           0.2 * apply(matrix(axes[,,i], nrow = m, ncol = n), MARGIN = (i-1) %% 2 + 1, FUN = any),
           1 * apply(matrix(!is.na(titles[,,i]), nrow = m, ncol = n), MARGIN = (i-1) %% 2 + 1, FUN = any))
  }

  # __Call to set_default_par() ####
  # set default graphic parameters
  opar = set_default_par(m,n)

  # __Call to start_plot() ####
  # start the plot
  start_plot(xlab = xlab, ylab = ylab, main = main, legend_args = legend_args)

  # __Call to get_subfigure_layout() ####
  # compute the corresponding subfigures layout
  subfigure_layout = get_subfigure_layout(margins)

  # print(x_range)
  # print(y_range)
  subfigure_fun = function(i = 1, j = 1) {
    # cat(i,j,'\n')
    # print(str(subfigure_layout))
    opar = graphics::par(omi = subfigure_layout$omi, mar = subfigure_layout$mar[i,j,],
                         fig = subfigure_layout$fig[i,j,], new = TRUE)
    # starting a new plot with plot.window() does not work ????
    graphics::plot(x_range[i,j,], y_range[i,j,], type = 'n', log = log,
                   axes = FALSE, frame.plot = FALSE,
                   xlab = NA, ylab = NA, main = NA, sub = NA)
    # cat('subfig_fun:', par('oma'), '\n', mar, '\n', fig, '\n', xlim, '\n', ylim, '\n')
    # graphics::box(which = 'figure', col = 'red')
    return(invisible(opar))
  }

  for (i in (1:m)) {
    for (j in (1:n)) {
      subfigure_fun(i,j)

      # __Call to plot_axes() ####
      plot_axes(axes = axes[i,j,], tick_labels = tick_labels[i,j,], x_date = 0,
                titles = titles[i,j,], style = style, parse_titles = parse_subfigure_main)

      for (l in (1:k)) graphics::lines(x[[l]][i,j,], y[[l]][i,j,],
                                       type = type[l], col = col[l], lty = lty[l],
                                       lwd = lwd[l], pch = pch[l], cex = cex.points[l],
                                       bg = bg.points[l])
    }
  }

  # reset graphic parameters
  graphics::par(opar)
  return(invisible(subfigure_fun))
}

# Level 1 functions ####

style_parameters = function(style = c('bw', 'bw2', 'gray', 'colored')) {
  style = match.arg(style)

  if (style == 'bw') {
    style = list(
      bg.grid = NA,
      border.grid = NA,
      col.grid = grDevices::gray(0.5),
      lty.grid = 'dotted',
      lwd.grid = 1,
      col.axis = 'black',
      col.box = 'black',
      lwd.box = 1,
      bg.labels = NA,
      border.labels = NA,
      col.labels = 'black'
    )
    return(style)
  }
  if (style == 'bw2') {
    style = list(
      bg.grid = NA,
      border.grid = NA,
      col.grid = grDevices::gray(0.5),
      lty.grid = 'dotted',
      lwd.grid = 1,
      col.axis = 'black',
      col.box = 'black',
      lwd.box = 1,
      bg.labels = grDevices::gray(0.75),
      border.labels = 'black',
      col.labels = 'black'
    )
    return(style)
  }
  if (style == 'gray') {
    style = list(
      bg.grid = grDevices::gray(0.9),
      border.grid = NA,
      col.grid = 'white',
      lty.grid = 'solid',
      lwd.grid = 1,
      col.axis = grDevices::gray(0.5),
      col.box = NA,
      lwd.box = 1,
      bg.labels = grDevices::gray(0.5),
      border.labels = NA,
      col.labels = 'white'
    )
    return(style)
  }
  if (style == 'colored') {
    style = list(
      bg.grid = grDevices::gray(0.9),
      border.grid = NA,
      col.grid = 'white',
      lty.grid = 'solid',
      lwd.grid = 1,
      col.axis = grDevices::gray(0.5),
      col.box = NA,
      lwd.box = 1,
      bg.labels = 'orange',
      border.labels = NA,
      col.labels = 'black'
    )
    return(style)
  }
  stop('this should not happen')
}

range_3D = function(x, MARGIN = integer(0)) {
  m = dim(x)[1]
  n = dim(x)[2]
  if (length(MARGIN) == 0) {
    r = suppressWarnings(range(x, na.rm = TRUE)) # 2
    r = array(r, dim = c(2,m,n))
    r = aperm(r, c(2,3,1))
    return(r)
  }
  if ((length(MARGIN) == 1) && (MARGIN == 1)) {
    r = apply(x, MARGIN = 1, FUN = function(m) {suppressWarnings(range(m, na.rm = TRUE))} ) # 2 x m
    r = array(r, dim = c(2,m,n))
    r = aperm(r, c(2,3,1))
    return(r)
  }
  if ((length(MARGIN) == 1) && (MARGIN == 2)) {
    r = apply(x, MARGIN = 2, FUN = function(m) {suppressWarnings(range(m, na.rm = TRUE))} ) # 2 x n
    r = array(r, dim = c(2,n,m))
    r = aperm(r, c(3,2,1))
    return(r)
  }
  if ((length(MARGIN) == 2) && all(MARGIN == c(1,2))) {
    r = apply(x, MARGIN = c(1,2), FUN = function(m) {suppressWarnings(range(m, na.rm = TRUE))} ) # 2 x m x n
    r = aperm(r, c(2,3,1))
    return(r)
  }
  stop('illegal MARGIN parameter')
}

default_colmap = function(n) {
  return( scales::hue_pal()(n) )
}

set_default_par = function(m=1,n=1){
  sc = sqrt((m^2+n^2)/2)
  sc = max(0.5, 0.5 + 0.5/sqrt(sc))
  opar = graphics::par(oma = c(0, 0, 0, 0),
                       fig = c(0,1,0,1), mar = c(0,0,0,0),
                       tcl = -0.1*sc, mgp = c(1, 0.25, 0)*sc, mex = 1,
                       cex = 1, mex = 1,
                       cex.axis = sc, cex.lab = sc, cex.main = sc,
                       xaxs = 'r', yaxs = 'r', col.axis = 'black')
  return(invisible(opar))
}


start_plot = function(xlab = NULL, ylab = NULL,
                      main = NULL, legend_args= NULL) {

  # set outer margins (place for xlab, ylab, main)
  opar = graphics::par(xaxs = 'i', yaxs = 'i')
  graphics::par(oma = c((!is.null(xlab))*1.5, (!is.null(ylab))*1.5,
                        (!is.null(main))*1.5*1.2, 0),
                fig = c(0,1,0,1), mar = c(0,0,0,0))
  omd = graphics::par('omd')

  # create a new figure which fills the region inside the outer margins
  graphics::plot.new()
  graphics::plot.window(c(0,1), c(0,1))

  if (!is.null(legend_args)) {
    # plot legend
    legend_args$x = 'right'
    legend_args$y = NULL
    out = do.call(graphics::legend, legend_args)
    # set right outer margin ( place for the legend )
    omd[2] = out$rect$left
    graphics::par(omd = omd)
  }

  # reset axis-styles
  graphics::par(opar)

  # plot main, xlab, ylab
  if (!is.null(xlab)) graphics::mtext(xlab, cex = 1,
                                      side = 1, line = 0, outer = TRUE)
  if (!is.null(ylab)) graphics::mtext(ylab, cex = 1,
                                      side = 2, line = 0, outer = TRUE)
  if (!is.null(main)) graphics::mtext(main, cex = 1.2,
                                      side = 3, line = 0.15, outer = TRUE)

  return(invisible(graphics::par('mar')))
}

get_subfigure_layout = function(margins) {
  bottom = rev(margins$bottom)
  top = rev(margins$top)
  left = margins$left
  right = margins$right

  dev_size = grDevices::dev.size('in') # width, height of device in inches

  omi = graphics::par('omi')  # outer margin in 'inches' and lines
  # cat('dev_size:', dev_size, '\n')
  l2in = line2inch()

  # convert margins from lines to inches
  bottom_in = bottom * l2in
  left_in = left * l2in
  top_in = top * l2in
  right_in = right * l2in
  m = length(bottom) # number of rows
  n = length(left)   # number of columns
  plot_width = ( dev_size[1] -
                   (omi[2] + sum(left_in) + sum(right_in) + omi[4]) ) / n
  plot_height = ( dev_size[2] -
                    (omi[1] + sum(bottom_in) + sum(top_in) + omi[3]) ) / m
  # cat('plot_size:', c(plot_width, plot_height), '\n')
  # x coordinate of the lower left corner of the sub figures
  x = cumsum(c(0, left_in + plot_width + right_in))
  # y coordinate of the lower left corner of the sub figures
  y = cumsum(c(0, bottom_in + plot_height + top_in))
  # cat('x:', x, '\n')
  # cat('diff(x):', diff(x),'\n')
  # cat('y:', y, '\n')
  # cat('diff(y):', diff(y),'\n')

  # convert to NDC coordinates
  x = x / (dev_size[1] - omi[2] - omi[4])
  x[x > 1] = 1
  y = y / (dev_size[2] - omi[1] - omi[3])
  y[y > 1] = 1

  mar = array(0, dim = c(m,n,4))
  fig = array(0, dim = c(m,n,4))

  for (i in (1:m)) {
    for (j in (1:n)) {
      mar[i,j, ] = c(bottom[i], left[j], top[i], right[j])
      fig[i,j, ] = c(x[j], x[j+1], y[i], y[i+1])
    }
  }
  mar = mar[m:1,,,drop = FALSE]
  fig = fig[m:1,,,drop=FALSE]
  return(list(omi = omi, mar = mar, fig = fig))
}

plot_axes = function(axes = c(TRUE, TRUE, FALSE, FALSE), tick_labels = axes, x_date = 0,
                     titles = rep(NA_character_, 4), style = style_parameters('bw'),
                     parse_titles = FALSE) {
  junk = graphics::par('cex.axis', 'cex.lab', 'usr', 'xlog', 'ylog')
  cex.axis = junk$cex.axis
  cex.lab = junk$cex.lab
  usr = junk$usr
  xlim = usr[1:2]
  if (junk$xlog) xlim = 10^xlim
  ylim = usr[3:4]
  if (junk$ylog) ylim = 10^ylim

  # print(style$bg.grid)
  # cat(xlim[1], ylim[1], xlim[2], ylim[2], style$bg.grid, style$border.grid, '\n')
  # graphics::rect(xlim[1], ylim[1], xlim[2], ylim[2],
  #                col = 'orange')
  # lines(xlim, ylim, xpd = NA)
  graphics::rect(xlim[1], ylim[1], xlim[2], ylim[2],
                 col = style$bg.grid, border = style$border.grid)

  at_x = NULL
  at_y = NULL
  for (i in (1:4)) {
    if (axes[i]) {
      if (i == 1) {
        a = 0.2
        delta = 0.25*cex.axis - (1-a)*(1-cex.axis)
      } else {
        delta = 0.25*cex.axis
      }
      if ( (i %% 2) == 1 ) {
        if (inherits(x_date, 'POSIXct')) {
          if (junk$xlog) message('logarithmic scale for "POSIXct" class does not work (well)')
          xlim = as.POSIXct(xlim, origin = x_date)
        }
        if (inherits(x_date, 'Date')) {
          if (junk$xlog) message('logarithmic scale for "Date" class does not work (well)')
          xlim = as.Date(xlim, origin = x_date)
        }
        at_x = graphics::Axis(x = xlim,  at = NULL,
                              lwd = 0, lwd.ticks = 1, col.axis = style$col.axis, mgp = c(0, delta, 0),
                              side = i, labels = tick_labels[i])
      } else {
        at_y = graphics::Axis(x = ylim,  at = NULL,
                              lwd = 0, lwd.ticks = 1, col.axis = style$col.axis, mgp = c(0, delta, 0),
                              side = i, labels = tick_labels[i])
      }
    }
  }

  if (!is.null(at_x)) graphics::abline(v = at_x, col = style$col.grid, lty = style$lty.grid, lwd = style$lwd.grid)
  if (!is.null(at_y)) graphics::abline(h = at_y, col = style$col.grid, lty = style$lty.grid, lwd = style$lwd.grid)

  if (!is.na(style$col.box)) graphics::box(which = 'plot', col = style$col.box, lwd = style$lwd.box)


  for (i in (1:4)) {
    if (!is.na(titles[i])) {
      btext(titles[i], side = i, cex = cex.lab, size = cex.lab,
            col = style$col.labels, bg = style$bg.labels, border = style$border.labels,
            parse_text = parse_titles)
    }
  }

  return(invisible(NULL))
}


# Level 2 functions ####
btext = function(text, side = 3, cex = 1, col = 'black',
                 size = cex, bg = 'lightgray', border = 'lightgray',
                 parse_text = FALSE) {
  l2in = line2inch()
  junk = graphics::par('usr','mai','fin','xlog','ylog')
  usr = junk$usr

  plt = junk$fin - c(junk$mai[2]+junk$mai[4],
                     junk$mai[1]+junk$mai[3]) # size of plot region in inches
  in2x = (usr[2]-usr[1])/plt[1] # lines to usr (x) coordinates
  in2y = (usr[4]-usr[3])/plt[2] # lines to usr (y) coordinates


  # cat(junk$fin, plt, usr, in2x, in2y, '\n')

  if ((side %% 2) == 1) {
    srt = 0
    size = size * l2in * in2y
    if (side == 1) {
      box = c(usr[1], usr[2], usr[3] - size, usr[3])
    } else {
      box = c(usr[1], usr[2], usr[4], usr[4] + size)
    }
  } else {
    size = size * l2in * in2x
    if (side == 2) {
      srt = 90
      box = c(usr[1] - size, usr[1], usr[3], usr[4])
    } else {
      srt = -90
      box = c(usr[2], usr[2] + size, usr[3], usr[4])
    }
  }
  center = c(box[1] + box[2], box[3] + box[4] ) /2
  if (junk$xlog) {
    box[1:2] = 10^box[1:2]
    center[1] = 10^center[1]
  }
  if (junk$ylog) {
    box[3:4] = 10^box[3:4]
    center[2] = 10^center[2]
  }

  if (is.expression(text)) parse_text = FALSE
  if (parse_text) {
    junk = try(parse(text = text), silent = TRUE)
    if (!inherits(junk, 'try-error')) text = junk
  }

  graphics::rect(box[1], box[3], box[2], box[4], col = bg, border = border, xpd = TRUE)
  text(center[1], center[2], text, cex = cex, srt = srt,
       adj = c(0.5, 0.5), col = col, xpd = TRUE)

  return(invisible(c(box, center)))
}


line2inch = function() {
  opar = graphics::par(mar = c(1,1,1,1))
  x = mean(graphics::par('mai'))
  graphics::par(opar)
  return(x)
}

# Not sure where this is needed ####

zoom_plot <- function(p, p0 = p, title = NULL, ...) {

  if (!requireNamespace("shiny", quietly = TRUE)) {
    message('"zoom_plot" needs the "shiny" package')
    return(invisible(NULL))
  }

  # create shiny - App
  app <- list(
    ui = shiny::fluidPage(
      shiny::titlePanel(title),
      shiny::fluidRow(style = "background-color: white;",
                      shiny::column(width = 12,
                                    shiny::plotOutput('main',height="500px") # main-plot
                      )
      ),
      shiny::hr(),
      shiny::fluidRow(style = "background-color: white;",
                      # control-plot with a "brush"
                      shiny::column(width = 12,
                                    shiny::plotOutput('control', height="150px",
                                                      brush = shiny::brushOpts(id = "zoom", direction = "x"))
                      )
      )
    ),
    server = function(input, output) {
      output$main <- shiny::renderPlot({
        if (is.null(input$zoom)) {
          p()
        } else {
          # new range of x-values
          xlim = c(input$zoom$xmin, input$zoom$xmax)
          p(xlim)
        }
      })
      output$control <- shiny::renderPlot({
        p0()
      })
      # output$plot_brushinfo <- renderPrint({
      #   cat("Brush (debounced):\n")
      #   str(input$zoom)
      # })
    }
  )
  # run the shiny-app
  shiny::runApp(app, ...)
}



# poles.___ and zeroes.___ method ############################################################
#
#
poles = function(x, tol = sqrt(.Machine$double.eps), print_message = TRUE, ...) {
  UseMethod("poles", x)
}

zeroes = function(x, tol = sqrt(.Machine$double.eps), print_message = TRUE, ...) {
  UseMethod("zeroes", x)
}

zeroes.polm = function(x, tol = sqrt(.Machine$double.eps), print_message = TRUE, ...) {
  x = prune(x, tol = 0) # cut off zero leading coefficients
  d = dim(unclass(x))
  if ((d[1] != d[2]) || (d[1] == 0) || (d[3] <= 0)) {
    stop('argument "x" must represent a non-empty, square and non-singular polynomial matrix')
  }

  # compute companion matrix
  A = try(companion_matrix(x), silent = TRUE)
  if (inherits(A,'try-error')) {
    stop('Could not generate companion matrix. Coefficient pertaining to smallest degree might be singular.')
  }

  # zero degree polynomial
  if (nrow(A) == 0) return(numeric(0))

  z = eigen(A, only.values=TRUE)$values

  if (any(abs(z) <= tol)){
    z <- z[!(abs(z) <= tol)]
    if (print_message){
      message("There are determinantal roots at (or close to) infinity.\nRoots close to infinity got discarded.")
    }
  }

  return(1/z)
}

zeroes.lpolm = function(x, tol = sqrt(.Machine$double.eps), print_message = TRUE, ...) {
  d = dim(x)
  if ((d[1] != d[2]) || (d[1] == 0)) {
    stop('Zeroes.lpolm(): Argument "x" must represent a non-empty, square and non-singular Laurent polynomial matrix.')
  }

  return(zeroes(polm(unclass(x)), tol, print_message))
}

zeroes.lmfd = function(x, tol = sqrt(.Machine$double.eps), print_message = TRUE, ...) {
  d = dim(x)
  if ((d[1] != d[2]) || (d[1] == 0)) {
    stop('argument "x" must represent a non-empty, square and non-singular rational matrix in LMFD form')
  }
  return(zeroes(x$b, tol, print_message))
}

zeroes.rmfd = function(x, tol = sqrt(.Machine$double.eps), print_message = TRUE, ...) {
  d = dim(x)
  if ((d[1] != d[2]) || (d[1] == 0)) {
    stop('argument "x" must represent a non-empty, square and non-singular rational matrix in RMFD form')
  }
  return(zeroes(x$d, tol, print_message))
}

zeroes.stsp = function(x, tol = sqrt(.Machine$double.eps), print_message = TRUE, ...) {
  d = dim(x)
  if ((d[1] != d[2]) || (d[1] ==0)) {
    stop('argument "x" must represent a non-empty, square and non-singular rational matrix in statespace form')
  }
  D = try(solve(x$D), silent = TRUE)
  if (inherits(D, 'try-error')) {
    stop('cannot compute zeroes, D is not invertibe')
  }

  # statespace dimension = 0, "static" matrix
  if ((d[3]) == 0) return(numeric(0))

  A = x$A - x$B %*% D %*% x$C
  z = eigen(A, only.values=TRUE)$values

  if (any(abs(z) <= tol)){
    z <- z[!(abs(z) <= tol)]
    if (print_message){
      message("There are eigenvalues at (or close to) zero.\nThe corresponding zeroes got discarded.")
    }
  }

  return(1/z)
}

poles.polm = function(x, ...) {
  # polynomial matrices have no poles
  return(numeric(0))
}

poles.lmfd = function(x, tol = sqrt(.Machine$double.eps), print_message = TRUE, ...) {
  d = dim(x)
  if ((d[1] == 0) || (d[3] < 0)) {
    stop('argument "x" does not represent a valid LMFD form')
  }
  return(zeroes(x$a, tol, print_message))
}

poles.rmfd = function(x, tol = sqrt(.Machine$double.eps), print_message = TRUE, ...) {
  d = dim(x)
  if ((d[2] == 0) || (d[3] < 0)) {
    stop('argument "x" does not represent a valid RMFD form')
  }
  return(zeroes(x$c, tol, print_message))
}

poles.stsp = function(x, tol = sqrt(.Machine$double.eps), print_message = TRUE, ...) {
  d = dim(x)

  # empty matrix, or "static" matrix (statespace dimension = 0)
  if (min(d) == 0) return(numeric(0))

  z = eigen(x$A, only.values=TRUE)$values

  if (any(abs(z) <= tol)){
    z <- z[!(abs(z) <= tol)]
    if (print_message){
      message("There are eigenvalues at (or close to) zero.\nThe corresponding poles got discarded.")
    }
  }

  return(1/z)
}# Essentials (computing zeros!) ####

is.coprime = function(a, b = NULL, tol = sqrt(.Machine$double.eps), only.answer = TRUE, debug = FALSE) {

  # check inputs polm objects a(z), b(z)
  if (!inherits(a, 'ratm')) {
    a = try(polm(a), silent = TRUE)
    if (inherits(a, 'try-error')) {
      stop('could not coerce input "a" to a "polm" object')
    }
  }
  if ( !( inherits(a,'polm') || inherits(a,'lmfd') || inherits(a, 'rmfd') ) ) {
    stop('input "a" must be a "polm", "lmfd", or "rmfd" object or an object which is coercible to a "polm" object')
  }

  # Concatenate/transpose MFDs
  if (inherits(a,'lmfd')) {
    c = cbind(a$a, a$b)
  } else if (inherits(a,'rmfd')){
    c = cbind(t(a$c), t(a$d))
  } else {
    if (is.null(b)) {
      c = a
    } else {
      c = try(cbind(a,b), silent = TRUE)
      if (inherits(c, 'try-error')) {
        stop('inputs "a" and "b" must be compatible "polm" objects')
      }
    }
  }

  # skip zero leading coefficients and "unclass"
  c = unclass(prune(c, tol = 0))
  d = dim(c)
  m = d[1]
  n = d[2]
  p = d[3] - 1

  if (m*n == 0) {
    # c is an empty polynomial
    stop('the polynomial "c=[a,b]" is empty')
  }

  if (p == (-1)) {
    # c is a zero polynomial
    if (only.answer) return(FALSE)
    return(list(answer = FALSE, A = NULL, B = NULL, zeroes = NA_real_, m = integer(0), n = integer(0)))
  }

  if (p == 0) {
    # c is a constant polynomial
    dim(c) = c(m,n)
    if (m > n) {
      answer = FALSE
    } else {
      # check singular values
      d = svd(c, nu = 0, nv = 0)$d
      answer = (d[m] >= tol)
    }
    if (only.answer) return(answer)
    if (answer) {
      zeroes = numeric(0)
    } else {
      zeroes = NA_real_
    }
    return(list(answer = answer, A = NULL, B = NULL, zeroes = zeroes, m = integer(0), n = integer(0)))
  }

  # construct Pencil A - Bz corresponding to the polynomial c(z)
  c0 = matrix(c[,,1], nrow = m, ncol = n)
  c = -c[,,-1,drop = FALSE]
  dim(c) = c(m, n*p)

  A = bdiag(c0, diag(n*(p-1)))
  B = rbind(c, diag(1, nrow = n*(p-1), ncol = n*p))

  if (m > n) {
    # c is a "tall" polynomial
    if (only.answer) return(FALSE)
    return(list(answer = FALSE, A = A, B = B, zeroes = NA_real_, m = m, n = n))
  }

  if (debug) {
    message('is.coprime() start:')
    cat('pencil (A - Bz):\n A=\n')
    print(round(A, 4))
    cat('B=\n')
    print(round(B, 4))
  }

  # dimensions of pencil (A - Bz)
  m = nrow(A)
  n = ncol(A)

  row = 1
  col = 1
  step = 1
  mm = integer(0)
  nn = integer(0)
  while ((row <= m) && (col<= n)) {
    # consider the lower, right block A22 = A[row:n, col:m], B22 = B[row:m, col:n]

    # Column Trafo: Separate columns of A22 and B22 ####
    # n1 = dimension of the right kernel of B22
    # make the first n1 columns of B22 equal to zero
    svd_x = svd(B[row:m, col:n, drop = FALSE], nv = (n-col+1), nu = 0)
    n1 = (n-col+1) - sum(svd_x$d > tol)
    svd_x$v = svd_x$v[,(n - col + 1):1, drop = FALSE]
    A[, col:n] = A[, col:n, drop = FALSE] %*% svd_x$v
    B[, col:n] = B[, col:n, drop = FALSE] %*% svd_x$v
    # impose zeroes
    if (n1 > 0) {
      B[row:m, col:(col+n1-1)] = 0
    }

    # Return early: If right-kernel of B22 is empty, i.e. n1 = 0  ####
    # (A22 - B22z) is a square or tall pencil (where B22 has full column rank) => pencil has zeroes
    if (n1 == 0) {
      if (only.answer) {
        return(FALSE)
      }
      mm = c(mm, m-row+1)
      nn = c(nn, n-col+1)

      # ( A22 - B22*z ) is a tall pencil => non trivial left kernel for all z!
      if ((m - row) > (n - col)) {
        return(list(answer = FALSE, A = A, B = B, zeroes = NA_real_, m = mm, n = nn))
      }

      # ( A22 - B22*z ) is a square, regular pencil
      zeroes = eigen(solve(B[row:m, col:n, drop = FALSE], A[row:m, col:n, drop = FALSE]), only.values = TRUE)$values
      return(list(answer = FALSE, A = A, B = B, zeroes = zeroes, m = mm, n = nn))
    }

    # Row Trafo: Separate rows of A and B ####
    # m1 = rank of A22
    # make the last ((m - row + 1) - m1) rows of A22 equal to zero
    # => the first m1 rows are linearly independent
    # Note: m1 may be zero
    svd_x = svd(A[row:m, col:(col+n1-1), drop = FALSE], nu = (m - row + 1), nv = 0)
    m1 = sum(svd_x$d > tol)
    A[row:m, col:n] = t(svd_x$u) %*% A[row:m, col:n, drop = FALSE]
    B[row:m, col:n] = t(svd_x$u) %*% B[row:m, col:n, drop = FALSE]
    # impose zeroes
    if (m1 < (m - row +1)) {
      A[(row+m1):m, col:(col+n1-1)] = 0
    }

    if (debug) {
      message(paste('is.coprime() step ', step, ' (',
                    row, ':', row + m1 - 1, ' x ',
                    col, ':', col + n1 - 1, '):', sep=''))
      cat('pencil (A - Bz):\n A=\n')
      print(round(A, 4))
      cat('B=\n')
      print(round(B, 4))
    }

    row = row + m1
    col = col + n1
    mm = c(mm, m1)
    nn = c(nn, n1)

    step = step + 1
  }

  if (row > m) {
    # coprime!
    if (only.answer) {
      return(TRUE)
    }
    nn[length(nn)] = nn[length(nn)] + (n - sum(nn))
    return(list(answer = TRUE, A = A, B = B, zeroes = numeric(0), m = mm, n = nn))
  }
  # not coprime
  if (only.answer) {
    return(FALSE)
  }
  mm[length(mm)] = mm[length(mm)] + (m - sum(mm))
  return(list(answer = FALSE, A = A, B = B, zeroes = NA_real_, m = mm, n = nn))
}

companion_matrix = function(a) {
  if (!inherits(a, 'polm')) {
    a = try(polm(a), silent = TRUE)
    if (inherits(a, 'try-error')) stop('argument "a" is not coercible to polm object!')
  }
  a = unclass(a)
  d = dim(a)
  if ((d[1] != d[2]) || (d[3] <= 0)) stop('argument "a" must represent a square, non-singular polynomial matrix')

  m = d[1]
  p = d[3] - 1

  if (m > 0) {
    # check a(0)
    a0 = try(solve(matrix(a[,,1], nrow = m, ncol = m)), silent = TRUE)
    if (inherits(a0, 'try-error')) stop('constant term a[0] is non invertible')
  }

  if ((m*p) == 0) return(matrix(0, nrow = 0, ncol = 0))

  # coerce to (m,m(p+1)) matrix
  dim(a) = c(m, m*(p+1))

  # normalize constant term a[0] -> I, a[i] -> - a[0]^{-1} a[i]
  a = (-a0) %*% a[, (m+1):(m*(p+1)), drop = FALSE]

  if (p == 1) {
    return(a)
  }
  return( rbind(a, diag(x = 1, nrow = m*(p-1), ncol = m*p)) )
}


# Degree related ####

degree = function(x, which = c('elements', 'rows', 'columns', 'matrix')) {
  if (!inherits(x, 'polm')) {
    stop('argument "x" must be a "polm" object!')
  }
  which = match.arg(which)
  x = unclass(x)
  # degree of a univariate polynomial (= vector):
  # length of vector - 1 - number of zero leading coefficients
  deg_scalar = function(x) {
    length(x) - sum(cumprod(rev(x) == 0)) - 1
  }
  deg = apply(x, MARGIN = c(1,2), FUN = deg_scalar)
  if (which == 'matrix') return(max(deg))
  if (which == 'columns') return(apply(deg, MARGIN = 2, FUN = max))
  if (which == 'rows') return(apply(deg, MARGIN = 1, FUN = max))
  return(deg)
}

col_end_matrix = function(x) {
  if (!inherits(x, 'polm')) {
    stop('argument "x" must be a "polm" object!')
  }
  d = dim(x)
  x = unclass(x)
  NAvalue = ifelse(is.complex(x), NA_complex_, NA_real_)
  m = matrix(NAvalue, nrow = d[1], d[2])
  if (length(x) == 0) {
    return(m)
  }

  # degree of a univariate polynomial (= vector):
  # length of vector - 1 - number of zero leading coefficients
  deg_scalar = function(x) {
    length(x) - sum(cumprod(rev(x) == 0)) - 1
  }
  deg = apply(x, MARGIN = c(1,2), FUN = deg_scalar)
  col_deg = apply(deg, MARGIN = 2, FUN = max)

  for (i in iseq(1, dim(x)[2])) {
    if (col_deg[i] >= 0) m[,i] = x[,i,col_deg[i]+1]
  }
  return(m)
}

prune = function(x, tol = sqrt(.Machine$double.eps), brutal = FALSE) {

  x_class = class(x)
  if (!inherits(x, c('polm', 'lpolm'))) {
    stop('argument "x" must be a polm or lpolm object!')
  }
  if ("lpolm" %in% x_class){
    min_deg = attr(x, "min_deg")
  }
  x = unclass(x)
  d = dim(x)
  if (min(d) <= 0) {
    # empty polynomial, or polynomial of degree (-1)
    return(polm(array(0, dim = c(d[1], d[2], 0))))
  }

  # step one: Set all small leading coefficients to zero
  issmall = ( (abs(Re(x)) <= tol) & (abs(Im(x)) <= tol) )

  issmall_polm = apply(issmall, MARGIN = c(1,2), FUN = function(x) { rev(cumprod(rev(x))) })

  # apply: returns an array of dim (d[3], d[1], d[2]) if d[3] > 1

  # make sure that issmall_polm is an array (also in the case where the matrix polynomial is constant)
  dim(issmall_polm) = d[c(3,1,2)]

  # permute the dimensions back to the polm form:
  # necessary because apply returns an array of dim (d[3], d[1], d[2]) if d[3] > 1
  issmall_polm = aperm(issmall_polm, c(2,3,1))
  issmall_polm = (issmall_polm == 1)

  # finish step one
  x[issmall_polm] = 0

  # Same steps in the other direction for lpolm
  if ("lpolm" %in% x_class){
    issmall_lpolm = apply(issmall, MARGIN = c(1,2), FUN = function(x) { cumprod(x) })
    dim(issmall_lpolm) = d[c(3,1,2)]
    issmall_lpolm = aperm(issmall_lpolm, c(2,3,1))
    issmall_lpolm = (issmall_lpolm == 1)
    x[issmall_lpolm] = 0
  }

  # step two: drop leading zero matrix coefficients
  keep = apply(!issmall_polm, MARGIN = 3, FUN = any)

  if ("polm" %in% x_class){
    # keep[1] = TRUE # keep the constant
    keep
    x = x[,, keep, drop = FALSE]
  } else if("lpolm" %in% x_class){
    keep_lpolm = apply(!issmall_lpolm, MARGIN = 3, FUN = any)
    keep_lpolm = keep_lpolm * keep
    keep_lpolm = (keep_lpolm == 1)
    x = x[,, keep_lpolm, drop = FALSE]
    # Adjust min_deg
    min_deg = min_deg + sum(cumprod(!keep_lpolm))
  }

  # step three: drop imaginary part if all imaginary parts are small
  if (is.complex(x)) {
    if ( all(abs(Im(x)) <= tol) ) {
      x = Re(x)
    }
  }

  # This option is provided to see, e.g., the lower triangularity of the zero power coefficient
  # matrix when using "transform_lower_triangular"
  if (brutal){
    issmall_brutal = ( (abs(Re(x)) <= tol) & (abs(Im(x)) <= tol) )
    x[issmall_brutal] = 0
  }

  if ("polm" %in% x_class){
    x = polm(x)
  } else if ("lpolm" %in% x_class){
    x = lpolm(x, min_deg = min_deg)
  }

  return(x)
}

# Small internal helpers for column reduction ####

# internal function
# l2 norm of a vector
l2_norm = function(x){
  return(sqrt(sum(x^2)))
}

# internal function
# consider a vector x = c(x[1], ..., x[k],  x[k+1], ..., x[n])
# return a logical  i = c(FALSE, .., FALSE, TRUE, .., TRuE)
# where k is the minimum integer such that | x[s] | <= tol for all s > k
is_small = function(x, tol = sqrt(.Machine$double.eps), count = TRUE) {
  if (length(x) == 0) {
    i = logical(0)
  } else {
    i = rev( cumprod( rev(abs(x) <= tol) ) == 1 )
  }
  if (count) {
    return(sum(i))
  } else {
    return(i)
  }
}

# internal function
# consider a vector x = c(x[1], ..., x[k],  x[k+1], ..., x[n])
# return a logical  i = c(TRUE, .., TRUE, FALSE, .., FALSE)
# where k is the minimum integer such that | x[s] | <= tol for all s > k
is_large = function(x, tol = sqrt(.Machine$double.eps), count = TRUE) {
  if (length(x) == 0) {
    i = logical(0)
  } else {
    i = rev( cumprod( rev(abs(x) <= tol) ) == 0 )
  }
  if (count) {
    return(sum(i))
  } else {
    return(i)
  }
}

# Normal forms (Hermite, Smith, Wiener-Hopf) and essential helpers ####

purge_rc = function(a, pivot = c(1,1), direction = c('down','up','left','right'),
                    permute = TRUE, tol = sqrt(.Machine$double.eps),
                    monic = FALSE, debug = FALSE) {
  direction = match.arg(direction)

  # Check argument 'a'
  stopifnot("purge_rc(): Argument *a* is not a polm object!" = inherits(a, 'polm'))

  # Dimensions of input
  d = dim(unclass(a))
  m = d[1]
  n = d[2]
  p0 = d[3] - 1
  stopifnot("purge_rc(): Argument *a* must have more than zero inputs and outputs!" = m*n != 0)

  # check pivot
  pivot = as.integer(as.vector(pivot))
  stopifnot("purge_rc(): Argument *pivot* must be an integer vector of length 2, 1 <= pivot <= dim(a)" = (length(pivot) == 2) && (min(pivot) > 0) && (pivot[1] <= m) && (pivot[2] <= n))

  i = pivot[1]
  j = pivot[2]

  # If direction is not "down", transform (transpose etc) the polm object
  if (direction == 'up') {
    a = a[m:1, ]
    i = m - (i - 1)
  }
  if (direction == 'right') {
    a = t(a)

    junk = i
    i = j
    j = junk

    junk = m
    m = n
    n = junk
  }
  if (direction == 'left') {
    a = t(a)

    junk = i
    i = j
    j = junk

    junk = m
    m = n
    n = junk

    a = a[m:1, ]
    i = m - (i - 1)
  }

  # Initialization of unimodular matrix
  u0 = polm(diag(m))
  u = u0
  u_inv = u0

  # (m x n) matrix of degrees of each element
  p = degree(a)

  # degrees of entries in the j-th column
  p_col = p[, j]

  # no permutations allowed, but pivot element is zero!
  if ( (i < m) && (!permute) && (p_col[i] == -1) && any(p_col[(i+1):m] > -1) ) {
    stop("purge_rc(): Pivot element is zero but permutation is not allowed. Purging not possible.")
  }

  # Main iteration ####
  iteration = 0

  # The column is not purged if
  # (in the case where permutations are allowed) any element below the pivot is non-zero
  # (in the case where permutations are not allowed) any element below the is non-zero is of equal or larger degree than the pivot
  not_purged = (i < m) && ( ( permute && any(p_col[(i+1):m] > -1) ) || ( (!permute) && any(p_col[(i+1):m] >= p_col[i]) ) )

  while ( not_purged )  {

    iteration = iteration + 1

    if (debug) {
      message('purge_rc: iteration=', iteration)
      print(a, format = 'i|zj', digits = 2)
      # print(a)
      print(p)
    }

    if (permute){
      # Permutation step

      # find (non zero) entry with smallest degree
      p_col[iseq(1, i-1)] = Inf   # ignore entries above the i-th row
      p_col[p_col == -1]  = Inf   # ignore zero entries
      k = which.min(p_col)

      # permute i-th row and k-th row
      perm = 1:m
      perm[c(k,i)] = c(i,k)

      a = a[perm, ]
      u = u[, perm]
      u_inv = u_inv[perm, ]
      p_col = p_col[perm]
    }

    # Division step

    q = a[(i+1):m, j] %/% a[i, j] # polynomial divison

    M = u0
    Mi = u0

    M[(i+1):m, i] = -q
    Mi[(i+1):m, i] = q

    a = prune(M %r% a, tol = tol)
    u = u %r% Mi
    u_inv = M %r% u_inv

    p = degree(a)

    # degrees of entries in the j-th column
    p_col = p[, j]
    stopifnot("purge_rc(): Reduction of degree failed! This should not happen." = all(p_col[(i+1):m] < p_col[i]))

    # The column is not purged if
    # (in the case where permutations are allowed) any element below the pivot is non-zero
    # (in the case where permutations are not allowed) any element below the is non-zero is of equal or larger degree than the pivot
    not_purged = (i < m) && ( ( permute && any(p_col[(i+1):m] > -1) ) || ( (!permute) && any(p_col[(i+1):m] >= p_col[i]) ) )
  }

  if (monic) {
    if (p[i,j] >= 0) {
      c = unclass(a)[i,j,p[i,j]+1]
      a[i, ] = a[i, ] %/% c # polynomial division, see ?Ops.ratm
      u[, i] = u[, i] * c
      u_inv[i, ] = u_inv[i, ] %/% c
    }
  }

  # If direction is not "down", transform (transpose etc) the polm object back
  if (direction == 'up') {
    a = a[m:1, ]
    u = u[m:1, m:1]
    u_inv = u_inv[m:1, m:1]
  }
  if (direction == 'right') {
    a = t(a)
    u = t(u)
    u_inv = t(u_inv)
  }
  if (direction == 'left') {
    a = t(a[m:1,])
    u = t(u[m:1,m:1])
    u_inv = t(u_inv[m:1,m:1])
  }

  return(list(h = a, u = u, u_inv = u_inv))
}

col_reduce = function(a, tol = sqrt(.Machine$double.eps), debug = FALSE) {

  # Check inputs
  stopifnot("col_reduce(): Input *a* must be a polm object!" = inherits(a, 'polm'))

  # Integer-valued parameters
  x = unclass(a)
  d = dim(x)
  m = d[1]
  n = d[2]
  p = d[3] - 1
  stopifnot("col_reduce(): Input *a* must be a square, non empty, non zero polynomial matrix." = (m*n != 0) && (m == n) && (p >= 0))

  # Initialize unimodular matrices
  v0 = polm(diag(n))
  v = v0
  v_inv = v0

  {# # balance the 'column norms'
  # col_norm = apply(x, MARGIN = 2, FUN = l2_norm)
  # a = a %r% diag(sqrt(mean(col_norm^2))/col_norm)
  # v_inv = v_inv %r% diag(sqrt(mean(col_norm^2))/col_norm)
  # v = diag(col_norm / sqrt(mean(col_norm^2))) %r% v
  }
  # Column degrees (taking rounding error into account)
  # output is a matrix! rows correspond to the norm of the respective column, columns correspond to degrees!
  col_norms = apply(x, MARGIN = c(2,3), FUN = l2_norm)
  col_degrees = apply(col_norms, MARGIN = 1, FUN = is_large, tol = tol, count = TRUE) - 1
  stopifnot("col_reduce(): The input *a* has (close to) zero columns." = min(col_degrees) >= 0)

  # set "small columns" to zero and retrieve column end matrix
  col_end_matrix = matrix(0, nrow = m, ncol = n)
  for (i in (1:n)) {
    x[, i, iseq(col_degrees[i]+2, p+1)] = 0
    col_end_matrix[,i] = x[, i, col_degrees[i]+1]
  }

  # reduce order of polynomial
  p = max(col_degrees)
  x = x[, , 1:(p+1), drop = FALSE]
  a = polm(x)

  # sort by column degrees
  o = order(col_degrees, decreasing = FALSE)
  col_degrees = col_degrees[o]
  col_end_matrix = col_end_matrix[,o,drop = FALSE]
  x = x[,o,,drop = FALSE]
  a = a[,o]
  v = v[o, ]
  v_inv = v_inv[, o]

  # SVD of column end matrix
  svd_x = svd(col_end_matrix, nv = n, nu = 0)

  if (debug) {
    message('col_reduce: initial matrix')
    cat('column degrees:', col_degrees,'\n')
    print(col_end_matrix)
    cat('singular values of column end matrix:', svd_x$d,'\n')
    print(svd_x$v)
  }

  z = polm(c(0,1))
  iter = 0
  while (min(svd_x$d) < tol) {
    iter = iter + 1

    # Skip small entries at the end of svd_x$v[,n] (last singular value)
    k = is_large(svd_x$v[,n])

    if (debug) {
      message('col_reduce: reduce degree of column ',k)
    }

    v_step = v0
    v_step_inv = v0

    b = numeric(n)
    b[1:k] = svd_x$v[1:k,n]/svd_x$v[k,n]
    for (i in (1:k)) {
      v_step_inv[i, k] = b[i] * z^(col_degrees[k] - col_degrees[i])
      v_step[i, k] = -v_step_inv[i, k]
    }
    v_step[k,k] = 1

    a[, k] = a %r% v_step_inv[, k]
    v_inv[, k] = v_inv %r% v_step_inv[, k]
    v = v_step %r% v

    {# balance the 'column norms'
    #     col_norm = apply(unclass(a), MARGIN = 2, FUN = l2_norm)
    # print(col_norm)
    #     a = a %r% diag(sqrt(mean(col_norm^2))/col_norm, nrow = n, ncol = n)
    #     v_inv = v_inv %r% diag(sqrt(mean(col_norm^2))/col_norm, nrow = n, ncol = n)
    #     v = diag(col_norm / sqrt(mean(col_norm^2)), nrow = n, ncol = n) %r% v
    }
    x = unclass(a)
    # eventually the degree of a has been reduced !?
    p = dim(x)[3] - 1

    # recompute col_degrees and col_end_matrix
    x[ , k, iseq(col_degrees[k]+1, p+1)] = 0 # this column has been purged!
    a = polm(x)

    col_norm = apply(matrix(x[, k, ], nrow = m, ncol = p+1), MARGIN = 2, FUN = l2_norm)
    col_degrees[k] = is_large(col_norm, tol = tol, count = TRUE) - 1
    if (col_degrees[k] < 0) {
      print(col_degrees)
      print(col_norm)
      print(x)
      stop('input "a" is (close to) singular')
    }
    x[, k, iseq(col_degrees[k]+2,p+1)] = 0
    col_end_matrix[,k] = x[, k, col_degrees[k] + 1]
    if (max(col_degrees) < p) {
      p = max(col_degrees)
      x = x[,,1:(p+1), drop = FALSE]
    }

    # (re) sort by column degrees
    o = order(col_degrees, decreasing = FALSE)
    col_degrees = col_degrees[o]
    col_end_matrix = col_end_matrix[,o,drop = FALSE]
    x = x[,o,,drop = FALSE]
    a = polm(x)
    v = v[o, ]
    v_inv = v_inv[, o]

    # iterate
    # SVD of column end matrix
    svd_x = svd(col_end_matrix, nv = n, nu = 0)

    if (debug) {
      message('col_reduce: step=', iter)
      cat('column degrees:', col_degrees,'\n')
      print(col_end_matrix)
      cat('singular values of column end matrix:', svd_x$d,'\n')
      print(svd_x$v)
    }
  }

  # Resort column degrees in non-increasing direction
  o = order(col_degrees, decreasing = TRUE)
  col_degrees = col_degrees[o]
  col_end_matrix = col_end_matrix[,o,drop = FALSE]
  x = x[,o,,drop = FALSE]
  a = polm(x)
  v = v[o, ]
  v_inv = v_inv[, o]

  return(list(a = a,
              v = v, v_inv = v_inv,
              col_degrees = col_degrees, col_end_matrix = col_end_matrix))
}

hnf = function(a, from_left = TRUE, tol = sqrt(.Machine$double.eps), debug = FALSE) {

  # Check inputs
  if (!inherits(a, 'polm')) {
    stop('input "a" is not a "polm" object!')
  }

  # skip zero leading coefficients
  a = prune(a, tol = tol)

  # Dimensions
  d = unname(dim(a))
  m = d[1]
  n = d[2]
  if (m*n == 0) {
    stop('input "a" is an empty polynomial matrix!')
  }


  # from_right
  if (!from_left) {
    a = t(a)

    # recompute dimensions
    d = unname(dim(a))
    m = d[1]
    n = d[2]
  }

  # Init of unimodular matrices
  u0 = polm(diag(m))
  u = u0
  u_inv = u0

  i = 0
  j = 0
  pivots = integer(0)
  while ((i<m) && (j<n))
  {
    if (debug) {
      message('hnf pivot: i=', i,' j=',j,'\n')
    }

    p = degree(a)
    # code zero elelements with Inf
    p[p == -1] = Inf

    if (all(is.infinite(p[(i+1):m,j+1]))) {
      # all remaining elements in (j+1)-th column are zero
      j = j+1
      next
    }

    # purge (j+1)-th column
    out = purge_rc(a, pivot = c(i+1, j+1), direction = "down", permute = TRUE,
                   tol = tol, monic = TRUE, debug = debug)
    a = out$h
    u = u %r% out$u
    u_inv = out$u_inv %r% u_inv

    # make sure that elements ABOVE the diagonal are smaller in degree than the diagonal element
    if (i > 0) {
      out = purge_rc(a, pivot = c(i+1,j+1), direction = "up", permute = FALSE,
                     tol = tol, monic = TRUE, debug = debug)
      a = out$h
      u = u %r% out$u
      u_inv = out$u_inv %r% u_inv
    }

    pivots = c(pivots, j+1)
    i = i+1
    j = j+1

  }

  if (from_left){
    return(list(h = a, u = u, u_inv = u_inv, pivots = pivots, rank = length(pivots)))
  } else {
    return(list(h = t(a), u = t(u), u_inv = t(u_inv), pivots = pivots, rank = length(pivots)))
  }
}

snf = function(a, tol = sqrt(.Machine$double.eps), debug = FALSE) {
  # Check inputs
  stopifnot("snf(): Input argument *a* must be a polm object" =  inherits(a, 'polm'))

  # skip "zero" leading coefficients
  a = prune(a, tol = tol)
  # a0 = a

  # Dimensions of a
  d = unname(dim(a))
  m = d[1]
  n = d[2]

  stopifnot("snf(): Input argument *a* must be a non-empty polynomial matrix!" = m*n != 0)

  # initialize unimodular matrices. a = u s v, u_inv is
  u = polm(diag(m))
  u_inv = u
  v = polm(diag(n))
  v_inv = v

  i = 1
  iteration = 0
  while (i <= min(m,n)) {
    if (iteration > 100) stop('iteration maximum reached')

    # a is block diagonal
    # a[1:(i-1),1:(i-1)] is already in SNF form
    # handle the lower, right block: a[i:m, i:n]
    iteration = iteration + 1

    p = degree(a) # degree of each element (i,j) of the polynomial matrix
    if (debug) {
      message('snf: i=', i, ', iteration=', iteration)
      print(a, format = 'i|zj', digits = 2)
      print(p)
    }
    {# print(all.equal(a0, u %r% a %r% v))
    # print(prune(u %r% u_inv, tol = tol))
    # print(prune(v %r% v_inv, tol = tol))
    }
    # code zero entries as Inf
    p[p == -1] = Inf

    # a[i:m, i:n] is zero => a is in SNF form!
    if (all(is.infinite(p[i:m, i:n]))) {
      return(list(s = a, u = u, u_inv = u_inv, v = v, v_inv = v_inv))
    }

    # a[i,i] is non zero and all entries to the right and below the (i,i)-th element are zero
    if (is.finite(p[i,i]) &&
        all(is.infinite(p[iseq(i+1,m), i])) &&
        all(is.infinite(p[i, iseq(i+1,n)])) ) {

      # At most bottom-right element: Make monic and return
      if (i == min(m,n)) {
        c = unclass(a)[i, i, p[i,i]+1]
        a[i, ] = a[i, ] %/% c
        u[, i] = u[, i] * c
        u_inv[i, ] = u_inv[i, ] %/% c
        return(list(s = a, u = u, u_inv = u_inv, v = v, v_inv = v_inv))
      }

      # Check remainder a[,] %% a[i,i] ####
      # (%% is polynomial remainder, %/%  division, see ?Ops.ratm)
      ra = test_polm(dim = c(m,n), degree = -1)
      ra[(i+1):m, (i+1):n] = prune(a[(i+1):m, (i+1):n] %% a[i,i], tol = tol)
      rp = degree(ra)
      rp[rp == -1] = Inf

      # all remainders are zero => next step i -> i+1
      if (all(is.infinite(rp))) {

        # first make diagonal element monic
        c = unclass(a)[i, i, p[i,i]+1]
        a[i, ] = a[i, ] %/% c
        u[, i] = u[, i] * c
        u_inv[i, ] = u_inv[i, ] %/% c

        i = i+1
        next
      }

      # find element with minimal (remainder) degree
      c = which.min(apply(rp, MARGIN = 2, FUN = min))
      r = which.min(rp[, c])
      f = a[r,c] %/% a[i,i] # element-wise polynomial division

      # a[r,c] <- (a[r,c] %% a[i,i]) = (a[r,c] - f * a[i,i]), (%% gives remainder, %/% divides polynomial, and discards remainder)
      if (debug) {
        message('snf: a[', r, ',', c, '] <- a[', r, ',', c, '] %% a[i,i]\n')
      }

      {# Go through th code below with a diagonal matrix containing two different factors. First, add the (i,i) element to the zero element in row r, then use a column transformation (Euclidean algorithm) to subtract a factor times the (i,i)-element from the (r,c)-element

      # add i-th row to r-th row
      a[r,] = a[r,] + a[i,]
      u_inv[r,] = u_inv[r,] + u_inv[i,]
      u[,i] = u[,i] - u[,r]
      # substract f*(i-th column) from c-th column
      a[,c] = a[,c] - f * a[,i]
      v_inv[,c] = v_inv[,c] - f * v_inv[,i]
      v[i,] = v[i,] + f * v[c,]
      a = prune(a, tol = tol)
      }

      # next iteration
      next
    }


    if (i > 1) {
      diag(p)[1:(i-1)] = Inf
    }

    # find element with minimal degree
    c = which.min(apply(p, MARGIN = 2, FUN = min))
    r = which.min(p[, c])

    # bring this element to position (i,i)
    if (debug) {
      message('snf: a[i,i] <- a[', r, ',', c, ']\n')
    }

    rperm = 1:m
    rperm[c(i, r)] = c(r, i)
    cperm = 1:n
    cperm[c(i, c)] = c(c, i)
    a = a[rperm, cperm]
    p = p[rperm, cperm]
    u = u[, rperm]
    u_inv = u_inv[rperm, ]
    v = v[cperm, ]
    v_inv = v_inv[, cperm]

    # apply column purge
    out = purge_rc(a, pivot = c(i,i), direction = 'right',
                   monic = FALSE, permute = TRUE, tol = tol, debug = debug)
    a = out$h
    v = out$u %r% v
    v_inv = v_inv %r% out$u_inv

    # apply row purge
    # note: this may generate non zero elements in the i-th column
    out = purge_rc(a, pivot = c(i,i), direction = 'down',
                   monic = FALSE, permute = TRUE, tol = tol, debug = debug)
    a = out$h
    u = u %r% out$u
    u_inv = out$u_inv %r% u_inv


    # next iteration
  }

  return(list(s = a, u = u, u_inv = u_inv, v = v, v_inv = v_inv))
}

# internal function
# factorize a scalar polynomial into an stable and an unstable part
whf_scalar = function(a, tol = sqrt(.Machine$double.eps)) {
  a = prune(a, tol = tol)
  z = polyroot(unclass(a))
  a_f = a       # forward part (zeroes |z| < 1))
  a_b = polm(1) # backward part (zeroes |z| > 1))
  for (i in iseq(1,length(z))) {
    if (abs(z[i]) == 1) stop('unit roots are not allowed')
    if (abs(z[i]) > 1) {
      a_i = polm(c(-z[i], 1))
      a_f = a_f %/% a_i
      a_b = a_b * a_i
    }
  }
  a_f = prune(a_f, tol = tol)
  a_b = prune(a_b, tol = tol)
  if (is.complex(c(unclass(a_f), unclass(a_b)))) {
    stop('factors "a_f", "a_b" are complex')
  }
  return(list(a_f = a_f, a_b = a_b))
}

whf = function(a, right_whf = TRUE, tol = sqrt(.Machine$double.eps), debug = FALSE) {

  # Check inputs ####
  d = dim(a)
  stopifnot("whf(): Input *a* must be a polynomial matrix" = inherits(a, 'polm'),
            "whf(): Input *a* must be a square, non singular, polynomial matrix" = (d[1]*d[2]*d[3] != 0) && (d[1] == d[2]))
  m = d[1]

  # Left or right WHF?
  if (!right_whf){
    a = t(a)
  }

  # Obtain Smith normal form (in particular the diagonal matrix)
  snf = snf(a, tol = tol, debug = debug)

  {# cat('**************** SNF\n')
  # print(snf$s, digits = 3, format = 'i|zj')
  # print(snf$u, digits = 3, format = 'i|zj')
  # print(snf$v, digits = 3, format = 'i|zj')
  }

  # Factorize the diagonal matrix into stable and unstable part
  s_f = polm(diag(m))
  s_b = s_f
  for (i in (1:m)) {
    out = whf_scalar(snf$s[i,i])
    s_f[i,i] = out$a_f
    s_b[i,i] = out$a_b
  }

  {# cat('************** diagonal\n')
  # print(snf$s, digits = 2, format = 'i|zj')
  # print(s_b*s_f, digits = 2, format = 'i|zj')
  }

  ar = snf$u %r% s_f
  ab = s_b %r% snf$v

  if (debug) {
    message('whf:')
    cat('col degree Ar', degree(ar, 'columns'),'\n')
    cat('row degree Ab', degree(ar, 'rows'),'\n')
  }

  {# cat('************** Ar Ab\n')
  # print(a, digits = 2, format = 'i|zj')
  # print(ar * ab, digits = 2, format = 'i|zj')
  }

  # column reduction of ar
  out = col_reduce(ar, tol = tol, debug = debug)
  # print(out)
  ar = out$a
  ab = out$v %r% ab
  idx = out$col_degrees

  # try to 'balance' Ar and Ab
  l2norm = function(x) sqrt(sum(x^2))
  nr = apply(unclass(ar), MARGIN = 2, FUN = l2norm)
  nb = apply(unclass(ab), MARGIN = 1, FUN = l2norm)
  nn = sqrt(nb/nr)
# print(rbind(nr,nb,nn))
  ar = prune(ar %r% diag(nn, nrow = m, ncol = m), tol = tol)
  ab = prune(diag(nn^{-1}, nrow = m, ncol = m) %r% ab, tol = tol)
{# nr = apply(unclass(ar), MARGIN = 2, FUN = l2norm)
# nb = apply(unclass(ab), MARGIN = 1, FUN = l2norm)
# nn = sqrt(nb/nr)
# print(rbind(nr,nb,nn))
  }
  a0 = polm(diag(m))
  z = polm(c(0,1))
  for (i in (1:m)) {
    if (idx[i] > 0) a0[i,i] = z^idx[i]
  }

  # create the forward part Af
  # multiply the j-th column with z^(-idx[j]) => polynomial in z^(-1)
  ar_tmp = unclass(ar)
  af = array(0, dim = dim(ar_tmp))
  idx_max = max(idx)
  for (j in (1:m)) {
    af[,j,(1+idx_max-idx[j]):(1+idx_max)] = ar_tmp[,j,1:(idx[j]+1), drop = FALSE]
  }
  af = lpolm(af, min_deg = -idx_max)
  af = prune(af, tol = tol)

  # if left- WHF transform elements back
  if (!right_whf){
    af = t(af)
    ab = t(ab)
    ar = t(ar)
  }

  return(list(af = af, a0 = a0, ab = ab, ar = ar, idx = idx))
}



# Laurent polynomial transformations ####

polm2fwd = function(polm_obj){
  x = unclass(polm_obj)
  d = dim(x)
  if (d[3] > 0){
    return(lpolm(x[,,d[3]:1], min_deg = -(d[3]-1)))
  } else {
    return(lpolm(x, min_deg = 0))
  }
}

get_fwd = function(lpolm_obj){
  md = attr(lpolm_obj, which = "min_deg")
  if (md >=0){
    return(lpolm_obj)
  } else {
    x = unclass(lpolm_obj)[,,1:(-md), drop = FALSE]
    return(lpolm(x, min_deg = md))
  }
}

get_bwd = function(lpolm_obj){
  attr_l = attributes(lpolm_obj)
  min_deg = attr_l$min_deg
  d = attr_l$dim
  attributes(lpolm_obj) = NULL
  dim(lpolm_obj) = d
  if (d[3]+1 < -min_deg){
    # No non-negative coefficients: Empty lpolm
    return(lpolm(array(0, c(d[1],d[2],0)), min_deg = 0))
  } else if (min_deg >= 0){
    polm_offset = array(0, dim = c(dim(lpolm_obj)[1:2], min_deg))
    return(lpolm(dbind(d = 3, polm_offset, lpolm_obj), min_deg = 0))
  } else {
    return(lpolm(lpolm_obj[,,(-min_deg+1):d[3], drop = FALSE], min_deg = 0))
  }
}



# Helpers ##############################################################

as_txt_scalarpoly = function(coefs,
                             syntax = c('txt', 'TeX', 'expression'), x = 'z',
                             laurent = FALSE) {
  coefs = as.vector(coefs)
  if (!is.numeric(coefs)) {
    stop('"coefs" must be a numeric vector')
  }
  syntax = match.arg(syntax)

  if ((length(coefs) == 0) || all(coefs == 0)) {
    return('0')
  }

  p = length(coefs)-1
  if (laurent){
    powers = laurent:(laurent+p)
  } else {
    powers = (0:p)
  }

  # skip zero coefficients
  non_zero = (coefs != 0)
  coefs = coefs[non_zero]
  powers = powers[non_zero]

  # convert powers to character strings
  if (syntax == 'txt') {
    # x^k
    powers_txt = paste(x, '^', powers, sep = '')
  } else {
    # x^{k}
    powers_txt = paste(x, '^{', powers, '}', sep = '')
    # fmt = 'x^{k}' # never used according to RStudio
  }
  powers_txt[powers == 0] = ''
  powers_txt[powers == 1] = x
  powers = powers_txt

  signs = ifelse(coefs < 0, '- ', '+ ')
  signs[1] = ifelse(coefs[1] < 0, '-', '')

  # convert coefficients to character strings
  coefs = paste(abs(coefs))
  coefs[ (coefs == '1') & (powers != '') ] = ''

  if (syntax == 'expression') {
    mults = rep('*', length(coefs))
    mults[ (coefs == '') | (powers == '') ] = ''
  } else {
    mults = rep('', length(coefs))
  }

  txt = paste(signs, coefs, mults, powers, sep = '', collapse = ' ')
  return(txt)

}

as_txt_scalarfilter = function(coefs, syntax = c('txt', 'TeX', 'expression'),
                               x = 'z', t = 't') {
  coefs = as.vector(coefs)
  if (!is.numeric(coefs)) {
    stop('"coefs" must be a numeric vector')
  }
  syntax = match.arg(syntax)

  if ((length(coefs) == 0) || all(coefs == 0)) {
    return('0')
  }
  lags = (0:(length(coefs)-1))

  # skip zero coefficients
  non_zero = (coefs != 0)
  coefs = coefs[non_zero]
  lags = lags[non_zero]

  # convert lags to character strings
  if (syntax == 'TeX') {
    # x_{t-k}
    lags_txt = paste(x, '_{', t, '-', lags, '}', sep = '')
    lags_txt[lags == 0] = paste(x, '_{', t, '}', sep = '')
  } else {
    # x[t-k]
    lags_txt = paste(x, '[', t, '-', lags, ']', sep = '')
    lags_txt[lags == 0] = paste(x, '[', t, ']', sep = '')
  }
  lags = lags_txt

  signs = ifelse(coefs < 0, '- ', '+ ')
  signs[1] = ifelse(coefs[1] < 0, '-', '')

  # convert coefficients to character strings
  coefs = paste(abs(coefs))
  coefs[ (coefs == '1') & (lags != '') ] = ''

  if (syntax == 'expression') {
    mults = rep('*', length(coefs))
    mults[ (coefs == '') | (lags == '') ] = ''
  } else {
    mults = rep('', length(coefs))
  }
  txt = paste(signs, coefs, mults, lags, sep = '', collapse = ' ')
  return(txt)

}

as_tex_matrix = function(x) {
  if ( !is.matrix(x) ) stop('"x" must be a matrix')

  m = nrow(x)
  n = ncol(x)

  if (length(x) == 0) return('\\begin{pmatrix}\n\\end{pmatrix}')

  tex = '\\begin{pmatrix}\n'
  for (i in (1:m)) {
    tex = paste(tex, paste(x[i,], collapse = ' & '), '\\\\\n', sep = '  ')
  }
  tex = paste(tex, '\\end{pmatrix}', sep = '')
  return(tex)
}


as_tex_matrixpoly = function(coefs, x = 'z', as_matrix_of_polynomials = TRUE) {
  # only some basic checks
  if ( (!is.array(coefs)) || (length(dim(coefs)) != 3) || (!is.numeric(coefs)) ) {
    stop('"coefs" must be 3-dimensional numeric array')
  }

  d = dim(coefs)
  m = d[1]
  n = d[2]
  p = d[3] - 1

  if ((m*n) == 0) {
    return('\\begin{pmatrix}\n\\end{pmatrix}')
  }

  if ((m*n) == 1) {
    return(as_txt_scalarpoly(coefs, syntax = 'TeX', x = x))
  }

  if ((p < 0) || all(coefs == 0)) {
    return(as_tex_matrix(matrix(0, nrow = m, ncol = n)))
  }

  if (as_matrix_of_polynomials) {
    tex = apply(coefs, MARGIN = c(1,2), FUN = as_txt_scalarpoly,
                syntax = 'TeX', x = x)
    tex = as_tex_matrix(tex)
    return(tex)
  }

  # print as polynomial with matrix coefficients

  powers = (0:p)
  # coerce powers to character strings of the form x^{k}
  powers_txt = paste(x, '^{', powers, '}', sep = '')
  powers_txt[powers == 0] = ''
  powers_txt[powers == 1] = x
  powers = powers_txt

  tex = ''
  for (k in (0:p)) {
    a = matrix(coefs[,,k+1], nrow = m, ncol = n)
    if ( !all(a == matrix(0, nrow = m, ncol = n)) ) {
      # non-zero coefficient matrix
      if (tex != '' ) tex = paste(tex, '+\n')
      if ( (m == n) && all(a == diag(m)) ) {
        # coefficient matrix is identity matrix
        tex = paste(tex, ' I_{', m, '} ', powers[k+1], sep = '')
      } else {
        tex = paste(tex, as_tex_matrix(a), powers[k+1])
      }
    }
  }

  return(tex)
}


as_tex_matrixfilter = function(coefs, x = 'z', t = 't') {
  # only some basic checks
  if ( (!is.array(coefs)) || (length(dim(coefs)) != 3) || (!is.numeric(coefs)) ) {
    stop('"coefs" must be 3-dimensional numeric array')
  }

  d = dim(coefs)
  m = d[1]
  n = d[2]
  p = d[3] - 1

  if ((m*n) == 0) {
    return('\\begin{pmatrix}\n\\end{pmatrix}')
  }

  if ((m*n) == 1) {
    return(as_txt_scalarfilter(coefs, syntax = 'TeX', x = x, t = t))
  }

  if ((p < 0) || all(coefs == 0)) {
    tex = as_tex_matrix(matrix(0, nrow = m, ncol = n))
    return( paste(tex, ' ', x, '_{', t, '}', sep = '') )
  }

  lags = (0:p)
  # coerce lags to character strings of the form x_{t-k}
  lags_txt = paste(x, '_{', t, '-', lags, '}', sep = '')
  lags_txt[lags == 0] = paste(x, '_{', t, '}', sep = '')
  lags = lags_txt

  tex = ''
  for (k in (0:p)) {
    a = matrix(coefs[,,k+1], nrow = m, ncol = n)
    if ( !all(a == matrix(0, nrow = m, ncol = n)) ) {
      # non-zero coefficient matrix
      if (tex != '' ) tex = paste(tex, '+\n')
      if ( (m==n) && all(a == diag(m)) ) {
        # coefficient matrix is identity matrix
        tex = paste(tex, ' I_{', m, '} ', lags[k+1], sep = '')
      } else {
        tex = paste(tex, as_tex_matrix(a), lags[k+1])
      }
    }
  }

  return(tex)
}

print_3D = function(a, digits = NULL,
                    format = c('i|jz', 'i|zj', 'iz|j', 'zi|j', 'i|j|z','character'),
                    laurent = FALSE) {
  dim = dim(a)
  m = dim[1]
  n = dim[2]
  p = dim[3]
  # empty array -> do nothing
  if (min(dim) == 0) return(invisible(NULL))

  # a must have full 'dimnames'
  names = dimnames(a)
  inames = names[[1]]
  jnames = names[[2]]
  znames = names[[3]]

  # round
  if (!is.null(digits)) a = round(a, digits)

  format = match.arg(format)

  if (format == 'character') {

    # convert vector of coefficients to character representation of a polynomial
    a = apply(a, MARGIN = c(1,2), FUN = as_txt_scalarpoly, syntax = "txt", x = "z",
              laurent = laurent)

    # add column names (jnames)
    a = rbind(jnames, a)

    # add row names (inames)
    a = cbind( c('',inames), a)

    # right justify columns
    w = nchar(a)
    w = apply(w, MARGIN = 2, FUN = max)
    for (j in (1:(n+1))) {
      fmt = paste('%', w[j], 's', sep='')
      pad = function(s) { sprintf(fmt, s) }
      a[,j] = apply(a[,j,drop = FALSE], MARGIN = 1, FUN = pad)
    }

    # convert matrix a to a string
    a = apply(a, MARGIN = 1, FUN = paste, collapse = '  ')
    a = paste(a, collapse = '\n')
    cat(a,'\n')
  }


  if (format == 'i|jz') {
      # create a vector of the form
      # j[1],...,j[n],j[1],...,j[n],...
      jnames = rep(jnames, p)
      # create a vector of the form
      # z[1],'',...,'',z[2],'',...,'',
      if (n > 1) {
        znames = as.vector(rbind(znames,
                                 matrix('', nrow = n-1, ncol = p)))
      }

      dim(a) = c(m,n*p)
      rownames(a) = inames
      colnames(a) = paste(znames, jnames, sep = ' ')
      print(a)
  }

  if (format == 'i|zj') {
    # create a vector of the form
    # z[1],...,z[p],z[1],...,z[p],...
    znames = rep(znames, n)
    # create a vector of the form
    # j[1],'',...,'',j[2],'',...,'',
    if (p > 1) {
      jnames = as.vector(rbind(jnames,
                               matrix('', nrow = p-1, ncol = n)))
    }

    a = aperm(a, c(1,3,2))
    dim(a) = c(m,p*n)
    rownames(a) = inames
    colnames(a) = paste(jnames, znames, sep = ' ')
    print(a)
  }

  if (format == 'iz|j')  {
    # create a vector of the form
    # i[1],...,i[m],i[1],...,i[m],...
    inames = rep(inames, p)
    # create a vector of the form
    # z[1],'  ',...,'  ',z[2],'  ',...,'  ',
    if (m > 1) {
      znames = as.vector(rbind( znames,
                                matrix(' ', nrow = m-1, ncol = p)))
    }
    # right justify
    fmt = paste('%', max(nchar(znames)), 's', sep='')
    pad = function(s) { sprintf(fmt, s) }
    znames = as.vector(apply(matrix(znames, ncol = 1), MARGIN = 1, FUN = pad))

    a = aperm(a, c(1,3,2))
    dim(a) = c(m*p, n)
    rownames(a) = paste(znames, inames, sep = ' ')
    colnames(a) = jnames
    print(a)
  }

  if (format == 'zi|j')  {
    # create a vector of the form
    # z[1],...,z[p],z[1],...,z[p],...
    znames = rep(znames, m)
    # create a vector of the form
    # i[1],'  ',...,'  ',i[2],'  ',...,'  ',
    if (p > 1) {
      inames = as.vector(rbind( inames,
                                matrix(' ',nrow = p-1, ncol = m)))
    }
    # right justify
    fmt = paste('%', max(nchar(inames)), 's', sep='')
    pad = function(s) { sprintf(fmt, s) }
    inames = as.vector(apply(matrix(inames, ncol = 1), MARGIN = 1, FUN = pad))

    a = aperm(a, c(3,1,2))
    dim(a) = c(p*m, n)
    rownames(a) = paste(inames, znames, sep = ' ')
    colnames(a) = jnames
    print(a)
  }

  if (format == 'i|j|z') {
    # the last case 'i|j|z' just uses the R default print of 3D array
    print(a)
  }

  return(invisible(NULL))
}

# print.___() for rationalmatrices objects ####

NULL


print.lpolm = function(x, digits = NULL,
        format = c('i|jz', 'i|zj', 'iz|j', 'zi|j', 'i|j|z','character'), ...) {
  if (!is.null(digits)) {
    digits = as.vector(as.integer(digits))[1]
  }
  format = match.arg(format)

  a = unclass(x)
  min_deg = attr(x, which = "min_deg")
  attr(a, 'min_deg') = NULL # remove 'min_deg' attribute
  m = dim(a)[1]
  n = dim(a)[2]
  p = dim(a)[3]-1+min_deg

  cat('( ', m, ' x ', n,' ) Laurent polynomial matrix with degree <= ', p,
      ', and minimal degree >= ', min_deg, '\n', sep = '')
  if ((m*n*(dim(a)[3])) == 0) {
    return(invisible(x))
  }

  if ((format == 'character') && (is.complex(a))) {
    stop(paste('the format option "character" is only implemented',
               'for Laurent polynomials with real coefficients'))
  }

  # use the above defined internal function print_3D
  dimnames(a) = list(paste('[', 1:m, ',]', sep = ''),
                     paste('[,', 1:n, ']', sep = ''),
                     paste('z^', min_deg:p, sep = ''))
  print_3D(a, digits, format, laurent = min_deg)

  invisible(x)
}

print.polm = function(x, digits = NULL,
        format = c('i|jz', 'i|zj', 'iz|j', 'zi|j', 'i|j|z','character'), ...) {
  if (!is.null(digits)) {
    digits = as.vector(as.integer(digits))[1]
  }
  format = match.arg(format)

  a = unclass(x)
  m = dim(a)[1]
  n = dim(a)[2]
  p = dim(a)[3]-1

  cat('(',m,'x',n,') matrix polynomial with degree <=', p,'\n')
  if ((m*n*(p+1)) == 0) {
    return(invisible(x))
  }

  if ((format == 'character') && (is.complex(a))) {
    stop(paste('the format option "character" is only implemented',
               'for polynomials with real coefficients'))
  }

  # use the above defined internal function print_3D
  dimnames(a) = list(paste('[', 1:m, ',]', sep = ''),
                     paste('[,', 1:n, ']', sep = ''),
                     paste('z^',0:p, sep = ''))
  print_3D(a, digits, format)

  invisible(x)
}


print.lmfd = function(x, digits = NULL, format = c('i|jz', 'i|zj', 'iz|j', 'zi|j', 'i|j|z','character'), ...) {
  if (!is.null(digits)) {
    digits = as.vector(as.integer(digits))[1]
  }
  format = match.arg(format)

  d = attr(x, 'order')
  m = d[1]
  n = d[2]
  p = d[3]
  q = d[4]

  cat('( ', m, ' x ', n,' ) left matrix fraction description a^(-1)(z) b(z) with degrees (p = ',
      p, ', q = ', q, ')\n', sep = '')

  if ((format == 'character') && (is.complex(unclass(x)))) {
    stop('the format option "character" is only implemented for LMFDs with real coefficients')
  }

  if ((m*m*(p+1)) > 0) {
    cat('left factor a(z):\n')

    a = unclass(x$a)

    # use the above defined internal function print_3D
    dimnames(a) = list(paste('[', 1:m, ',]', sep = ''),
                       paste('[,', 1:m, ']', sep = ''),
                       paste('z^',0:p, sep = ''))
    print_3D(a, digits, format)
  }

  if ((m*n*(q+1)) > 0) {
    cat('right factor b(z):\n')

    a = unclass(x$b)

    # use the above defined internal function print_3D
    dimnames(a) = list(paste('[', 1:m, ',]', sep = ''),
                       paste('[,', 1:n, ']', sep = ''),
                       paste('z^',0:q, sep = ''))
    print_3D(a, digits, format)
  }

  invisible(x)
}

print.rmfd = function(x, digits = NULL, format = c('i|jz', 'i|zj', 'iz|j', 'zi|j', 'i|j|z','character'), ...) {
  if (!is.null(digits)) {
    digits = as.vector(as.integer(digits))[1]
  }
  format = match.arg(format)

  d = attr(x, 'order')
  m = d[1]
  n = d[2]
  p = d[3]
  q = d[4]

  cat('( ', m, ' x ', n,' ) right matrix fraction description d(z) c^(-1)(z) with degrees deg(c(z)) = p = ',
      p, ', deg(d(z)) = q = ', q, '\n', sep = '')

  if ((format == 'character') && (is.complex(unclass(x)))) {
    stop('the format option "character" is only implemented for RMFDs with real coefficients')
  }


  if ((m*n*(q+1)) > 0) {
    cat('left factor d(z):\n')

    d = unclass(x$d)

    # use the above defined internal function print_3D
    dimnames(d) = list(paste('[', 1:m, ',]', sep = ''),
                       paste('[,', 1:n, ']', sep = ''),
                       paste('z^',0:q, sep = ''))
    print_3D(d, digits, format)
  }

  if ((n*n*(p+1)) > 0) {
    cat('right factor c(z):\n')

    c = unclass(x$c)

    # use the above defined internal function print_3D
    dimnames(c) = list(paste('[', 1:n, ',]', sep = ''),
                       paste('[,', 1:n, ']', sep = ''),
                       paste('z^',0:p, sep = ''))
    print_3D(c, digits, format)
  }

  invisible(x)
}


print.stsp = function(x, digits = NULL, ...) {
  if (!is.null(digits)) {
    digits = as.vector(as.integer(digits))[1]
  }

  d = attr(x, 'order')
  m = d[1]
  n = d[2]
  s = d[3]

  cat('statespace realization [', m, ',', n, '] with s = ', s, ' states\n', sep = '')

  a = unclass(x)
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
  if (m > 0) xnames = paste('x[',1:m,']',sep = '')
  unames = character(n)
  if (n > 0) unames = paste('u[',1:n,']',sep = '')

  rownames(a) = c(snames, xnames)
  colnames(a) = c(snames, unames)
  print(a)

  invisible(x)
}


print.pseries = function(x, digits = NULL, format = c('i|jz', 'i|zj', 'iz|j', 'zi|j', 'i|j|z'), ...) {
  if (!is.null(digits)) {
    digits = as.vector(as.integer(digits))[1]
  }
  format = match.arg(format)

  a = unclass(x)
  m = dim(a)[1]
  n = dim(a)[2]
  lag.max = dim(a)[3]-1

  cat('(',m,'x',n,') impulse response with maximum lag =', lag.max,'\n')
  if ((m*n*(lag.max+1)) == 0) {
    return(invisible(x))
  }

  # use the above defined internal function print_3D
  dimnames(a) = list(paste('[', 1:m, ',]', sep = ''),
                     paste('[,', 1:n, ']', sep = ''),
                     paste('lag=',0:lag.max, sep = ''))
  print_3D(a, digits, format)

  invisible(x)
}

print.zvalues = function(x, digits = NULL, format = c('i|jz', 'i|zj', 'iz|j', 'zi|j', 'i|j|z'), ...) {
  if (!is.null(digits)) {
    digits = as.vector(as.integer(digits))[1]
  }
  format = match.arg(format)

  z = attr(x, 'z')
  n.z = length(z)

  a = unclass(x)
  m = dim(a)[1]
  n = dim(a)[2]
  attr(a, 'z') = NULL # remove 'z' attribute

  cat('(',m,'x',n,') frequency response\n')
  if ((m*n*n.z) == 0) {
    return(invisible(x))
  }

  # use the above defined internal function print_3D
  if ((format == 'i|jz') || (format == 'i|zj')) {
    dimnames(a) = list(paste('[', 1:m, ',]', sep = ''),
                       paste('[,', 1:n, ']', sep = ''),
                       paste('z[',1:n.z, ']', sep = ''))
  } else {
    dimnames(a) = list(paste('[', 1:m, ',]', sep = ''),
                       paste('[,', 1:n, ']', sep = ''),
                       paste('z=', round(z,3), sep = ''))
  }
  print_3D(a, digits, format)

  invisible(x)
}
# Power series parameters ######################################################

pseries = function(obj, lag.max, ...){
  UseMethod("pseries", obj)
}

pseries.default = function(obj, lag.max = 5, ...) {
  # try to coerce to polm
  obj = try(polm(obj), silent = TRUE)
  if (inherits(obj, 'try-error')) stop('could not coerce argument "obj" to "polm" object')

  lag.max = as.integer(lag.max[1])
  if (lag.max < 0) {
    stop('lag.max must be non-negative.')
  }

  obj = unclass(obj)
  m = dim(obj)[1]
  n = dim(obj)[2]
  p = dim(obj)[3] - 1

  ir = array(0, dim = unname(c(m,n,lag.max + 1)))

  if (p > -1) ir[,,1:min(lag.max+1, p+1)] = obj[,, 1:min(lag.max+1,p+1)]

  ir = structure(ir, class = c('pseries','ratm'))
  return(ir)
}

pseries.polm = function(obj, lag.max = 5, ...) {
  lag.max = as.integer(lag.max[1])
  if (lag.max < 0) {
    stop('lag.max must be non-negative.')
  }

  obj = unclass(obj)
  m = dim(obj)[1]
  n = dim(obj)[2]
  p = dim(obj)[3] - 1

  ir = array(0, dim = unname(c(m,n,lag.max + 1)))

  if (p > -1) ir[,,1:min(lag.max+1, p+1)] = obj[,, 1:min(lag.max+1,p+1)]

  ir = structure(ir, class = c('pseries','ratm'))
  return(ir)
}

pseries.lmfd = function(obj, lag.max = 5, ...) {
  lag.max = as.integer(lag.max[1])
  if (lag.max < 0) {
    stop('lag.max must be non-negative.')
  }

  d = attr(obj, 'order')
  m = d[1]
  n = d[2]
  p = d[3]
  q = d[4]

  ir = array(0, dim = unname(c(m,n,lag.max +1)))

  if ((m*n*(q+1)) == 0) {
    ir = structure(ir, class = c('pseries','ratm'))
    return(ir)
  }

  ab = unclass(obj)
  a = ab[,1:(m*(p+1)), drop = FALSE]
  b = ab[,(m*(p+1)+1):(m*(p+1)+n*(q+1)), drop = FALSE]

  # check a(z)
  if (p < 0) stop('left factor "a(z)" is not a valid polynomial matrix')

  a0 = a[,1:m, drop = FALSE]
  junk = try(solve(a0))
  if (inherits(junk, 'try-error')) stop('left factor "a(0)" is not invertible')

  # b[,,i] -> a0^{-1} b[,,i]
  b = solve(a0, b)
  dim(b) = c(m, n, q+1)

  ir[ , , 1:min(lag.max+1, q+1)] = b[ , , 1:min(lag.max+1, q+1)]

  if ((p == 0) || (lag.max == 0)) {
    ir = structure(ir, class = c('pseries','ratm'))
    return(ir)
  }

  # a[,,i] -> a0^{-1} a[,,i]
  a = solve(a0, a[,(m+1):ncol(a), drop = FALSE])
  dim(a) = c(m, m, p)

  for (lag in (1:lag.max)) {
    for (i in (1:min(p,lag))) {
      ir[,,lag+1] = matrix(ir[,,lag+1], nrow = m, ncol = n) -
        matrix(a[,,i], nrow = m, ncol = m) %*% matrix(ir[,,lag+1-i], nrow = m, ncol = n)
    }
  }

  ir = structure(ir, class = c('pseries','ratm'))
  return(ir)
}


pseries.rmfd = function(obj, lag.max = 5, ...) {
  lag.max = as.integer(lag.max[1])
  if (lag.max < 0) {
    stop('lag.max must be non-negative.')
  }

  d = attr(obj, 'order')
  m = d[1]
  n = d[2]
  p = d[3]
  q = d[4]

  ir = array(0, dim = unname(c(m,n,lag.max +1)))

  if ((m*n*(q+1)) == 0) {
    ir = structure(ir, class = c('pseries','ratm'))
    return(ir)
  }

  # check c(z)
  if (p < 0) stop('right factor "c(z)" is not a valid polynomial matrix')

  # Separate c(z) and d(z)
  cd = unclass(obj)
  c = cd[1:(n*(p+1)), , drop = FALSE]
  d = cd[(n*(p+1)+1):(n*(p+1)+m*(q+1)), , drop = FALSE]

  c0 = matrix(c[1:n, , drop = FALSE], nrow = n, ncol = n)
  c0inv = tryCatch(solve(c0),
                   error = function(cnd) stop(' "c(0)" is not invertible'))

  c = c %*% c0inv
  d = d %*% c0inv

  # Convert d(z) to an array
  d = t(d)
  dim(d) = c(n, m, q+1)
  d = aperm(d, perm = c(2,1,3))

  # Initialize impulse response with d(z)
  ir[ , , 1:min(lag.max+1, q+1)] = d[ , , 1:min(lag.max+1, q+1)]

  if ((p == 0) || (lag.max == 0)) {
    ir = structure(ir, class = c('pseries','ratm'))
    return(ir)
  }

  # Convert c(z) an array and discard zero-lag coefficient (equal to I_q)
  c = c[(n+1):nrow(c),, drop = FALSE]
  c = t(c)
  dim(c) = c(n, n, p)
  c = aperm(c, perm = c(2,1,3))

  # Calculate impulse response
  for (lag in (1:lag.max)) {
    for (i in (1:min(p,lag))) {
      ir[,,lag+1] = matrix(ir[,,lag+1], nrow = m, ncol = n) - matrix(ir[,,lag+1-i], nrow = m, ncol = n) %*% matrix(c[,,i], nrow = n, ncol = n)
    }
  }

  ir = structure(ir, class = c('pseries','ratm'))
  return(ir)
}



pseries.stsp = function(obj, lag.max = 5, ...) {
  lag.max = as.integer(lag.max[1])
  if (lag.max < 0) {
    stop('lag.max must be non-negative.')
  }

  d = attr(obj, 'order')
  m = d[1]
  n = d[2]
  s = d[3]

  ir = array(0, dim = unname(c(m,n,lag.max +1)))

  if ((m*n) == 0) {
    ir = structure(ir, class = c('pseries','ratm'))
    return(ir)
  }

  ABCD = unclass(obj)
  A = ABCD[iseq(1,s), iseq(1,s), drop = FALSE]
  B = ABCD[iseq(1,s), iseq(s+1,s+n), drop = FALSE]
  C = ABCD[iseq(s+1,s+m), iseq(1,s), drop = FALSE]
  D = ABCD[iseq(s+1,s+m), iseq(s+1,s+n), drop = FALSE]

  ir = array(0, dim = unname(c(m,n,lag.max + 1)))
  ir[,,1] = D

  if ((s == 0) || (lag.max == 0)) {
    ir = structure(ir, class = c('pseries','ratm'))
    return(ir)
  }

  for (lag in (1:lag.max)) {
    ir[,,lag+1] = C %*% B
    B = A %*% B
  }

  ir = structure(ir, class = c('pseries','ratm'))
  return(ir)
}

pseries.lpolm = function(obj, ...) {
  min_deg = attr(obj, which = "min_deg")
  if (min_deg >= 0){
    obj = lpolm(obj) # returns a polm object
    return(pseries(obj))
  } else {
    stop("A lpolm object cannot be coerced to a pseries object.")
  }
}

# reflect_poles_zeroes.R

## //' \ifelse{html}{\figure{internal_Rcpp.svg}
## {options: alt='Internal (Rcpp) function'}}{\strong{Internal (Rcpp)} function}

polm_div = function(a, b) {
  # a,b must be 'polm'-objects
  if ((!inherits(a, 'polm')) || (!inherits(b, 'polm'))) {
    stop('The input parameter "a", "b" must be "polm" objects!')
  }
  a = unclass(a)
  dim.a = dim(a)
  b = unclass(b)
  dim.b = dim(b)

  if ((length(dim.a) !=3) || (length(dim.b) != 3) ||
      (dim.a[2] != dim.b[1]) || (dim.b[1] != dim.b[2]) || (dim.b[1] == 0)) {
    stop('The polynomial matrices a, b are not compatible!')
  }

  m = dim.a[1]
  n = dim.a[2]
  p = dim.a[3]-1
  q = dim.b[3]-1

  if (q > p) {
    c = polm(array(0, dim = c(m, n, 1)))
    d = polm(a)
    return(list(qu = c, rem = d))
  }
  if (q == -1) {
    stop('The polynomial matrix b is zero!')
  }
  a = matrix(a, nrow = m, ncol = n*(p+1))
  b = matrix(b, nrow = n, ncol = n*(q+1))
  c = matrix(0, nrow = m, ncol = n*((p-q)+1))

  inv_bq = try(solve(b[ , (q*n+1):((q+1)*n), drop = FALSE]))
  if (inherits(inv_bq, 'try-error')) {
    # print(b[ , (q*n+1):((q+1)*n), drop = FALSE])
    stop('The leading coefficient of the polynomial matrix b is singular!')
  }

  for (i in ((p-q):0)) {
    # print(i)
    c[ , (i*n+1):((i+1)*n)] =
      a[ , ((q+i)*n+1):((q+i+1)*n), drop = FALSE] %*% inv_bq
    if (q > 0) {
      # print(a[ , (i*n+1):((i+q)*n), drop = FALSE])
      # print(c[ , (i*n+1):((i+1)*n), drop = FALSE])
      # print(b[, 1:(q*n), drop = FALSE])
      a[ , (i*n+1):((i+q)*n)] = a[ , (i*n+1):((i+q)*n), drop = FALSE] -
        c[ , (i*n+1):((i+1)*n), drop = FALSE] %*% b[, 1:(q*n), drop = FALSE]
    }
  }

  c = polm(array(c, dim = c(m,n,p-q+1)))
  if (q > 0) {
    d = polm(array(a[,1:(q*n)], dim = c(m,n,q)))
  } else {
    d = polm(array(0, dim = c(m,n,1)))
  }

  return(list(qu = c, rem = d))
}


blaschke = function(alpha) {
  alpha = as.vector(alpha)[1]
  BM = lmfd(polm(c(-alpha,1)), polm(c(1, -Conj(alpha))))
  return(BM)
}

blaschke2 = function(alpha, w = NULL, tol = 100*.Machine$double.eps) {
  alpha = as.vector(alpha)[1]
  if (Im(alpha) == 0) {
    stop('"alpha" must be complex with a non zero imaginary part!')
  }

  # Univariate case
  if (is.null(w)) {
    BM = lmfd(polm(c(Mod(alpha)^2, -2*Re(alpha), 1)),
              polm(c(1, -2*Re(alpha), Mod(alpha)^2)))
    return(BM)
  }

  # root 'alpha' is on the unit circle: simply return the identity!
  if (abs( Mod(alpha) - 1 ) < tol) {
    return(lmfd(polm(diag(2)), polm(diag(2))))
  }

  # w must be two dimensional complex vector and
  # w and Conj(w) must be linearly independent
  if (Mod(alpha) < 1) {
    # root alpha is inside the unit circle
    # a(z) = (-A + Iz)
    A = try(t(solve( t( cbind(w, Conj(w)) ),
                     t( cbind(alpha*w, Conj(alpha*w))) )), silent = TRUE)
    if (inherits(A, 'try-error')) {
      stop('w and Conj(w) are linearly dependent!')
    }
    A = Re(A)
    Gamma0 = lyapunov(t(A), diag(2), non_stable = 'ignore')
    # b(z) = (I - Bz)T^{-1}
    B = solve(Gamma0, t(A) %*% Gamma0)
    T = chol(Gamma0 - t(B) %*% Gamma0 %*% B)

    BM = lmfd(polm(array(cbind(-A, diag(2)), dim = c(2,2,2))),
              polm(array(cbind(diag(2), -B), dim = c(2,2,2))) %r% solve(T))
    return(BM)

  } else {
    # root alpha is outside the unit circle
    # a(z) = (I - Az)
    A = try( t(solve( t( cbind(alpha*w, Conj(alpha*w)) ),
                      t( cbind(w, Conj(w))) )), silent = TRUE)
    if (inherits(A, 'try-error')) {
      stop('w and Conj(w) are linearly dependent!')
    }
    A = Re(A)
    Gamma0 = lyapunov(t(A), diag(2), non_stable = 'ignore')
    # b(z) = (B - Iz)T^{-1}
    B = solve(Gamma0, t(A) %*% Gamma0)
    T = chol(Gamma0 - t(B) %*% Gamma0 %*% B)

    BM = lmfd(polm(array(cbind(diag(2), -A), dim = c(2,2,2))),
              polm(array(cbind(B, -diag(2)), dim = c(2,2,2))) %r% solve(T))
    return(BM)
  }
}

make_allpass = function(A, B) {
  m = ncol(B)
  s = nrow(B)

  P = lyapunov(A, B %*% t(B))
  C = -t(solve(A %*% P, B))
  #
  M = diag(m) - t(solve(A, B)) %*% t(C)
  # make sure that M is symmetric
  M = (M + t(M))/2
  D = solve(t(chol(M)))
  C = D %*% C
  # print(D %*% M %*% t(D))

  x = stsp(A, B, C, D)
  return(x)
}


reflect_zeroes = function(x, zeroes, ...) {
  UseMethod("reflect_zeroes")
}

reflect_zeroes.polm = function(x, zeroes, tol = sqrt(.Machine$double.eps),
                               check_zeroes = TRUE, ...) {
  if (!is.numeric(x)) {
    stop('The argument "x" must be a polynomial with real coefficients.')
  }
  d = dim(x)
  if (d[1] != d[2]) {
    stop('The argument "x" must be a square polynomial matrix (in polm form).')
  }
  m = d[1]

  zeroes = as.vector(zeroes)
  k = length(zeroes)
  if (k == 0) {
    # nothing to do
    return(x)
  }

  if (check_zeroes) {
    z = zeroes(x, tol = tol, print_message = FALSE)
    zz = c(zeroes, Conj(zeroes[Im(zeroes) != 0]))
    if (length(z) < length(zz)) {
      stop(paste('The polynomial "x" has less zeroes than specified',
                 'in the argument "zeroes"!'))
    }
    j = match_vectors(zz, z)
    if (max(abs(zz - z[j])) > tol) {
      print(cbind(zz, z[j], zz - z[j]))
      stop(paste('The specified zeroes do not match',
                 'with the zeroes of the polynomial matrix "x"!'))
    }
  }

  # should we sort/order the zeroes?
  for (i in (1:k)) {
    z0 = zeroes[i]
    # cat('zero:', z0, '\n')

    if (Im(z0) == 0) {
      # real zero
      z0 = Re(z0)
      x0 = zvalue(x, z0)

      U = svd(x0, nu = 0, nv = m)$v
      x = x %r% U
      B = blaschke(z0)
      x[, m] = (polm_div(x[, m], B$a)$qu ) %r% B$b
    } else {
      # complex zero
      x0 = zvalue(x, z0)

      # flip z0
      U = svd(x0, nu = 0, nv = m)$v
      x = x %r% U
      B = blaschke(z0)
      # print(zvalue(x, z0)[, m])
      x[, m] = (polm_div(x[, m], B$a)$qu) %r% B$b

      # flip conjugate zeroe Conj(z0) and the corresponding null vector Conj(w)
      w0 = t(Conj(U)) %*% Conj(U[, m])
      w0[m] = w0[m] / zvalue(B, Conj(z0))

      U = svd(w0, nu = m, nv = 0)$u
      x = x %r% U
      B = blaschke(Conj(z0))
      # print(zvalue(x, Conj(z0))[, 1])
      x[, 1] = (polm_div(x[, 1], B$a)$qu) %r% B$b

      # make real!
      # since we use x(z0) w0 = 0 <=> x(Conj(z0)) Conj(w0) = 0
      x0 = zvalue(x, 0)
      L = t(chol(Re(x0 %*% t(Conj(x0)))))
      x = x %r% solve(x0, L)
      x = Re(x)
    }
  }
  return(x)
}


reflect_zeroes.lmfd = function(x, zeroes, tol = sqrt(.Machine$double.eps),
                               check_zeroes = TRUE, ...) {
  if (!is.numeric(x)) {
    stop('The argument "x" must be a rational matrix with real coefficients.')
  }
  d = dim(x)
  if (d[1] != d[2]) {
    stop('The argument "x" must be a square polynomial matrix (in polm form).')
  }
  m = d[1]

  zeroes = as.vector(zeroes)
  k = length(zeroes)

  if (k == 0) {
    # nothing to do
    return(x)
  }

  x = lmfd(a = x$a,
           b = reflect_zeroes(x$b, zeroes, check_zeroes = check_zeroes,
                              tol = tol))
  return(x)
}


reflect_zeroes.stsp = function(x, zeroes, tol = sqrt(.Machine$double.eps), ...) {
  # check inputs ....
  if (!is.numeric(x)) {
    stop('The argument "x" must be a polynomial with real coefficients.')
  }
  d = dim(x)
  if (d[1] != d[2]) {
    stop('argument "x" must be a square rational matrix (in stsp form).')
  }

  zeroes = as.vector(zeroes)
  k = length(zeroes)
  if (k == 0) {
    # nothing to do
    return(x)
  }

  if (min(abs(abs(zeroes) - 1)) < tol) {
    stop('one of the selected zeroes has modulus close to one.')
  }

  # append complex conjugates
  zeroes = c(zeroes, Conj(zeroes[Im(zeroes) != 0]))

  A = x$A
  B = x$B
  C = x$C
  D = x$D
  s = nrow(B)
  m = ncol(B)

  # transform (A - BD^{-1} C) to upper block diagonal matrix,
  # where the top block corresponds to the selected zeroes
  # schur decomposition of (A - BD^{-1} C)
  out = try(schur(A - B %*% solve(D, C), 1/zeroes))
  if (inherits(out, 'try-error')) stop('ordered schur decomposition failed.')

  # (A - B D^{-1} C) = U S U'
  A = t(out$U) %*% A %*% out$U
  B = t(out$U) %*% B
  C = C %*% out$U
  k = out$k

  i1 = (1:k)
  i2 = iseq((k+1), s)

  # create all-pass function with given A,C
  U = t(make_allpass(t(out$S[i1, i1, drop = FALSE]),
                     t(solve(D, C[, i1,drop = FALSE]))))

  Ah = rbind( cbind(A[i2, i2, drop = FALSE], A[i2, i1, drop = FALSE]),
              cbind(A[i1, i2, drop = FALSE], A[i1, i1, drop = FALSE]) )
  Bh = rbind( B[i2, , drop = FALSE] %*% U$D, B[i1, , drop = FALSE] %*% U$D + U$B )
  Ch = cbind( C[, i2, drop = FALSE], C[, i1, drop = FALSE] )
  Dh = D %*% U$D

  return(stsp(Ah, Bh, Ch, Dh))
}



reflect_poles = function(x, poles, ...) {
  UseMethod("reflect_poles")
}

reflect_poles.stsp = function(x, poles, tol = sqrt(.Machine$double.eps), ...) {
  # check inputs ....
  if (!is.numeric(x)) {
    stop('The argument "x" must be a polynomial with real coefficients.')
  }
  d = dim(x)
  if (d[1] != d[2]) {
    stop('argument "x" must be a square rational matrix (in stsp form).')
  }

  poles = as.vector(poles)
  k = length(poles)
  if (k == 0) {
    # nothing to do
    return(x)
  }

  if (min(abs(abs(poles) - 1)) < tol) {
    stop('one of the selected poles has modulus close to one.')
  }

  # append complex conjugates
  poles = c(poles, Conj(poles[Im(poles) != 0]))

  A = x$A
  B = x$B
  C = x$C
  D = x$D
  s = nrow(B)
  m = ncol(B)

  # transform A matrix to lower block diagonal matrix,
  # where the top diagonal block corresponds to the selected poles!
  # ordered schur decomposition of A'
  out = try(schur(t(A), 1/poles))
  if (inherits(out, 'try-error')) stop('ordered schur decomposition failed.')

  # A' = U S U' => A = U S' U'
  A = t(out$S)
  B = t(out$U) %*% B
  C = C %*% out$U
  k = out$k

  i1 = (1:k)
  i2 = iseq((k+1), s)

  # create all-pass function
  U = make_allpass(A[i1, i1, drop = FALSE], B[i1,,drop = FALSE])
  U = U^{-1}

  Ah = rbind( cbind(A[i2, i2, drop = FALSE],
                    A[i2, i1, drop = FALSE] + B[i2, ,drop = FALSE] %*% U$C),
              cbind(matrix(0, nrow = k, ncol = s-k), U$A) )
  Bh = rbind( B[i2, , drop = FALSE] %*% U$D, U$B )
  Ch = cbind( C[, i2, drop = FALSE], C[ , i1, drop = FALSE] + D %*% U$C )
  Dh = D %*% U$D

  return(stsp(Ah, Bh, Ch, Dh))
}

reflect_poles.rmfd = function(x, poles, tol = sqrt(.Machine$double.eps),
                              check_poles = TRUE, ...) {
  # check inputs ....
  if (!is.numeric(x)) {
    stop('The argument "x" must be a polynomial with real coefficients.')
  }
  d = dim(x)
  if (d[1] != d[2]) {
    stop('argument "x" must be a square rational matrix (in stsp form).')
  }

  poles = as.vector(poles)
  k = length(poles)
  if (k == 0) {
    # nothing to do
    return(x)
  }

  c = x$c
  d = x$d

  c = t(reflect_zeroes(t(c), poles, check_zeroes = check_poles,
                              tol = tol))
  return(rmfd(c = c, d= d))
}

roots_as_list = function(roots, tol = sqrt(.Machine$double.eps)) {
  roots = as.vector(roots)
  if (is.numeric(roots)) {
    return(as.list(roots))
  }
  if (is.complex(roots)) {
    p = length(roots)
    if (p == 0) return(as.list(roots))

    j = match_vectors(roots, Conj(roots))
    i = (1:p)
    if (max(abs(roots[i] - Conj(roots[j]))) > tol) {
      stop('could not match pairs of complex conjugate roots')
    }

    # if roots[k] is a "real" root,
    #   then j[k] = k should hold
    # if roots[k1], roots[k2] is a pair of complex conjugate roots,
    #   then j[k1] = k2 and j[k2] = k1 should hold.
    # Together this means jj = j[j] = (1:p) must hold!
    jj = j[j]

    if (any(jj != i)) {
      stop('could not match pairs of complex conjugate roots')
    }

    # make sure that "real" roots have imaginary part = 0
    # and that pairs of roots are complex conjugates of each other
    roots = (roots + Conj(roots[j]))/2

    # skip roots with negative imaginary part
    roots = roots[Im(roots) >= 0]
    # order by imaginary part, => real roots come first
    roots = roots[order(Im(roots), Re(roots))]

    # convert to list
    roots = as.list(roots)

    # convert roots with zero imaginary part to 'numeric'
    roots = lapply(roots, function(x) ifelse(Im(x) > 0, x, Re(x)))

    return(roots)
  }
  stop('"roots" must be a vector of class "numeric" or "complex"!')
}
# schur.R
#
# Schur decomposition and related tools and methods

schur = function(A, select = NULL, tol = sqrt(.Machine$double.eps)) {
  if ( !is.numeric(A) || !is.matrix(A) || (nrow(A) != ncol(A)) ) {
    stop('input "A" must be a square, complex or real valued matrix.')
  }
  if ( any(!is.finite(A)) ) stop('input "A" contains missing or infinite elements.')

  m = nrow(A)
  if (m == 0) stop('input "A" has zero rows/columns.')

  out = QZ::qz.dgees(A)
  if (out$INFO != 0) stop('Schur decomposition failed. Error code of "dgees.f": ', out$INFO)
  lambda = complex(real = out$WR, imaginary = out$WI) # eigenvalues

  selected = logical(m)
  if (!is.null(select)) {
    if (is.character(select)) {
      selected = switch(select,
                        iuc = (abs(lambda) < 1),
                        ouc = (abs(lambda) > 1),
                        lhp = (out$WR < 0),
                        rhp = (out$WR > 0),
                        real = (out$WI == 0),
                        cplx = (out$WI != 0))
    } else {
      select = as.vector(select)
      if ( length(select) > m ) stop('The input vector "select" has more than ', m, ' entries.')
      if (length(select) > 0) {
        C = abs(matrix(select, nrow = length(select), ncol = m) -
                  matrix(lambda, nrow = length(select), ncol = m, byrow = TRUE))
        match = munkres(C)
        if ( match$c > length(select)*tol ) stop('could not match the target eigenvalues to the eigenvalues of "A"!')
        selected[match$a[,2]] = TRUE

        # make sure that complex conjugated pairs are both selected
        i = which( (out$WI > 0) & (selected) ) # positive imaginary parts come first
        selected[i+1] = TRUE
        i = which( (out$WI < 0) & (selected) ) # positive imaginary parts come first
        selected[i-1] = TRUE
      }
    }
  }
  k = sum(selected)
  if ((k > 0) && (k < m)) {
    # reorder diagonal blocks
    out = QZ::qz.dtrsen(out$T, out$Q, selected, job = "N", want.Q = TRUE)
    if (out$INFO != 0) stop('Reordering of Schur decomposition failed. Error code of "dtrsen.f": ', out$INFO)
    lambda = complex(real = out$WR, imaginary = out$WI) # eigenvalues
  }
  return(list(S = out$T, U = out$Q, lambda = lambda, k = k))
}
# str.____ methods ##############################################################

str.lpolm = function(object, ...) {
  d = dim(unclass(object))
  min_deg = attr(object, "min_deg")
  cat('( ',d[1],' x ',d[2],' ) Laurent polynomial matrix with degree <= ',
      d[3]-1+min_deg, ', and minimum degree >= ', min_deg, '\n', sep = '')
  return(invisible(NULL))
}

str.polm = function(object, ...) {
  d = dim(unclass(object))
  cat('( ',d[1],' x ',d[2],' ) matrix polynomial with degree <= ', d[3]-1,'\n', sep = '')
  return(invisible(NULL))
}

str.lmfd = function(object, ...) {
  d = attr(object, 'order')
  cat('( ',d[1],' x ',d[2],' ) left matrix fraction description with degrees (p = ',
      d[3], ', q = ', d[4],')\n', sep = '')
  return(invisible(NULL))
}

str.rmfd = function(object, ...) {
  d = attr(object, 'order')
  cat('( ',d[1],' x ',d[2],' ) right matrix fraction description with degrees (deg(c(z)) = p = ',
      d[3], ', deg(d(z)) = q = ', d[4],')\n', sep = '')
  return(invisible(NULL))
}


str.stsp = function(object, ...) {
  d = attr(object, 'order')
  cat('( ',d[1],' x ',d[2],' ) statespace realization with s = ', d[3], ' states\n', sep = '')
  return(invisible(NULL))
}

str.pseries = function(object, ...) {
  d = dim(unclass(object))
  cat('( ',d[1],' x ',d[2],' ) power series parameters with maximum lag = ', d[3]-1, '\n', sep = '')
  return(invisible(NULL))
}

str.zvalues = function(object, ...) {
  d = dim(unclass(object))
  cat('( ',d[1],' x ',d[2],' ) functional values with ', d[3], ' frequencies/points\n', sep = '')
  return(invisible(NULL))
}
# stsp_methods. R #########################################
# special methods/operations on statespace realizations


ctr_matrix = function(A, B, o = NULL) {
  if (missing(A)) stop('parameter A is missing')
  if (inherits(A, 'stsp')) {
    B = A$B
    A = A$A
    s = nrow(B)
    n = ncol(B)
  } else {
    if ( !( is.numeric(A) || is.complex(A) ) ) stop('parameter A is not numeric or complex')
    if ( (!is.matrix(A)) || (nrow(A) != ncol(A)) ) stop('parameter A is not a square matrix')
    s = nrow(A)
    if (missing(B)) stop('parameter B is missing')
    if ( !( is.numeric(B) || is.complex(B) ) ) stop('parameter B is not numeric or complex')
    if ( (!is.matrix(B)) || (nrow(B) != s) ) stop('parameters A,B are not compatible')
    n = ncol(B)
  }
  if (is.null(o)) o = s
  o = as.integer(o)[1]
  if (o < 0) stop('o must be a non negative integer')

  Cm = array(0, dim = c(s,n,o))
  Cm[,,1] = B
  for (i in iseq(2,o)) {
    B = A %*% B
    Cm[,,i] = B
  }
  dim(Cm) = c(s,o*n)
  return(Cm)
}


obs_matrix = function(A, C, o = NULL) {
  if (missing(A)) stop('parameter A is missing')
  if (inherits(A, 'stsp')) {
    C = A$C
    A = A$A
    s = ncol(C)
    m = nrow(C)
  } else {
    if ( !( is.numeric(A) || is.complex(A) ) ) stop('parameter A is not numeric or complex')
    if ( (!is.matrix(A)) || (nrow(A) != ncol(A)) ) stop('parameter A is not a square matrix')
    s = nrow(A)
    if (missing(C)) stop('parameter C is missing')
    if ( !( is.numeric(C) || is.complex(C) ) ) stop('parameter C is not numeric or complex')
    if ( (!is.matrix(C)) || (ncol(C) != s) ) stop('parameters A,C are not compatible')
    m = nrow(C)
  }
  if (is.null(o)) o = s
  o = as.integer(o)[1]
  if (o < 0) stop('o must be a non negative integer')

  Om = array(0, dim = c(s,m,o))
  A = t(A)
  C = t(C)
  Om[,,1] = C
  for (i in iseq(2,o)) {
    C = A %*% C
    Om[,,i] = C
  }
  dim(Om) = c(s,o*m)
  return(t(Om))
}


is.minimal = function(x, ...) {
  UseMethod("is.minimal", x)
}

is.minimal.stsp = function(x, tol = sqrt(.Machine$double.eps), only.answer = TRUE, ...) {
  d = dim(x)
  s = unname(d[3])
  if (prod(d) == 0) {
    H = matrix(0, nrow = d[1]*s, ncol = d[2]*s)
    svH = numeric(0)
    s0 = 0
  } else {
    ir = unclass(pseries(x, lag.max = 2*s-1))[,,-1]
    H = bhankel(ir)
    svH = svd(H, nu = 0, nv = 0)$d
    s0 = sum(svH > tol)
  }
  is_minimal = (s == s0)

  if (only.answer) return(is_minimal)

  return(list(answer = is_minimal, H = H, sv = svH, s0 = s0))
}


state_trafo = function(obj, T, inverse = FALSE) {
  if (!inherits(obj, 'stsp')) stop('argument "obj" must be "stsp" objcet')

  d = unname(dim(obj))
  m = d[1]
  n = d[2]
  s = d[3]

  if (s == 0) {
    if (length(T) == 0 ) return(obj)
    stop('The argument "T" is not compatible with "obj"')
  }

  if ( !(is.numeric(T) || is.complex(T)) ) stop('T must be a numeric (or complex) vector or a matrix!')
  if (is.vector(T)) {
    if (length(T) == s) {
      T = diag(x = T, nrow = s)
    } else {
      if (length(T) == (s^2)) {
        T = matrix(T, nrow = s, ncol = s)
      } else stop('T is not a compatible vector')
    }
  }
  if ( (!is.matrix(T)) || (ncol(T) != s) || (nrow(T) != s) ) stop('T must be a square, non-singular and compatible matrix!')

  obj = unclass(obj)

  if (inverse) {
    junk = try(solve(T, obj[1:s,,drop = FALSE]))
    if (inherits(junk, 'try-error')) stop('T is singular')
    obj[1:s,] = junk
    obj[,1:s] = obj[,1:s,drop = FALSE] %*% T
  } else {
    obj[1:s,] = T %*% obj[1:s,,drop = FALSE]
    junk = try(t(solve(t(T), t(obj[,1:s, drop = FALSE]))))
    if (inherits(junk, 'try-error')) stop('T is singular')
    obj[,1:s] = junk
  }
  obj = structure(obj, order = as.integer(c(m,n,s)),  class = c('stsp','ratm'))

  return(obj)
}

grammians = function(obj, which = c('lyapunov','minimum phase',
                                    'ctr','obs','obs_inv','ctr_inv')) {
  if (!inherits(obj, 'stsp')) stop('argument "obj" must a "stsp" object')

  which = match.arg(which)
  if (which == 'ctr') {
    P = try(lyapunov(obj$A, obj$B %*% t(obj$B), non_stable = 'stop'))
    if (inherits(P, 'try-error')) {
      stop('statespace realization is not stable')
    }
    return(P)
  }
  if (which == 'obs') {
    Q = try(lyapunov(t(obj$A), t(obj$C)%*% obj$C, non_stable = 'stop'))
    if (inherits(Q, 'try-error')) {
      stop('statespace realization is not stable')
    }
    return(Q)
  }
  if (which == 'ctr_inv') {
    d = dim(obj)
    if (d[1] != d[2]) stop('the rational matrix must be square')

    # B matrix of inverse
    B = t(solve(t(obj$D), t(obj$B)))
    if (inherits(B, 'try-error')) stop('obj is not invertible')
    # A matrix of inverse
    A = obj$A - B %*% obj$C

    P = try(lyapunov(A, B %*% t(B), non_stable = 'stop'))
    if (inherits(P, 'try-error')) {
      stop('statespace realization is not minimum phase')
    }
    return(P)
  }
  if (which == 'obs_inv') {
    d = dim(obj)
    if (d[1] != d[2]) stop('the rational matrix must be square')

    # C matrix of inverse
    C = -solve(obj$D, obj$C)
    if (inherits(C, 'try-error')) stop('obj is not invertible')
    # A matrix of inverse
    A = obj$A + obj$B %*% C

    Q = try(lyapunov(t(A), t(C) %*% C, non_stable = 'stop'))
    if (inherits(Q, 'try-error')) {
      stop('statespace realization is not minimum phase')
    }
    return(Q)
  }
  if (which == 'lyapunov') {
    P = try(lyapunov(obj$A, obj$B %*% t(obj$B), non_stable = 'stop'))
    if (inherits(P, 'try-error')) {
      stop('statespace realization is not stable')
    }
    # this is not efficient, since we compute the schur decomposition of A (above)
    # and then the schur decomposition of A' (below)
    Q = try(lyapunov(t(obj$A), t(obj$C) %*% obj$C, non_stable = 'stop'))
    if (inherits(Q, 'try-error')) {
      # this should not happen
      stop('statespace realization is not stable')
    }
    return(list(P = P, Q = Q))
  }
  if (which == 'minimum phase') {
    d = dim(obj)
    if (d[1] != d[2]) stop('the rational matrix must be square')

    P = try(lyapunov(obj$A, obj$B %*% t(obj$B), non_stable = 'stop'))
    if (inherits(P, 'try-error')) {
      stop('statespace realization is not stable')
    }
    # C matrix of inverse
    C = -solve(obj$D, obj$C)
    if (inherits(C, 'try-error')) stop('obj is not invertible')
    # A matrix of inverse
    A = obj$A + obj$B %*% C

    Q = try(lyapunov(t(A), t(C) %*% C, non_stable = 'stop'))
    if (inherits(Q, 'try-error')) {
      stop('statespace realization is not minimum phase')
    }
    return(list(P = P, Q = Q))
  }

  stop('this should not happen!')
}


balance = function(obj, gr, tol = 10*sqrt(.Machine$double.eps), s0 = NULL, truncate = TRUE) {
  # check inputs
  if (!inherits(obj,'stsp')) stop('obj must be an "stsp" object!')

  d = unname(dim(obj))
  m = d[1]
  n = d[2]
  s = d[3] # statespace dimension of the given statespace realization "obj"

  P = gr$P
  if ( (!is.numeric(P)) || (!is.matrix(P)) || (any(dim(P) != s)) ) {
    stop('argument "P" is not a compatible matrix')
  }
  Q = gr$Q
  if ( (!is.numeric(Q)) || (!is.matrix(Q)) || (any(dim(Q) != s)) ) {
    stop('argument "Q" is not a compatible matrix')
  }

  if (s == 0) {
    return(list(obj = obj, T = matrix(0, nrow = 0, ncol = 0),
                Tinv = matrix(0, nrow = 0, ncol = 0),
                P = P, Q = Q, sigma = numeric(0)))
  }

  # construct square roots of P and Q
  # is it better to use eigen()?
  # with eigen we could check that P,Q are positive semidefinite

  out = svd(P, nv = 0)
  P2 = out$u %*% diag(x = sqrt(out$d), nrow = s, ncol = s) %*% t(out$u)
  # print(all.equal(P, P2 %*% P2))

  out = svd(Q, nv = 0)
  Q2 = out$u %*% diag(x = sqrt(out$d), nrow = s, ncol = s) %*% t(out$u)
  # print(all.equal(Q, Q2 %*% Q2))

  # svd of the product P2*Q2
  out = svd(P2 %*% Q2)
  # (Hankel) singular values
  sigma = out$d

  if (!is.null(s0)) {
    s0 = as.integer(s0)[1]
    if ((s0 < 0) || (s0 > s)) stop('illegal (target) statespace dimension "s0"')
  } else {
    s0 = sum(sigma > (tol * sigma[1]))
  }

  junk = matrix(1/sqrt(sigma[1:s0]), nrow = s0, ncol = s)
  T = (t(out$v[, 1:s0, drop = FALSE]) %*% Q2) * junk
  S = (P2 %*% out$u[, 1:s0, drop = FALSE]) * t(junk)
  # print(all.equal(diag(s0), T %*% S))

  if ((truncate) || (s0 == s)) {
    # balance and truncate
    if (s0 == 0) {
      return(list(obj = stsp(D = obj$D), T = matrix(0, nrow = 0, ncol = s0),
                  Tinv = matrix(0, nrow = s0, ncol = 0),
                  P = matrix(0, nrow = 0, ncol = 0), Q = matrix(0, nrow = 0, ncol = 0),
                  sigma = sigma))
    }

    obj = stsp(T %*% obj$A %*% S, T %*% obj$B, obj$C %*% S, obj$D)
    return(list(obj = obj, T = T, Tinv = S, P = diag(x = sigma[1:s0], nrow = s0, ncol = s0),
                Q = diag(sigma[1:s0], nrow = s0, ncol = s0), sigma = sigma))
  }

  # just balance

  # extend T and S to square matrices
  out = svd(S %*% T)
  # print(out)
  T2 = t(out$u[, (s0+1):s, drop = FALSE])
  S2 = out$v[, (s0+1):s, drop = FALSE]
  # print(T %*% S2)
  # print(T2 %*% S)

  out = svd(T2 %*% S2)
  junk = matrix(1/sqrt(out$d), nrow = (s-s0), ncol = s)
  T2 = (t(out$u) %*% T2) * junk
  S2 = (S2 %*% out$v) * t(junk)
  T = rbind(T, T2)
  S = cbind(S, S2)
  # print( T %*% S)

  obj = stsp(T %*% obj$A %*% S, T %*% obj$B, obj$C %*% S, obj$D)
  Pb = diag(x = sigma, nrow = s, ncol = s)
  Qb = Pb


  # print(P)
  # print(Pb)
  # print(Pb[(s0+1):s, (s0+1:s)])
  # print(T2)
  Pb[(s0+1):s, (s0+1):s] = T2 %*% P %*% t(T2)
  Qb[(s0+1):s, (s0+1):s] = t(S2) %*% Q %*% S2

  return(list(obj = obj, T = T, Tinv = S, P = Pb, Q = Qb, sigma = sigma))
}

# Kronecker indices and echelon form ---------------------------------------------------------

#
NULL

basis2nu = function(basis, m) {
  # no basis elements => rank of H is zero
  if (length(basis)==0) {
    return(integer(m))
  }

  # What is the highest Kronecker indices?
  p = ceiling(max(basis)/m)

  # Create a matrix with one row for each variable.
  # The element in the j-th column is one if the (n*(j-1)+i)-th row of the Hankel matrix is in the basis
  in_basis = matrix(0, nrow=m, ncol=p)
  in_basis[basis] = 1

  # Calculate the Kronecker indices by summing across columns
  nu  = apply(in_basis, MARGIN=1, FUN=sum)

  # Check if there are holes in the basis, i.e. the (n*j+i)-th row is in the basis while while the (n*(j-1)+i)-th is not
  nu1 = apply(in_basis, MARGIN=1, FUN=function (x) sum(cumprod(x)))
  if (any (nu != nu1)){
    stop('This is not a nice basis, i.e. there are holes.')
  }

  return(nu)
}


nu2basis = function(nu) {
  m = length(nu) # number of variables
  p = max(nu)    # largest Kronecker index

  # all Kronecker indices are equal to zero => rank is zero
  if (p==0) return(integer(0))

  # Create as many ones (boolean) in a row as the Kronecker index indicates
  # (apply returns columns and cbinds them, therefore t())
  in_basis = t( apply(matrix(nu), MARGIN = 1, FUN = function(x) c(rep(TRUE,x),rep(FALSE,p-x)) ) )

  # Where are the ones?
  basis = which(in_basis==1)
  return(basis)
}

pseries2hankel = function(obj, Hsize = NULL) {

  # check input parameter "obj"
  k = unclass(obj)
  if ( (!(is.numeric(k) || is.complex(k))) || (!is.array(k)) || (length(dim(k)) !=3) ) {
    stop('parameter "obj" must be an "pseries" object or a 3-D array')
  }
  d = dim(k)
  m = d[1]
  n = d[2]
  lag.max = d[3] - 1
  if (lag.max < 2) {
    stop('the impulse response contains less than 2 lags')
  }

  # check size of Hankel matrix
  if (is.null(Hsize)) {
    Hsize = ceiling((lag.max+2)/2)
  }
  Hsize = as.vector(as.integer(Hsize))
  if (length(Hsize) == 1) {
    Hsize[2] = lag.max + 1 - Hsize
  }
  f = Hsize[1] # number of block rows of the Hankel matrix <=> future
  p = Hsize[2] # number of block columns of the Hankel matrix <=> past
  # the default choice is:
  # (f+p) = lag.max + 1 and f >= p+1

  if ((f < 2) || (p < 1) || (lag.max < (f+p-1)) ) {
    stop('the conditions (f>1), (p>0) and (lag.max >= f+p-1) are not satisfied')
  }

  k0 = matrix(k[,,1], nrow = m, ncol = n)

  # Hankel matrix
  H = bhankel(k[,,-1,drop=FALSE], d = c(f,p))
  attr(H,'order') = c(m,n,f,p)
  attr(H,'k0') = k0

  return(H)
}

hankel2nu = function(H, tol = sqrt(.Machine$double.eps)) {

  d = attr(H, 'order')
  m = d[1]
  n = d[2]
  f = d[3]
  p = d[4]

  if ((m*n) == 0) return(integer(m))

  # compute 'nice' basis for row space of H, via QR decomposition
  # of the transposed matrix.
  qr.H = qr(t(H), tol = tol)

  # rank of H is zero!
  if (qr.H$rank ==0) {
    nu = integer(m)
  } else {
    basis = qr.H$pivot[1:qr.H$rank]
    nu = try(basis2nu(basis, m))
    if (inherits(nu,'try-error')) {
      stop(paste(paste(basis,  collapse=' '),' is not a nice basis for the (',
                 m*f,',',n*p,') Hankel matrix', sep = ''))
    }
  }

  return(nu)
}

hankel2mu = function(H, tol = sqrt(.Machine$double.eps)) {

  d = attr(H, 'order')
  m = d[1]
  n = d[2]
  f = d[3]
  p = d[4]

  if ((m*n) == 0) return(integer(n))

  # compute 'nice' basis for column space of H, via QR decomposition
  qr.H = qr(H, tol = tol)

  if (qr.H$rank == 0) {   # rank of H is zero!
    mu = integer(n)
  } else {
    basis = qr.H$pivot[1:qr.H$rank]
    mu = try(basis2nu(basis, n)) # basis2nu does what would be needed for basis2mu, so no additional function is necessary
    if (inherits(mu,'try-error')) {
      stop(paste(paste(basis,  collapse=' '),' is not a nice basis for the (',
                 m*f,',',n*p,') Hankel matrix', sep = ''))
    }
  }

  return(mu)
}


pseries2nu = function(obj, Hsize = NULL, tol = sqrt(.Machine$double.eps)) {

  # call the helper functions pseries2hankel and hankel2nu
  nu = try(hankel2nu(pseries2hankel(obj, Hsize = Hsize), tol = tol))
  if (inherits(nu,'try-error')) {
    stop('computation of Kronecker indices failed')
  }
  return(nu)
}


# internal function
# return statespace realization template
nu2stsp_template = function(nu, D) {
  m = nrow(D)
  n = ncol(D)
  s = sum(nu)
  if ((length(nu) != m) || ( (s > 0) && (n==0) )) {
    # notw: n = 0 implies nu = (0,...,0)!
    stop('the parameters "nu" and "dim" are not compatible')
  }

  if (s == 0) {
    A = matrix(0, nrow = 0, ncol = 0)
    B = matrix(0, nrow = 0, ncol = n)
    C = matrix(0, nrow = m, ncol = 0)
    return(stsp(A, B, C, D))
  }

  basis = nu2basis(nu)
  AC = matrix(0, nrow = s + m, ncol = s)
  dependent = c(basis + m, 1:m)
  for (i in (1:length(dependent))) {
    d = abs(basis-dependent[i])
    if (min(d) == 0) {
      # dependent[i]-th row is in basis
      j = which(d == 0)
      AC[i,j] = 1
    } else {
      j = which(basis < dependent[i])
      AC[i,j] = NA_real_
    }
  }
  A = AC[1:s,,drop = FALSE]
  C = AC[(s+1):(s+m), , drop = FALSE]
  B = matrix(NA_real_, nrow = s, ncol = n)
  return(stsp(A, B, C, D))
}
# Misc Tools and Utilities ###############################################

iseq = function(from = 1, to = 1) {
  if (to<from) return(integer(0))
  return( seq(from = from, to = to) )
}


# expand_letters is an internal function
expand_letters = function(n, l = letters) {
  if (n == 0) return(character(0))
  if ((n>1) && (length(l) <= 1)) stop('"l" must have at least 2 entries!')
  l0 = l
  while (length(l) < n) {
    l = outer(l,l0,FUN = function(a,b) {paste(b,a,sep='')})
  }
  l[1:n]
}

match_vectors = function(x, y = Conj(x)) {
  x = as.vector(x)
  y = as.vector(y)
  p = length(x)
  q = length(y)
  if (q < p) stop('The length of vector "x" is larger than the length of "y"!')

  match = munkres(abs(matrix(x, nrow = p, ncol = q) -
                        matrix(y, nrow = p, ncol = q, byrow = TRUE)))
  j = integer(p)
  j[match$a[,1]] = match$a[,2]

  return(j)
}

# Array Tools ##################################################################

dbind = function(d = 1, ...) {
  d = as.integer(as.vector(d))[1]
  if (d < 1) stop('"d" must be a positive integer!')

  args = list(...)
  n.args = length(args)
  if (n.args == 0) {
    stop('no arrays to bind!')
  }

  a = args[[1]]
  if (is.vector(a)) {
    dimnames.a = names(a)
    dim(a) = length(a)
    dimnames(a) = list(dimnames.a)
  }
  dim.a = dim(a)
  dimnames.a = dimnames(a)
  if (d > length(dim.a)) {
    dim(a) = c(dim(a),rep(1, d - length(dim.a)))
    dimnames(a) = c(dimnames.a, vector('list', d - length(dim.a)))
  }
  dim.a = dim(a)
  dimnames.a = dimnames(a)

  if (n.args == 1) {
    return(a)
  }

  b = args[[2]]
  if (is.vector(b)) {
    dimnames.b = names(b)
    dim(b) = length(b)
    dimnames(b) = list(dimnames.b)
  }
  dim.b = dim(b)
  dimnames.b = dimnames(b)
  if (d > length(dim.b)) {
    dim(b) = c(dim(b),rep(1, d - length(dim.b)))
    dimnames(b) = c(dimnames.b, vector('list', d - length(dim.b)))
  }
  dim.b = dim(b)
  dimnames.b = dimnames(b)

  n.dim = length(dim.a)
  if ( (length(dim.b) != n.dim) || ( (n.dim > 1) && (any(dim.a[-d]!=dim.b[-d]))) ) {
    stop('arrays are not compatible!')
  }

  # bind arrays a and b -> c
  dim.c = dim.a # dimension of the result
  dim.c[d] = dim.a[d] + dim.b[d]
  p = c((1:n.dim)[-d],d)
  pp = (1:n.dim)
  pp[p] = 1:n.dim
  a = aperm(a, p)
  b = aperm(b, p)
  dim.c = dim.c[p]
  c = array(c(a,b), dim = dim.c)
  c = aperm(c, pp)

  # take care of dimnames
  if ( (is.null(dimnames.a)) && (!is.null(dimnames.b)) ) {
    dimnames.a = dimnames.b
    dimnames.a[[d]] = character(0)
  }
  if ( (is.null(dimnames.b)) && (!is.null(dimnames.a)) ) {
    dimnames.b = dimnames.a
    dimnames.b[[d]] = character(0)
  }

  if ( !is.null(dimnames.a) ) {
    dimnames.c = dimnames.a
    if (is.null(names(dimnames.c))) names(dimnames.c) = names(dimnames.b)

    for (i in (1:n.dim)) {
      if ( (names(dimnames.c)[i] == '') && (!is.null(names(dimnames.b))) ) {
        names(dimnames.c)[i] = names(dimnames.b)[i]
      }

      if (i != d) {
        if ( (length(dimnames.c[[i]]) == 0) && (length(dimnames.b[[i]]) > 0) ) {
          dimnames.c[[i]] = dimnames.b[[i]]
        }
      } else {
        if ( (length(dimnames.a[[i]]) == 0) || (length(dimnames.b[[i]]) == 0) ) {
          dimnames.c[[i]] = character(0)
        } else {
          dimnames.c[[i]] = c(dimnames.a[[i]], dimnames.b[[i]])
        }
      }

    }
    dimnames(c) = dimnames.c
  }

  if (n.args == 2) {
    return(c)
  }

  args[[2]] = c
  cc = do.call(dbind, args = c(list(d = d), args[-1]))
  return(cc)
}




array2data.frame = function(x, rows = NULL, cols = NULL) {

  # check parameter "x"
  if ( (!is.array(x)) || (min(dim(x)) <= 0) || (length(dim(x)) == 0) ) {
    stop('"x" must be a (non empty) array!')
  }
  dim.x = dim(x)
  n.dims = length(dim.x)
  dims = (1:n.dims)
  dimnames.x = dimnames(x)
  if (is.null(dimnames.x)) stop('"x" must have a complete "dimnames" attribute!')
  names.dimnames.x = names(dimnames.x)
  if ( (is.null(names.dimnames.x)) || (any(names.dimnames.x =='')) ||
       (length(unique(names.dimnames.x)) < n.dims ) ) {
    stop('"x" must have a complete "dimnames" attribute!')
  }
  for (i in dims) {
    if ( (is.null(dimnames.x[[i]])) || (any(dimnames.x[[i]] =='')) ||
         (length(unique(dimnames.x[[i]])) < dim.x[i]) ) {
      stop('"x" must have a complete "dimnames" attribute!')
    }
  }

  # check input parameters "rows" and "cols"
  if (is.null(rows) && is.null(cols)) {
    stop('missing parameters "rows" and "cols"!')
  }

  # check "rows"
  if (!is.null(rows)) {
    rows = as.integer(as.vector(rows))
    # rows must be subset of dim.x
    if ( (length(rows) > 0) && ( (min(rows) < 1) || (max(rows) > n.dims) ) ) {
      stop('parameter "rows" does not correspond to a selection of dimensions (of "x")!')
    }
  }

  # check "cols"
  if (!is.null(cols)) {
    cols = as.integer(as.vector(cols))
    # cols must be subset of dim.x
    if ( (length(cols) > 0) && ( (min(cols) < 1) || (max(cols) > n.dims) ) ) {
      stop('parameter "cols" does not correspond to a selection of dimensions (of "x")!')
    }
  }

  if (is.null(rows)) {
    if (length(cols) > 0) {
      rows = dims[-cols]
    } else {
      rows = dims
    }
  }
  if (is.null(cols)) {
    if (length(rows) > 0) {
      cols = dims[-rows]
    } else {
      cols = dims
    }
  }

  # cat('rows:',rows,'cols:',cols,'\n')

  # check that the union of "rows" and "cols" = "dim.x"
  if (!isTRUE(all.equal(dims, sort(c(rows,cols))))) {
    stop('the parameters "rows" and "cols" do not correspond to a partition of the set of dimensions (of "x")!')
  }

  # coerce x to a matrix
  x = aperm(x, c(rows,cols))
  dim(x) = c(prod(dim.x[rows]),prod(dim.x[cols]))

  if (length(rows) > 0) {
    cases = expand.grid(dimnames.x[rows], KEEP.OUT.ATTRS = FALSE, stringsAsFactors = TRUE)
    # print(cases)
    # print(x)
    x = cbind(cases, x)
    # print(x)
    # print(str(x))
  } else {
    x = data.frame(x)
  }

  if (length(cols) > 0) {
    variables = expand.grid(dimnames.x[cols], KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    # print(variables)
    # print(str(variables))
    variables = apply(variables, MARGIN = 1, FUN = paste, collapse='.')
    # print(variables)
  } else {
    variables = 'value'
  }

  if (length(rows)>0) {
    colnames(x) = c(names.dimnames.x[rows],variables)
  } else {
    colnames(x) = variables
  }

  return(x)
}


# Block Matrices #################################################################

bdiag = function(...) {

  args = list(...)
  n_args = length(args)

  # no inputs -> return NULL
  if (n_args == 0) return(NULL)

  # skip 'NULL' arguments
  i = sapply(args, FUN = function(x) (!is.null(x)))

  if (sum(i) == 0) return(NULL)
  args = args[i]
  n_args = length(args)

  # check type of arguments
  i = sapply(args, FUN = function(x) ( !( is.logical(x) || is.integer(x) || is.numeric(x) || is.complex(x) ) ) )
  if (any(i)) stop('inputs must be of type logical, integer, numeric or complex!')
  i = sapply(args, FUN = function(x) ( !( is.vector(x) || is.matrix(x) ) ) )
  if (any(i)) stop('inputs must be vectors or matrices!')

  A = args[[1]]
  if (is.vector(A)) {
    rownamesA = names(A)
    A = diag(A, nrow = length(A), ncol = length(A))
    if (!is.null(rownamesA)) {
      rownames(A) = rownamesA
      colnames(A) = rownamesA
    }
  }
  rownamesA = rownames(A)
  colnamesA = colnames(A)

  if (n_args == 1) return( A )

  for (k in (2:n_args)) {
    B = args[[k]]

    if (is.vector(B)) {
      rownamesB = names(B)
      B = diag(B, nrow = length(B), ncol = length(B))
      if (!is.null(rownamesB)) {
        rownames(B) = rownamesB
        colnames(B) = rownamesB
      }
    }
    rownamesB = rownames(B)
    colnamesB = colnames(B)

    if ( !((is.null(rownamesA)) && (is.null(rownamesB))) ) {
      if (is.null(rownamesA)) rownamesA = character(nrow(A))
      if (is.null(rownamesB)) rownamesB = character(nrow(B))
      rownamesA = c(rownamesA, rownamesB)
    }
    if ( !((is.null(colnamesA)) && (is.null(colnamesB))) ) {
      if (is.null(colnamesA)) colnamesA = character(ncol(A))
      if (is.null(colnamesB)) colnamesB = character(ncol(B))
      colnamesA = c(colnamesA, colnamesB)
    }
    A = rbind( cbind(A, matrix(FALSE, nrow = nrow(A), ncol = ncol(B))),
               cbind(matrix(FALSE, nrow = nrow(B), ncol = ncol(A)), B) )

  }

  rownames(A) = rownamesA
  colnames(A) = colnamesA
  return(A)
}

bmatrix = function(x, rows = NULL, cols = NULL) {
  # check input parameters
  if (is.null(rows) && is.null(cols)) {
    stop('missing parameters "rows" and "cols"!')
  }
  if (is.vector(x)) {
    x = array(x, dim = length(x))
  }
  if (!is.array(x)) {
    stop('"x" must be a vector, matrix or array!')
  }
  dim.x = dim(x)
  n.dims = length(dim.x)
  dims = (1:n.dims)

  # check "rows"
  if (!is.null(rows)) {
    rows = as.integer(as.vector(rows))
    # rows must be subset of dims
    if ( (length(rows) > 0) && ( (min(rows) < 1) || (max(rows) > n.dims) ) ) {
      stop('parameter "rows" does not correspond to a selection of dimensions (of "x")!')
    }
  }

  # check "cols"
  if (!is.null(cols)) {
    cols = as.integer(as.vector(cols))
    # cols must be subset of dims
    if ( (length(cols) > 0) && ( (min(cols) < 1) || (max(cols) > n.dims) ) ) {
      stop('parameter "cols" does not correspond to a selection of dimensions (of "x")!')
    }
  }

  if (is.null(rows)) {
    if (length(cols) > 0) {
      rows = dims[-cols]
    } else {
      rows = dims
    }
  }
  if (is.null(cols)) {
    if (length(rows) > 0) {
      cols = dims[-rows]
    } else {
      cols = dims
    }
  }

  #  cat('rows:',rows,'cols:',cols,'\n')

  # check that the union of "rows" and "cols" = "dims"
  if (!isTRUE(all.equal(dims, sort(c(rows,cols))))) {
    stop('the parameters "rows" and "cols" do not correspond to a partition of the set of dimensions (of "x")!')
  }

  # coerce x to a matrix
  x = aperm(x, c(rows,cols))
  dim(x) = c(prod(dim.x[rows]),prod(dim.x[cols]))

  return(x)
}


btoeplitz = function(R, C) {

  # Check for correct inputs: Vector, matrix, or array ####
  # Both missing?
  if ((missing(R)) && (missing(C))) {
    stop('At least one parameter "R" or "C" has to be supplied.')
  }

  # R correct?
  if (!missing(R)) {
    if (is.vector(R)) {
      R = array(R, dim = c(1, 1, length(R)))
    }
    if (is.matrix(R)) {
      dim(R) = c(nrow(R), ncol(R), 1)
    }
    if ((!is.array(R)) || (length(dim(R)) != 3)) {
      stop('R must be a vector, a matrix or a 3-dim array!')
    }
    dimR = dim(R)
  }

  # C correct?
  if (!missing(C)) {
    if (is.vector(C)) {
      C = array(C, dim = c(1, 1, length(C)))
    }
    if (is.matrix(C)) {
      dim(C) = c(nrow(C), ncol(C), 1)
    }
    if ((!is.array(C)) || (length(dim(C)) != 3)) {
      stop('C must be a vector, a matrix or a 3-dim array!')
    }
    dimC = dim(C)
  }

  # If one argument is missing, replace it such that we obtain a symmetric Toeplitz matrix ####

  if (missing(R)) {
    if (dim(C)[1] != dim(C)[2]) stop('if "R" is missing then dim(C)[1]=dim(C)[2] must hold!')
    # create a 'symmetric' block toeplitz matrix
    R = aperm(C, c(2, 1, 3))
    dimR = dim(R)
  }

  if (missing(C)) {
    if (dim(R)[1] != dim(R)[2]) stop('if "C" is missing then dim(R)[1]=dim(R)[2] must hold!')
    # create a 'symmetric' block toeplitz matrix
    C = aperm(R, c(2, 1, 3))
    dimC = dim(C)
  }

  # Check for compatible dimensions ####
  if (any(dimR[1:2] != dimC[1:2])) {
    stop('Dimensions of "R" and "C" are not compatible.')
  }

  # Define integer valued parameters (dimensions) ####
  n_rows = dimR[1] # p
  n_cols = dimR[2] # q
  n_depth_R = dimR[3] # n
  n_depth_C = dimC[3] # m

  # Deal with trivial cases (one array consists of 0 or 1 matrix, i.e. there is no third dimension) ####
  if (n_depth_C == 0)
    return(matrix(0, nrow = 0, ncol = n_cols * n_depth_R))
  if (n_depth_C == 1)
    return(matrix(R, nrow = n_rows, ncol = n_cols * n_depth_R))

  if (n_depth_R == 0)
    return(matrix(0, nrow = n_rows * n_depth_C, ncol = 0))
  if (n_depth_R == 1) {
    C[, , 1] = R[, , 1]
    return(matrix(aperm(C, c(1, 3, 2)), nrow = n_rows * n_depth_C, n_cols))
  }

  # Concatenate C and R, and put them in the correct order (to fill block Toeplitz matrix in column-major order) ####
  CR = array(0, dim = c(n_rows, n_cols, n_depth_C + n_depth_R - 1))
  CR[, , 1:(n_depth_C - 1)] = C[, , n_depth_C:2] # order needs to be inverted because we fill from top to bottom
  CR[, , n_depth_C:(n_depth_C + n_depth_R - 1)] = R

  # Create indices
  mat_tmp = matrix(raw(), nrow = n_depth_C, ncol = n_depth_R)
  idx_mat = col(mat_tmp) - row(mat_tmp) + n_depth_C

  # Create Block Toeplitz matrix
  T = CR[, , idx_mat, drop = FALSE]
  dim(T) = c(n_rows, n_cols, n_depth_C, n_depth_R)
  T = aperm(T, c(1, 3, 2, 4))
  dim(T) = c(n_depth_C * n_rows, n_depth_R * n_cols)
  return(T)
}


bhankel = function(R, d = NULL) {
  if (is.vector(R)) R = array(R,dim=c(1,1,length(R)))
  if (is.matrix(R)) R = array(R,dim=c(nrow(R),ncol(R),1))
  dim = dim(R)
  if (length(dim)!=3) stop('"R" must be a vector, matrix or 3-dimensional array!')
  p = dim[1]
  q = dim[2]
  k = dim[3]
  # if (k==0) stop('"R" is an empty array!')

  if ( is.null(d) ) {
    d = ceiling((k+1)/2)
  }
  d = as.integer(as.vector(d))
  if (min(d)<0) stop('"d" has negative entries!')
  if (length(d)==1) d = c(d,max(k+1-d,1))
  d = d[1:2]
  m = d[1]
  n = d[2]

  if (min(c(m,n,p,q))==0) return(matrix(0,nrow = m*p,ncol=n*q))

  # pad with zeros
  if (k < (m+n-1)) R = array(c(R,double((m+n-1-k)*p*q)),dim=c(p,q,m+n-1))

  if (m==1) return(matrix(R[,,1:n],nrow=p,ncol=q*n))

  if (n==1) {
    return(matrix(aperm(R[,,1:m,drop=FALSE],c(1,3,2)),nrow=m*p,q))
  }

  # j+i-1
  ji = matrix(1:n,nrow=m,ncol=n,byrow = TRUE) + matrix(0:(m-1),nrow=m,ncol=n)
  T = array(0,dim=c(p,q,m*n))
  T = R[,,ji,drop=FALSE]
  dim(T) = c(p,q,m,n)
  T = aperm(T,c(1,3,2,4))
  dim(T) = c(m*p,n*q)
  return(T)
}



# Linear Indices vs Matrix indices ##############################################

ind2sub = function(dim, ind){
  row = ((ind-1) %% dim[1]) + 1
  col = floor((ind-1) / dim[1]) + 1
  return(c(row, col))
}


sub2ind = function(dim, row, col){
  ind = (col-1)*dim[1] + row
  return(ind)
}

# QL and LQ decomposition ####
ql_decomposition = function(x,...) {
  # Check inputs
  if (length(x) == 0) {
    stop('x is an empty matrix!')
  }
  # coerce vectors to matrices
  x = as.matrix(x)
  m = nrow(x)
  n = ncol(x)

  x = x[, n:1, drop = FALSE]

  qr_x = qr(x,...)
  if (any(qr_x$pivot != 1:n)){
    stop("QL implementation via QR does not work if QR decomposition pivots columns.")
  }

  r = qr.R(qr_x)
  q = qr.Q(qr_x)

  # take care of signs of diagonal elements of r (> 0)
  i = (diag(r) < 0)
  if (any(i)) {
    q[, i] = -q[, i]
    r[i, ] = -r[i, ]
  }

  return(list(l = r[nrow(r):1, ncol(r):1, drop = FALSE],
              q = q[, ncol(q):1, drop = FALSE]))
}



lq_decomposition = function(x, ...){

  # Check inputs
  if (length(x) == 0) {
    stop('x is an empty matrix!')
  }

  # QR of transpose of x (in order to obtain the LQ of x)
  tx = t(x) # note t() converts vectors to matrices!
  qr_tx = qr(tx, ...)

  if (any(qr_tx$pivot != 1:ncol(tx))){
    stop("This function does not work if pivoting is used by the base::qr() function.
         Only way to resolve this error is using a version without pivoting
         (of which BF is not aware to exist in R; but it does so in scipy library)")
  }
  r = qr.R(qr_tx)
  q = qr.Q(qr_tx)

  # take care of signs of diagonal elements of r (> 0)
  i = (diag(r) < 0)
  if (any(i)) {
    q[, i] = -q[, i]
    r[i, ] = -r[i, ]
  }

  return(
    list(l = t(r),
         q = t(q))
  )
}

# test_ objects (lmfd, rmfd, polm, stsp) ####


test_array = function(dim, random = FALSE, dimnames = FALSE) {
  # check input parameter "dim"
  dim = as.vector(as.integer(dim))
  if ((length(dim)==0) || (min(dim) < 0)) {
    stop('"dim" must be a (non empty) vector of non negative integers!')
  }
  n.dims = length(dim)

  if (min(dim)==0) {
    # empty array
    x = array(0, dim=dim)
  } else {
    if (random) {
      x = array(stats::rnorm(prod(dim)), dim = dim)
    } else {
      x = (1:dim[1])
      for (i in (iseq(2,n.dims))) {
        x = outer(10*x,(1:dim[i]),FUN = '+')
      }
      x = array(x, dim = dim)
    }
  }

  # set dimnames
  if (dimnames) {
    dimnames.x = as.list(1:n.dims)
    names(dimnames.x) = expand_letters(n.dims, l = LETTERS)
    for (i in (1:n.dims)) {
      if (dim[i]>0) {
        dimnames.x[[i]] = paste(names(dimnames.x)[i],1:dim[i],sep = '=')
      } else {
        dimnames.x[[i]] = character(0)
      }
    }
    dimnames(x) = dimnames.x
  }

  return(x)
}


test_rmfd = function(dim = c(1,1), degrees = c(1,1), digits = NULL,
                     bpoles = NULL, bzeroes = NULL, n.trials = 100) {
  # check input parameter "dim"
  dim = as.integer(dim) # note: as.integer converts to vector!
  if ((length(dim) != 2) || (dim[1] <= 0) || (dim[2] < 0)) {
    stop('argument "dim" must be a vector of integers with length 2, dim[1] > 0 and dim[2] >= 0!')
  }
  # check input parameter "degree"
  degrees = as.integer(degrees) # note: as.integer converts to vector!
  if ((length(degrees) != 2) || (degrees[1] < 0) || (degrees[2] < (-1))) {
    stop('argument "degrees" must be a vector of integers with length 2, degrees[1] >= 0 and degrees[2] >= -1!')
  }

  m = dim[1]
  n = dim[2]
  p = degrees[1]
  q = degrees[2]

  if (p == 0) bpoles = NULL
  if ( (m != n) || (q == 0) ) bzeroes = NULL

  i.trial = 0
  err = TRUE
  sd = 1
  while ( (i.trial < n.trials) && (err) ) {
    c = cbind(diag(n), matrix(stats::rnorm(n*n*p, sd = sd), nrow = n, ncol = n*p))
    dim(c) = c(n,n,p+1)
    d = matrix(stats::rnorm(m*n*(q+1), sd = sd), nrow = m*(q+1), ncol = n)
    dim(d) = c(m,n,q+1)
    if (!is.null(digits)) {
      c = round(c, digits)
      d = round(d, digits)
    }
    x = rmfd(c,d)

    err = FALSE
    if ( !is.null(bpoles) ) {
      err = try(min(abs(poles(x, print_message = FALSE))) <= bpoles, silent = TRUE)
      if (inherits(err, 'try-error')) err = TRUE
    }
    if ( (!err) && (!is.null(bzeroes)) ) {
      err = try((min(abs(zeroes(x, print_message = FALSE))) <= bzeroes), silent = TRUE)
      if (inherits(err, 'try-error')) err = TRUE
    }
    i.trial = i.trial + 1
    sd = sd/1.1
  }
  if (err) {
    stop('Could not generate a suitable rational matrix with ', n.trials, ' trials!')
  }
  return(x)
}


test_lmfd = function(dim = c(1,1), degrees = c(1,1), digits = NULL,
                     bpoles = NULL, bzeroes = NULL, n.trials = 100) {
  # check input parameter "dim"
  dim = as.integer(dim) # note: as.integer converts to vector!
  if ((length(dim) != 2) || (dim[1] <= 0) || (dim[2] < 0)) {
    stop('argument "dim" must be a vector of integers with length 2, dim[1] > 0 and dim[2] >= 0!')
  }
  # check input parameter "degree"
  degrees = as.integer(degrees) # note: as.integer converts to vector!
  if ((length(degrees) != 2) || (degrees[1] < 0) || (degrees[2] < (-1))) {
    stop('argument "degrees" must be a vector of integers with length 2, degrees[1] >= 0 and degrees[2] >= -1!')
  }

  m = dim[1]
  n = dim[2]
  p = degrees[1]
  q = degrees[2]

  if (p == 0) bpoles = NULL
  if ( (m != n) || (q == 0) ) bzeroes = NULL

  i.trial = 0
  err = TRUE
  sd = 1
  while ( (i.trial < n.trials) && (err) ) {
    a = cbind(diag(m), matrix(stats::rnorm(m*m*p, sd = sd), nrow = m, ncol = m*p))
    dim(a) = c(m,m,p+1)
    b = matrix(stats::rnorm(m*n*(q+1), sd = sd), nrow = m, ncol = n*(q+1))
    dim(b) = c(m,n,q+1)
    if (!is.null(digits)) {
      a = round(a, digits)
      b = round(b, digits)
    }
    x = lmfd(a,b)

    err = FALSE
    if ( !is.null(bpoles) ) {
      err = try(min(abs(poles(x, print_message = FALSE))) <= bpoles, silent = TRUE)
      if (inherits(err, 'try-error')) err = TRUE
    }
    if ( (!err) && (!is.null(bzeroes)) ) {
      err = try((min(abs(zeroes(x, print_message = FALSE))) <= bzeroes), silent = TRUE)
      if (inherits(err, 'try-error')) err = TRUE
    }
    i.trial = i.trial + 1
    sd = sd/1.1
  }
  if (err) {
    stop('Could not generate a suitable rational matrix with ', n.trials, ' trials!')
  }
  return(x)
}

# only real matrices
# col_end_matrix changes the degree
# col_end_matrix for degree -1 columns is ignored

test_polm = function(dim = c(1,1), degree = 0, random = FALSE, digits = NULL, col_end_matrix = NULL,
                     value_at_0 = NULL, bzeroes = NULL, n.trials = 100) {

  # Check inputs: "dim" ####
  dim = as.integer(dim) # note: as.integer converts to vector!
  if ((length(dim) != 2) || (min(dim) < 0)) {
    stop('argument "dim" must be a vector of non negative integers with length 2!')
  }

  # Empty matrix: ignore all other parameters ####
  if (prod(dim) == 0) {
    return(polm(array(0, dim = c(dim,0))))
  }

  # Check inputs: "degree" ####
  # degree can be (in addition to an integer) a vector (prescribing column degrees) or a matrix
  if (is.vector(degree)) {
    degree = as.integer(degree)
    if (length(degree) == 1) degree = rep(degree, dim[2])
    if (length(degree) != dim[2]) stop('parameters "degree" and "dim" are not compatible.')
    degree = matrix(degree, nrow = dim[1], ncol = dim[2], byrow = TRUE)
  }
  if (is.matrix(degree)) {
    # coerce to integer matrix. note that as.integer() returns a vector!
    degree = matrix(as.integer(degree), nrow = nrow(degree), ncol = ncol(degree))
  }
  if ( (!is.integer(degree)) || (!is.matrix(degree)) ||
       any(dim(degree) != dim) || any(is.na(degree)) || any(degree < -1) ) {
    stop('argument "degree" must be a scalar, vector or matrix of integers (>= -1), compatible to "dim".')
  }

  p = max(degree)
  if (p == (-1)) {
    # zero polynomial
    return(polm(array(0, dim = c(dim,0))))
  }
  # compute column degrees
  col_degree = apply(degree, MARGIN = 2, FUN = max)

  # Create polm: Fixed entries ####
  if (!random) {
    x = test_array(dim = c(dim, p+1)) - 1
    # impose the desired degrees!
    for (i in (1:dim[1])) {
      for (j in (1:dim[2])) {
        if (degree[i,j] < p) {
          x[i,j,(degree[i,j]+2):(p+1)] = 0
        }
      }
    }
    return(polm(x))
  }

  # Create polm: Random entries ####

  # create an array with NA's for the "free" parameters
  x0 = array(NA_real_, dim = c(dim, p+1))

  # impose the desired degrees!
  for (i in (1:dim[1])) {
    for (j in (1:dim[2])) {
      if (degree[i,j] < p) {
        x0[i, j, (degree[i,j]+2):(p+1)] = 0
      }
    }
  }

  # impose column end matrix
  if (!is.null(col_end_matrix)) {
    # check input parameter col_end_matrix
    if ( !is.numeric(col_end_matrix) || !is.matrix(col_end_matrix) ||
         any(dim(col_end_matrix) != dim) ) {
      stop('argument "col_end_matrix"  is not compatible')
    }
    if (min(col_degree) < 0) {
      stop('some column degrees are negative, but column end matrix has been given')
    }
    for (i in (1:dim[2])) {
      x0[,i,col_degree[i]+1] = col_end_matrix[,i]
    }
  }

  # impose value at z=0
  if (!is.null(value_at_0)) {
    # check input parameter col_end_matrix
    if ( !is.numeric(value_at_0) || !is.matrix(value_at_0) ||
         any(dim(value_at_0) != dim) ) {
      stop('argument "value_at_0"  is not compatible')
    }
    if (min(col_degree) < 0) {
      stop('some column degrees are negative, but value at z=0 has been given')
    }
    x0[,,1] = value_at_0
  }

  i = is.na(x0)
  n_theta = sum(i) # number of "free" parameters

  # for non-square polynomials, or polynomials with degree p=0, ignore the parameter "bzeroes"
  if ( (dim[1] != dim[2]) || (p <= 0) ) bzeroes = NULL

  if (n_theta == 0) {
    x = polm(x0)
    if (!is.null(bzeroes)) {
      err = try(min(abs(zeroes(x, print_message = FALSE))) <= bzeroes, silent = TRUE)
      if ( inherits(err, 'try-error') ) err = TRUE
      if (err) {
        stop('the zeroes of the generated polynomial do not satisfy the constraint that their absolute values are larger than (bzeroes)!')
      }
    }
    return(x)
  }

  i.trial = 0
  err = TRUE
  sd = 1
  while ( (i.trial < n.trials) && (err) ) {
    theta = stats::rnorm(n_theta, sd = sd)
    if (!is.null(digits)) theta = round(theta, digits)
    x0[i] = theta
    x = polm(x0)
    err = FALSE

    if ( !is.null(bzeroes) ) {
      err = try(min(abs(zeroes(x, print_message = FALSE))) <= bzeroes, silent = TRUE)
      if ( inherits(err, 'try-error') ) err = TRUE
    }
    sd = sd/1.1
    i.trial = i.trial + 1
  }

  if (err) {
    stop('Could not generate a suitable rational matrix with ', n.trials, ' trials!')
  }
  return(x)
}



test_lpolm = function(dim = c(1,1), degree_max = 1, degree_min = -1,
                      random = FALSE,
                      col_start_matrix = NULL, value_at_0 = NULL, col_end_matrix = NULL){

  # Check inputs: "dim" ####
  dim = as.integer(dim) # note: as.integer converts to vector!
  if ((length(dim) != 2) || (min(dim) < 0)) {
    stop('argument "dim" must be a vector of non negative integers with length 2!')
  }

  # Empty matrix: ignore all other parameters ####
  if (prod(dim) == 0) {
    return(lpolm(array(0, dim = c(dim,0)), min_deg = 0))
  }

  # Check inputs: degree_max and degree_min ####
  # Degree can be (in addition to an integer) a vector (prescribing column degrees) or a matrix
  # Eventually, it will be transformed to a matrix (scalar or vector fill up matrix)

  # __ degree_max ####
  if (is.vector(degree_max)) {
    degree_max = as.integer(degree_max)
    if (length(degree_max) == 1) degree_max = rep(degree_max, dim[2])
    if (length(degree_max) != dim[2]) stop('Parameters "degree_max" and "dim" are not compatible. If "degree_max" is a vector it must correspond to the (maximal) column degrees.')
    degree_max = matrix(degree_max, nrow = dim[1], ncol = dim[2], byrow = TRUE)
  }
  if (is.matrix(degree_max)) {
    # coerce to integer matrix. note that as.integer() returns a vector!
    degree_max = matrix(as.integer(degree_max), nrow = nrow(degree_max), ncol = ncol(degree_max))
  }
  if ( (!is.integer(degree_max)) ||
       (!is.matrix(degree_max)) ||
       any(dim(degree_max) != dim) ||
       any(is.na(degree_max)) ||
       any(degree_max < -1) ) {
    stop('argument "degree_max" must be a scalar, vector or matrix of integers (>= -1), compatible to "dim".')
  }

  # compute column degrees and overall max degree
  col_degree_max = apply(degree_max, MARGIN = 2, FUN = max)
  p = max(degree_max)

  # __ degree_min ####
  if (is.vector(degree_min)) {
    degree_min = as.integer(degree_min)
    if (length(degree_min) == 1) degree_min = rep(degree_min, dim[2])
    if (length(degree_min) != dim[2]) stop('Parameters "degree_min" and "dim" are not compatible. If "degree_min" is a vector, it corresponds to minimal column degrees.')
    degree_min = matrix(degree_min, nrow = dim[1], ncol = dim[2], byrow = TRUE)
  }
  if (is.matrix(degree_min)) {
    # coerce to integer matrix. note that as.integer() returns a vector!
    degree_min = matrix(as.integer(degree_min), nrow = nrow(degree_min), ncol = ncol(degree_min))
  }
  if ( (!is.integer(degree_min)) ||
       (!is.matrix(degree_min)) ||
       any(dim(degree_min) != dim) ||
       any(is.na(degree_min)) ||
       any(degree_min > 0) ) {
    stop('Argument "degree_min" must be a scalar, vector or matrix of integers (<= 0), compatible with "dim". \n The constant term must be given through the polynomial part.')
  }
  # Handling of zero polynomial: This would correspond to degree_min = 0 and degree_max = -1, i.e. degree_min > degree_max.
  # This case should be rather handled with degree_min = 0 and degree_max = 0, and setting value_at_0 such that the element is indeed zero
  if (any(degree_min > degree_max)){
    stop("degree_min is not allowed to be larger than degree_max. Handle the zero polynomial by setting value_at_0 or by using test_polm().")
  }

  # compute column degrees and min overall degree
  col_degree_min = apply(degree_min, MARGIN = 2, FUN = min)
  q = min(degree_min)

  # Construct Laurent polynomial ####
  if ((p == -1) && (q == 0)) {
    # zero polynomial (this should not happen)
    return(lpolm(array(0, dim = c(dim,0)), min_deg = 0))
  }

  # [lp(z)]_{+} part ####
  aa_p = test_array(dim = c(dim[1], dim[2], p+1)) - 1

  # impose the desired degrees!
  for (i in (1:dim[1])) {
    for (j in (1:dim[2])) {
      if (degree_max[i,j] < p) {
        aa_p[i,j,(degree_max[i,j]+2):(p+1)] = 0 # here is the advantage for defining deg(zero poly) = -1...
      }
    }
  }

  # [lp(z)]_{-} part ####
  aa_m = -test_array(dim = c(dim[1], dim[2], -q))

  # impose the desired degrees!
  for (i in (1:dim[1])) {
    for (j in (1:dim[2])) {
      if (-degree_min[i,j] < -q) {
        aa_m[i,j,(-degree_min[i,j]+1):(-q)] = 0 # here is the advantage of defining deg(zero poly) = -1...
      }
    }
  }

  if (q != 0){
    aa_m = aa_m[,,(-q):1, drop = FALSE]
  }

  # Concatenate negative and non-negative partpositive
  aa = dbind(d = 3, aa_m, aa_p)

  # Fill random coefficients ####
  n_par = sum(aa != 0)
  aa[aa != 0] = stats::rnorm(n_par, sd = 1)

  # Set column_start_matrix, value_at_0, column_end_matrix (in this order) ####
  # __Column start matrix ####
  if (!is.null(col_start_matrix)) {
    # check input parameter col_start_matrix
    if ( !is.numeric(col_start_matrix) || !is.matrix(col_start_matrix) ||
         any(dim(col_start_matrix) != dim) ) {
      stop('argument "col_start_matrix"  is not compatible')
    }
    if (min(col_degree_min) > -1) {
      stop('some column degrees are non-negative, but column start matrix has been given')
    }
    for (ix_col in (1:dim[2])) {
      aa[,ix_col,q+1+col_degree_min[ix_col]] = col_start_matrix[,ix_col]
    }
  }

  # __Value at zero ####
  if (!is.null(value_at_0)) {
    # check input parameter col_end_matrix
    if ( !is.numeric(value_at_0) || !is.matrix(value_at_0) ||
         any(dim(value_at_0) != dim) ) {
      stop('argument "value_at_0"  is not compatible')
    }
    if (min(col_degree_max) < 0) {
      stop('some column degrees are negative, but value at z=0 has been given')
    }
    aa[,,q+1] = value_at_0
  }

  # __Column end matrix ####
  if (!is.null(col_end_matrix)) {
    # check input parameter col_end_matrix
    if ( !is.numeric(col_end_matrix) || !is.matrix(col_end_matrix) ||
         any(dim(col_end_matrix) != dim) ) {
      stop('argument "col_end_matrix"  is not compatible')
    }
    if (min(col_degree_max) < 0) {
      stop('some column degrees are negative, but column end matrix has been given')
    }
    for (ix_col in (1:dim[2])) {
      aa[,i,q+1+col_degree_max[ix_col]] = col_end_matrix[,ix_col]
    }
  }

  return(lpolm(aa, min_deg = q))
}


test_stsp = function(dim = c(1,1), s = NULL, nu = NULL, D = NULL,
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

  # check input parameter "nu"
  if (!is.null(nu)) {
    nu = as.integer(nu) # as.integer converts to vector
    if ((length(nu) != m) || (min(nu) < 0)) {
      stop('"nu" must be a vector of non negative integers with length equal to dim[1]!')
    }
    tmpl = unclass(nu2stsp_template(nu, D))
    s = sum(nu)
  } else {
    if (is.null(s)) stop('either "s" or "nu" must be specified')
    s = as.integer(s)[1]
    if (s < 0) stop('parameter "s" must be a nonnegative integer')
    tmpl = stsp(A = matrix(NA_real_, nrow = s, ncol = s),
                B = matrix(NA_real_, nrow = s, ncol = n),
                C = matrix(NA_real_, nrow = m, ncol = s),
                D = D)
    tmpl = unclass(tmpl)
  }
  order = c(m,n,s)

  if (m != n) bzeroes = NULL # ignore bound for zeroes for non-square matrices
  if ( (m*n == 0) || (s == 0) ) {
    # ignore bounds for poles and zeroes for empty matrices or state dimension = 0
    bpoles = NULL
    bzeroes = NULL
  }

  i = which(is.na(tmpl))
  if (length(i) == 0) {
    x = structure(tmpl, class = c('stsp','ratm'), order = order)
    # check poles and zeroes
    if (!is.null(bpoles)) {
      err = try((min(abs(poles(x, print_message = FALSE))) <= bpoles), silent = TRUE)
      if ( inherits(err, 'try-error') || err ) {
        stop('the poles of the generated rational matrix do not satisfy the prescribed bound.')
      }
    }
    if (!is.null(bzeroes)) {
      err = try((min(abs(zeroes(x, print_message = FALSE))) <= bzeroes), silent = TRUE)
      if ( inherits(err, 'try-error') || err ) {
        stop('the zeroes of the generated rational matrix do not satisfy the prescribed bound.')
      }
    }
    return(x)
  }

  i.trial = 0
  err = TRUE
  sd = 1
  while ( (i.trial < n.trials) && (err) ) {
    theta = stats::rnorm(length(i), sd = sd)
    if (!is.null(digits)) theta = round(theta, digits)
    tmpl[i] = theta
    x = structure(tmpl, class = c('stsp','ratm'), order = order)
    err = FALSE

    if ( !is.null(bpoles) ) {
      err = try((min(abs(poles(x, print_message = FALSE))) <= bpoles), silent = TRUE)
      if (inherits(err, 'try-error')) err = TRUE
    }
    if ( (!err) && (!is.null(bzeroes)) ) {
      err = try((min(abs(zeroes(x, print_message = FALSE))) <= bzeroes), silent = TRUE)
      if (inherits(err, 'try-error')) err = TRUE
    }
    i.trial = i.trial + 1
    sd = sd/1.1
  }
  if (err) {
    stop('Could not generate a suitable rational matrix with ', n.trials, ' trials!')
  }
  return(x)
}
# transpose for rational matrices ##############################################


t.polm = function(x) {
  x = unclass(x)
  x = aperm(x, c(2,1,3))
  return(polm(x))
}

t.lpolm = function(x) {
  x = unclass(x)
  min_deg = attr(x, which = "min_deg")
  x = aperm(x, c(2,1,3))
  return(lpolm(x, min_deg = min_deg))
}


t.lmfd = function(x) {
  rmfd(c = t(x$a), d = t(x$b))
}

t.rmfd = function(x) {
  lmfd(a = t(x$c), b = t(x$d))
}


t.stsp = function(x) {

  d = attr(x, 'order')
  m = d[1]
  n = d[2]
  s = d[3]

  x = unclass(x)
  A = t(x[iseq(1,s), iseq(1,s), drop = FALSE])
  C = t(x[iseq(1,s), iseq(s+1,s+n), drop = FALSE])
  B = t(x[iseq(s+1,s+m), iseq(1,s), drop = FALSE])
  D = t(x[iseq(s+1,s+m), iseq(s+1,s+n), drop = FALSE])

  x = stsp(A = A, B = B, C = C, D = D)
  return(x)
}

t.pseries = function(x) {

  x = aperm(unclass(x), c(2,1,3))
  x = structure(x, class = c('pseries','ratm'))
  return(x)
}

t.zvalues = function(x) {

  z = attr(x, 'z')
  x = aperm(unclass(x), c(2,1,3))
  x = structure(x, z = z, class = c('zvalues','ratm'))
  return(x)
}

# Hermitian Transpose ####

Ht = function(x) {
  UseMethod("Ht", x)
}

Ht.polm = function(x) {
  x = unclass(x)
  max_deg = dim(x)[3] - 1
  x = aperm(Conj(x), c(2, 1, 3))
  if (dim(x)[3] > 1) x = x[ , , (dim(x)[3]):1, drop = FALSE]
  return(structure(x, class = c("lpolm", 'ratm'), min_deg = -max_deg))
}


Ht.lpolm = function(x) {
  x = unclass(x)
  min_deg = attr(x, which = "min_deg")
  max_deg = dim(x)[3] - 1 + min_deg
  x = aperm(Conj(x), c(2, 1, 3))
  if (dim(x)[3] > 1) x = x[ , , (dim(x)[3]):1, drop = FALSE]
  return(structure(x, class = c("lpolm", 'ratm'), min_deg = -max_deg))
}

Ht.stsp = function(x) {
  d = attr(x, 'order')
  m = d[1]
  n = d[2]
  s = d[3]

  if ((m*n) == 0) {
    x = stsp(A = matrix(0, nrow = 0, ncol = 0),
             B = matrix(0, nrow = 0, ncol = m),
             C = matrix(0, nrow = n, ncol = 0),
             D = matrix(0, nrow = n, ncol = m))
    return(x)
  }

  ABCD = Conj(unclass(x))
  # if (is.complex(ABCD)) {
  #   stop('Hermitean transpose is only implemented for "real" state space models')
  # }
  A = ABCD[iseq(1,s), iseq(1,s), drop = FALSE]
  B = ABCD[iseq(1,s), iseq(s+1,s+n), drop = FALSE]
  C = ABCD[iseq(s+1,s+m), iseq(1,s), drop = FALSE]
  D = ABCD[iseq(s+1,s+m), iseq(s+1,s+n), drop = FALSE]

  if (s == 0) {
    x = stsp(A, B = t(C), C = t(B), D = t(D))
    return(x)
  }

  A = try(solve(A), silent = TRUE)
  if (inherits(A, 'try-error')) {
    stop('The state transition matrix "A" is singular.')
  }

  junk = -t( A %*% B)         # C -> - B' A^{-T}
  D = t(D) + junk %*% t(C)    # D -> D' - B' A^{-T} C'
  B = t( C %*% A )            # B -> A^{-T} C'
  A = t(A)                    # A -> A^{-T}
  C = junk

  x = stsp(A = A, B = B, C = C, D = D)
  return(x)
}


Ht.zvalues = function(x) {
  z = 1/Conj(attr(x,'z'))

  x = Conj(aperm(unclass(x), c(2,1,3)))
  x = structure(x, z = z, class = c('zvalues','ratm'))

  return(x)
}

# frequency response function ######################################################

zvalues = function(obj, z, f, n.f, sort.frequencies, ...){
  UseMethod("zvalues", obj)
}

# Internal function used by zvalues.____ methods
# Intended to check the inputs to zvalues() and bring them in a normalized form
get_z = function(z = NULL, f = NULL, n.f = 5, sort.frequencies = FALSE) {

  # If there are no values in z,
  # then either take the frequencies directly,
  # or create them from the "number of frequencies" n.f
  if (is.null(z)) {
    if (is.null(f)) {
      n.f = as.vector(as.integer(n.f))[1]
      stopifnot("Number of frequencies must be a non negative integer" = (n.f >= 0))
      if (n.f == 0) {
        f = double(0)
      } else {
        f = (0:(n.f-1))/n.f
      }
    }

    if (!is.numeric(f)) {
      stop('"f" must be a vector of (real) frequencies')
    } else {
      f = as.vector(f)
      n.f = length(f)
      z = exp(complex(imaginary = (-2*pi)*f))
    }
  }

  # Check non-NULL inputs z
  if (! (is.numeric(z) || is.complex(z))) stop('"z" must be a numeric or complex valued vector')
  z = as.vector(z)

  # Bring potential other arguments in line with z values
  n.f = length(z)
  f = -Arg(z)/(2*pi)

  # Sort if requested
  if ((sort.frequencies) && (n.f > 0)) {
    i = order(f)
    z = z[i]
    f = f[i]
  }

  return(list(z = z, f = f, n.f = n.f))
}

zvalues.default = function(obj, z = NULL, f = NULL, n.f = NULL, sort.frequencies = FALSE,  ...) {

  stopifnot("For the default method (constructing a zvalues object), *obj* needs to be an array, matrix, or vector." = any(class(obj) %in% c("matrix", "array")) || is.vector(obj),
            "Input *obj* must be numeric or complex" = (is.numeric(obj) || is.complex(obj)),
            "For the default method (constructing a zvalues object), *z*, *f*, or *n.f* needs to be supplied. Use of *z* is recommended." = any(!is.null(z), !is.null(f), !is.null(f)),
            "For the default method (constructing a zvalues object), *sort.frequencies* needs to be FALSE" = !sort.frequencies)

  if (is.vector(obj)) {
    dim(obj) = c(1,1,length(obj))
  }
  if (is.matrix(obj)) {
    dim(obj) = c(dim(obj),1)
  }

  stopifnot("could not coerce input parameter *obj* to a valid 3-D array!" = (is.array(obj)) && (length(dim(obj)) == 3))

  zz = get_z(z, f, n.f, sort.frequencies)
  z = zz$z
  f = zz$f
  n.f = zz$n.f

  d = dim(obj)
  m = d[1]
  n = d[2]

  stopifnot("Length of *obj* (i.e. length of vector, third dimension of array) needs to be equal to *n.f* or length of *z*, *f*." = (d[3] == n.f))

  # Empty zvalues object
  if (prod(d)*n.f == 0) {
    fr = array(0, dim = unname(c(m, n, n.f)))
    fr = structure(fr, z = z, class = c('zvalues','ratm'))
    return(fr)
  }

  fr = structure(obj, z = z, class = c('zvalues','ratm'))
  return(fr)
}


zvalues.polm = function(obj, z = NULL, f = NULL, n.f = 5, sort.frequencies = FALSE,  ...) {

  # Obtain checked and normalized values
  zz = get_z(z, f, n.f, sort.frequencies)
  z = zz$z # values on the unit circle (if z,f NULL)
  f = zz$f # between -0.5 and 0.5
  n.f = zz$n.f

  # Integer-valued parameters
  obj = unclass(obj)
  d = dim(obj)
  m = d[1]
  n = d[2]
  p = d[3] - 1

  # Check if dimension, degree, or number of frequencies is zero
  if ((m*n*(p+1)*n.f) == 0) {
    fr = array(0, dim = unname(c(m, n, n.f)))
    fr = structure(fr, z = z, class = c('zvalues','ratm'))
    return(fr)
  }

  # Copy highest coefficient n.f times in an array as initial value
  fr = array(obj[,, p+1], dim = unname(c(m, n, n.f)))

  if (p == 0) {
    fr = structure(fr, z = z, class = c('zvalues','ratm'))
    return(fr)
  }

  # Array of dim ( out, in, n_points): every element (i,j) of a polynomial will be evaluated at n_f values
  zz = aperm(array(z, dim = c(n.f, m, n)), c(2,3,1))

  # Horner-scheme for evaluation
  for (j in (p:1)) {
    fr = fr*zz + array(obj[,,j], dim = unname(c(m,n,n.f)))
  }

  fr = structure(fr, z = z, class = c('zvalues','ratm'))
  return(fr)
}

zvalues.lpolm = function(obj, z = NULL, f = NULL, n.f = 5, sort.frequencies = FALSE,  ...) {

  # Obtain checked and normalized values
  zz = get_z(z, f, n.f, sort.frequencies)
  z = zz$z # values on the unit circle (if z,f NULL)
  f = zz$f # between -0.5 and 0.5
  n.f = zz$n.f

  # Integer-valued parameters
  d = dim(obj)
  m = d[1]
  n = d[2]
  p = d[3] # Dimensions are retrieved from polm object, NOT from unclassed polm obj = array
  min_deg = d[4]

  obj = obj %>% unclass()

  # Check if a dimension or number of frequencies is zero (degree not part of this in contrast to ___.polm() function)
  if ((m*n*n.f) == 0) {
    fr = array(0, dim = unname(c(m, n, n.f)))
    fr = structure(fr, z = z, class = c('zvalues','ratm'))
    return(fr)
  }

  # Copy highest coefficient n.f times in an array as initial value
  fr = array(obj[,, p+1-min_deg], dim = unname(c(m, n, n.f)))

  if (p == 0 && min_deg == 0) {
    fr = structure(fr, z = z, class = c('zvalues','ratm'))
    return(fr)
  }

  # Array of dim ( out, in, n_points): every element (i,j) of a polynomial will be evaluated at n_f values
  zz = aperm(array(z, dim = c(n.f, m, n)), c(2,3,1))

  # Horner-scheme for evaluation
  for (j in ((p-min_deg):1)) {
    fr = fr*zz + array(obj[,,j], dim = unname(c(m,n,n.f)))
  }

  # Every polynomial (i,j) - evaluated at n_f points - needs to be premultiplied by z_0^{min_deg} (where z_0 denotes one particular point of evaluation)
  zz_factor_min_deg = aperm(array(z, dim = c(n.f, m, n)), perm = c(2,3,1))^min_deg
  fr = fr * zz_factor_min_deg

  fr = structure(fr, z = z, class = c('zvalues','ratm'))
  return(fr)
}

zvalues.lmfd = function(obj, z = NULL, f = NULL, n.f = 5, sort.frequencies = FALSE,  ...) {
  zz = get_z(z, f, n.f, sort.frequencies)
  z = zz$z
  f = zz$f
  n.f = zz$n.f

  d = attr(obj,'order')
  m = d[1]
  n = d[2]
  p = d[3]
  q = d[4]

  if ((m<1) || (p<0)) stop('obj is not a valid "rldm" object')

  if ((m*n*(q+1)*n.f) == 0) {
    fr = array(0, dim = unname(c(m,n,n.f)))
    fr = structure(fr, z = z, class = c('zvalues','ratm'))
    return(fr)
  }

  fr_a = unclass(zvalues(obj$a, z = z))
  fr = unclass(zvalues(obj$b, z = z))

  for (i in (1:n.f)) {
    # x = try(solve(matrix(fr_a[,,i], nrow = m, ncol = m),
    #               matrix(fr[,,i], nrow = m, ncol = n)), silent = TRUE)
    # # print(x)
    # if (!inherits(x, 'try-error')) {
    #   fr[,,i] = x
    # } else {
    #   fr[,,i] = NA
    # }
    fr[,,i] = tryCatch(solve(matrix(fr_a[,,i], nrow = m, ncol = m),
                             matrix(fr[,,i], nrow = m, ncol = n)),
                       error = function(e) NA_real_)
  }

  fr = structure(fr, z = z, class = c('zvalues','ratm'))
  return(fr)
}

zvalues.rmfd = function(obj, z = NULL, f = NULL, n.f = 5, sort.frequencies = FALSE,  ...) {
  zz = get_z(z, f, n.f, sort.frequencies)
  z = zz$z
  f = zz$f
  n.f = zz$n.f

  d = attr(obj,'order')
  m = d[1]
  n = d[2]
  p = d[3]
  q = d[4]

  if ((m<1) || (p<0)) stop('obj is not a valid "rldm" object')

  if ((m*n*(q+1)*n.f) == 0) {
    fr = array(0, dim = unname(c(m,n,n.f)))
    fr = structure(fr, z = z, class = c('zvalues','ratm'))
    return(fr)
  }

  # compute first transposed values: t(c(z))^(-1) t(d(z))
  fr_c = aperm(unclass(zvalues(obj$c, z = z)), c(2,1,3))
  fr = aperm(unclass(zvalues(obj$d, z = z)), c(2,1,3))

  for (i in (1:n.f)) {
    # fr[,,i] = tryCatch(matrix(fr[,,i], nrow = m, ncol = n) %*%
    #                      solve(matrix(fr_c[,,i], nrow = n, ncol = n)),
    #                    error = function(e) NA)
    fr[,,i] = tryCatch(solve(matrix(fr_c[,,i], nrow = n, ncol = n),
                             matrix(fr[,,i], nrow = n, ncol = m)),
                       error = function(e) NA_real_)
  }
  # undo transposition
  fr = aperm(fr, c(2,1,3))

  fr = structure(fr, z = z, class = c('zvalues','ratm'))
  return(fr)
}

zvalues.stsp = function(obj, z = NULL, f = NULL, n.f = 5, sort.frequencies = FALSE,  ...) {
  zz = get_z(z, f, n.f, sort.frequencies)
  z = zz$z
  f = zz$f
  n.f = zz$n.f

  d = attr(obj, 'order')
  m = d[1]
  n = d[2]
  s = d[3]

  if ((m*n*n.f) == 0) {
    fr = array(0, dim = unname(c(m,n,n.f)))
    fr = structure(fr, z = z, class = c('zvalues','ratm'))
    return(fr)
  }

  ABCD = unclass(obj)
  A = ABCD[iseq(1,s), iseq(1,s), drop = FALSE]
  B = ABCD[iseq(1,s), iseq(s+1,s+n), drop = FALSE]
  C = ABCD[iseq(s+1,s+m), iseq(1,s), drop = FALSE]
  D = ABCD[iseq(s+1,s+m), iseq(s+1,s+n), drop = FALSE]

  fr = array(D, dim = unname(c(m,n,n.f)))

  if (s == 0) {
    fr = structure(fr, z = z, class = c('zvalues','ratm'))
    return(fr)
  }

  Is = diag(s)
  for (i in (1:n.f)) {
    # zz = z[i]^{-1}
    # x = try(C %*% solve(diag(zz, nrow = s, ncol = s) - A, B), silent = TRUE)
    # if (!inherits(x, 'try-error')) {
    #   fr[,,i] = fr[,,i] + x
    # } else {
    #   fr[,,i] = NA
    # }
    fr[,,i] = fr[,,i] + z[i]*tryCatch(C %*% solve(Is - z[i]*A, B),
                                 error = function(e) NA_real_)
  }

  fr = structure(fr, z = z, class = c('zvalues','ratm'))
  return(fr)
}

# Evaluate at one point only ####

zvalue = function(obj, z, ...){
  UseMethod("zvalue", obj)
}

zvalue.lpolm = function(obj, z = NULL, ...) {

  stopifnot(length(z) == 1)
  stopifnot(inherits(obj, "ratm"))

  # Integer-valued parameters
  d = dim(obj)
  m = d[1]
  n = d[2]
  p = d[3] # Dimensions are retrieved from polm object, NOT from unclassed polm obj = array
  min_deg = d[4]

  obj = obj %>% unclass()

  # Check if a dimension or number of frequencies is zero (degree not part of this in contrast to ___.polm() function)
  if ((m*n) == 0) {
    return(matrix(0, m, n))
  }

  # Copy highest coefficient in a matrix as initial value
  out = matrix(obj[,, p+1-min_deg], m, n)

  if (p == 0 && min_deg == 0) {
    return(out)
  }

  # Horner-scheme for evaluation
  for (j in ((p-min_deg):1)) {
    out = out*z + matrix(obj[,,j], m, n)
  }

  # Every polynomial (i,j) needs to be premultiplied by z_0^{min_deg} (where z_0 denotes one particular point of evaluation)
  factor_min_deg = z^min_deg
  return(out * factor_min_deg)
}

zvalue.polm = function(obj, z = NULL, ...) {

  stopifnot(length(z) == 1)
  stopifnot(inherits(obj, "ratm"))

  # Integer-valued parameters
  d = dim(obj)
  m = d[1]
  n = d[2]
  p = d[3] # Dimensions are retrieved from polm object, NOT from unclassed polm obj = array

  obj = obj %>% unclass()

  # Check if a dimension or number of frequencies is zero (degree not part of this in contrast to ___.polm() function)
  if ((m*n*(p+1)) == 0) {
    return(matrix(0, m, n))
  }

  # Copy highest coefficient in a matrix as initial value
  out = matrix(obj[,, p+1], m, n)

  if (p == 0) {
    return(out)
  }

  # Horner-scheme for evaluation
  for (j in (p:1)) {
    out = out*z + matrix(obj[,,j], m, n)
  }

  return(out)
}

zvalue.lmfd = function(obj, z = NULL, ...) {

  d = attr(obj,'order')
  m = d[1]
  n = d[2]
  p = d[3]
  q = d[4]

  if (p == -1) {
    stop("a(z) matrix is identically zero.")
  }

  if ((m*n*(q+1)) == 0) {
    return(matrix(0, m, n))
  }

  out_a = unclass(zvalue(obj$a, z = z))
  out_b = unclass(zvalue(obj$b, z = z))

  out_a_inv = try(solve(out_a))
  if(inherits(out_a_inv, "try-error")) {
    stop("a(z) not invertible for input z.")
  }

  return(out_a_inv %*% out_b)
}

zvalue.rmfd = function(obj, z = NULL, ...) {
  d = attr(obj,'order')
  m = d[1]
  n = d[2]
  p = d[3]
  q = d[4]

  if (p == -1) {
    stop("c(z) matrix is identically zero.")
  }

  if ((m*n*(q+1)) == 0) {
    return(matrix(0, m, n))
  }

  out_c = unclass(zvalue(obj$c, z = z))
  out_d = unclass(zvalue(obj$d, z = z))

  out_c_inv = try(solve(out_c))
  if(inherits(out_c_inv, "try-error")) {
    stop("c(z) not invertible for input z.")
  }

  return(out_d %*% out_c_inv)
}

zvalue.stsp = function(obj, z = NULL,  ...) {

  d = attr(obj, 'order')
  m = d[1]
  n = d[2]
  s = d[3]

  if ((m*n) == 0) {
    return(matrix(0, m, n))
  }

  ABCD = unclass(obj)
  A = ABCD[iseq(1,s), iseq(1,s), drop = FALSE]
  B = ABCD[iseq(1,s), iseq(s+1,s+n), drop = FALSE]
  C = ABCD[iseq(s+1,s+m), iseq(1,s), drop = FALSE]
  D = ABCD[iseq(s+1,s+m), iseq(s+1,s+n), drop = FALSE]

  out = matrix(D, m, n)

  if (s == 0 || z == 0) {
    return(out)
  }

  z = z^{-1}
  x = try(C %*% solve(diag(z, nrow = s, ncol = s) - A, B), silent = TRUE)
  if (!inherits(x, 'try-error')) {
    out = out + x
  } else {
    stop("State space object has a pole at input z.")
  }

  return(out)
}
