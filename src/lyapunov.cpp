// lyapunov.cpp


// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>

using namespace Rcpp;

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// To make your C++ code callable from C++ code in other packages.
// This will generate a header file, inst/include/mypackage.h that
// can be included by other packages
// [[Rcpp::interfaces(r, cpp)]]

// via the exports attribute we tell Rcpp to make this function
// available from R

// solve lyapunov equation P = A * P * A' + Q
// just for internal use, no checks!
// if A is not stable, return FALSE
//
// RcppArmadillo does not fully implement "pass by reference for complex vectors"
// therefore I simply split the eigenvalues (lambda) into the real and imaginary part

// [[Rcpp::export]]
bool lyapunov_cpp(const arma::mat& A, const arma::mat& Q, arma::mat& P,
              arma::vec& lambda_r, arma::vec& lambda_i, bool stop_if_non_stable) {

  int m = A.n_rows;

  arma::cx_mat cA = arma::cx_mat(A, arma::mat(m,m,arma::fill::zeros));
  arma::cx_mat cQ = arma::cx_mat(Q, arma::mat(m,m,arma::fill::zeros));
  arma::cx_mat cU = cQ;
  arma::cx_mat cS = cQ;

  // call Schur decomposition of cA, such that cA = cU*cS*cU.t() and cS is upper triangular!
  bool ok = arma::schur(cU, cS, cA);
  if (!ok) {
    stop("RcppArmadillo \"schur\" algorithm failed");
  }
  bool is_stable = (max(abs(cS.diag())) < 1);
  lambda_r = real(cS.diag());
  lambda_i = imag(cS.diag());

  // Rcout << std::endl << stop_if_non_stable << " " << is_stable << " " << std::endl;
  // Rcout << std::endl << "P" << std::endl << P << std::endl;

  // if A is not stable, return FALSE
  if ( (stop_if_non_stable)  and  (!is_stable) ) {
    return FALSE;
  }

  // transform Q
  cQ = cU.t() * cQ * cU;
  
  // Rcout << std::endl << "Q" << std::endl << cQ << std::endl << std::endl << "S" << std::endl << cS << std::endl;

  
  // RcppArmadillo syntax for submatrices:
  //   X.submat( first_row, first_col, last_row, last_col )

  // cQ is recursively overwritten with the desired solution P
  int i;
  for (i = (m-1); i>0; i--) {
    // P22 = S22 * P22 * conj(S22) + Q22 
    // Q22 <- P22 = Q22 / ( 1- S22*conj(S22)) )
    cQ(i,i) = cQ(i,i) / (1.0 - cS(i,i) * conj(cS(i,i)));

    // [Q11, Q12] <- [Q11, Q12] + S12 * P22 * [conj(S21), conj(S22)]  
    cQ.submat(0,0,i-1,i) = cQ.submat(0,0,i-1,i) + ( cS.submat(0,i,i-1,i) * cQ(i,i) ) * cS.submat(0,i,i,i).t();

    // (I - S11 *conj(S22))
    cA.submat(0,0,i-1,i-1) = (-conj(cS(i,i))) * cS.submat(0,0,i-1,i-1);
    cA.submat(0,0,i-1,i-1).diag() += 1;
    //    cA.print();

    // Q12 <- P12 = (I - S11 B22)^-1 * Q12
    cQ.submat(0,i,i-1,i) = solve(trimatu(cA.submat(0,0,i-1,i-1)), cQ.submat(0,i,i-1,i));
    // Q21 <- P21 = conj(P12)
    cQ.submat(i,0,i,i-1) = cQ.submat(0,i,i-1,i).t();

    // S11 * P12 * conj(S12)
    cA.submat(0,0,i-1,i-1) = ( cS.submat(0,0,i-1,i-1) * cQ.submat(0,i,i-1,i) ) * ( cS.submat(0,i,i-1,i).t() );

    // Q11 <- Q11 + S11 * P12 * conj(S12) + conj(S11 * P12 * conj(S12))
    cQ.submat(0,0,i-1,i-1) = cQ.submat(0,0,i-1,i-1) + cA.submat(0,0,i-1,i-1) + cA.submat(0,0,i-1,i-1).t();
  }
  i = 0;
  cQ(i,i) = cQ(i,i) / (1.0 - cS(i,i) * conj(cS(i,i)));

  // Rcout << std::endl << "Q" << std::endl << cQ << std::endl;
  
  // retransform Q
  cQ = cU * cQ * cU.t();
  // make sure that Q is Hermitean
  cQ = (cQ + cQ.t())/2;
  // make real
  // the conversion conv_to<arma::mat>::from(cQ) does not work !?!
  // P = arma::conv_to<arma::mat>::from(cQ);
  P = real(cQ);

  // Rcout << std::endl << "P" << std::endl << P << std::endl;

  return is_stable;
}


// [[Rcpp::export]]
bool lyapunov_Jacobian_cpp(const arma::mat& A, const arma::mat& Q, arma::mat& P,
                           const arma::mat& dA, const arma::mat& dQ, arma::mat& J, 
                           arma::vec& lambda_r, arma::vec& lambda_i, bool stop_if_non_stable) {
  
  int m = A.n_rows;
  int n = dA.n_cols;
  
  arma::cx_mat cA = arma::cx_mat(A, arma::mat(m,m,arma::fill::zeros));
  arma::cx_mat cQ = arma::cx_mat(Q, arma::mat(m,m,arma::fill::zeros));
  arma::cx_mat cU = cQ;
  arma::cx_mat cS = cQ;
  
  // call Schur decomposition of cA, such that cA = cU*cS*cU.t() and cS is upper triangular!
  bool ok = arma::schur(cU, cS, cA);
  if (!ok) {
    stop("RcppArmadillo \"schur\" algorithm failed");
  }
  bool is_stable = (max(abs(cS.diag())) < 1);
  lambda_r = real(cS.diag());
  lambda_i = imag(cS.diag());
  
  // Rcout << std::endl << stop_if_non_stable << " " << is_stable << " " << std::endl;
  // Rcout << std::endl << "P" << std::endl << P << std::endl;
  
  // if A is not stable, return FALSE
  if ( (stop_if_non_stable)  and  (!is_stable) ) {
    return FALSE;
  }
  
  // compute solution P of the lyapunov equation //////////////////////////////
  
  // transform Q
  cQ = cU.t() * cQ * cU;
  
  // Rcout << std::endl << "Q" << std::endl << cQ << std::endl << std::endl << "S" << std::endl << cS << std::endl;

  // RcppArmadillo syntax for submatrices:
  //   X.submat( first_row, first_col, last_row, last_col )
  
  // cQ is recursively overwritten with the desired solution P
  int i;
  for (i = (m-1); i>0; i--) {
    // P22 = S22 * P22 * conj(S22) + Q22 
    // Q22 <- P22 = Q22 / ( 1- S22*conj(S22)) )
    cQ(i,i) = cQ(i,i) / (1.0 - cS(i,i) * conj(cS(i,i)));
    
    // [Q11, Q12] <- [Q11, Q12] + S12 * P22 * [conj(S21), conj(S22)]  
    cQ.submat(0,0,i-1,i) = cQ.submat(0,0,i-1,i) + ( cS.submat(0,i,i-1,i) * cQ(i,i) ) * cS.submat(0,i,i,i).t();
    
    // (I - S11 *conj(S22))
    cA.submat(0,0,i-1,i-1) = (-conj(cS(i,i))) * cS.submat(0,0,i-1,i-1);
    cA.submat(0,0,i-1,i-1).diag() += 1;
    //    cA.print();
    
    // Q12 <- P12 = (I - S11 B22)^-1 * Q12
    cQ.submat(0,i,i-1,i) = solve(trimatu(cA.submat(0,0,i-1,i-1)), cQ.submat(0,i,i-1,i));
    // Q21 <- P21 = conj(P12)
    cQ.submat(i,0,i,i-1) = cQ.submat(0,i,i-1,i).t();
    
    // S11 * P12 * conj(S12)
    cA.submat(0,0,i-1,i-1) = ( cS.submat(0,0,i-1,i-1) * cQ.submat(0,i,i-1,i) ) * ( cS.submat(0,i,i-1,i).t() );
    
    // Q11 <- Q11 + S11 * P12 * conj(S12) + conj(S11 * P12 * conj(S12))
    cQ.submat(0,0,i-1,i-1) = cQ.submat(0,0,i-1,i-1) + cA.submat(0,0,i-1,i-1) + cA.submat(0,0,i-1,i-1).t();
  }
  i = 0;
  cQ(i,i) = cQ(i,i) / (1.0 - cS(i,i) * conj(cS(i,i)));
  
  // Rcout << std::endl << "Q" << std::endl << cQ << std::endl;
  
  // retransform Q
  cQ = cU * cQ * cU.t();
  // make sure that Q is Hermitean
  cQ = (cQ + cQ.t())/2;
  // make real
  // the conversion conv_to<arma::mat>::from(cQ) does not work !?!
  // P = arma::conv_to<arma::mat>::from(cQ);
  P = real(cQ);
  // compute solution P of the lyapunov equation //////////////////////////////
  

  // compute derivatives of the solution P of the lyapunov equation ///////////
  // dP = A dP A' + dA P A' + A P dA' + dQ   
  
  arma::mat dQk = arma::mat(m,m);
  
  for (int k = 0; k < n; k++) {
    
    dQk = reshape(dA.col(k), m, m);
    dQk = dQk * P * A.t();
    dQk = dQk + dQk.t() + reshape(dQ.col(k), m, m);
      
    cQ = arma::cx_mat(dQk, arma::mat(m,m,arma::fill::zeros));
    
    // transform Q
    cQ = cU.t() * cQ * cU;

    // cQ is recursively overwritten with the desired solution P
    for (i = (m-1); i>0; i--) {
      // P22 = S22 * P22 * conj(S22) + Q22 
      // Q22 <- P22 = Q22 / ( 1- S22*conj(S22)) )
      cQ(i,i) = cQ(i,i) / (1.0 - cS(i,i) * conj(cS(i,i)));
      
      // [Q11, Q12] <- [Q11, Q12] + S12 * P22 * [conj(S21), conj(S22)]  
      cQ.submat(0,0,i-1,i) = cQ.submat(0,0,i-1,i) + ( cS.submat(0,i,i-1,i) * cQ(i,i) ) * cS.submat(0,i,i,i).t();
      
      // (I - S11 *conj(S22))
      cA.submat(0,0,i-1,i-1) = (-conj(cS(i,i))) * cS.submat(0,0,i-1,i-1);
      cA.submat(0,0,i-1,i-1).diag() += 1;
      //    cA.print();
      
      // Q12 <- P12 = (I - S11 B22)^-1 * Q12
      cQ.submat(0,i,i-1,i) = solve(trimatu(cA.submat(0,0,i-1,i-1)), cQ.submat(0,i,i-1,i));
      // Q21 <- P21 = conj(P12)
      cQ.submat(i,0,i,i-1) = cQ.submat(0,i,i-1,i).t();
      
      // S11 * P12 * conj(S12)
      cA.submat(0,0,i-1,i-1) = ( cS.submat(0,0,i-1,i-1) * cQ.submat(0,i,i-1,i) ) * ( cS.submat(0,i,i-1,i).t() );
      
      // Q11 <- Q11 + S11 * P12 * conj(S12) + conj(S11 * P12 * conj(S12))
      cQ.submat(0,0,i-1,i-1) = cQ.submat(0,0,i-1,i-1) + cA.submat(0,0,i-1,i-1) + cA.submat(0,0,i-1,i-1).t();
    }
    i = 0;
    cQ(i,i) = cQ(i,i) / (1.0 - cS(i,i) * conj(cS(i,i)));
    
    // Rcout << std::endl << "Q" << std::endl << cQ << std::endl;
    
    // retransform Q
    cQ = cU * cQ * cU.t();
    // make sure that Q is Hermitean
    cQ = (cQ + cQ.t())/2;
    // make real
    // the conversion conv_to<arma::mat>::from(cQ) does not work !?!
    // P = arma::conv_to<arma::mat>::from(cQ);
    J.col(k) = vectorise(real(cQ));
  }
  
  
  return is_stable;
}

