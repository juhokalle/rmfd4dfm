#include <RcppArmadillo.h>
using namespace Rcpp;
//' Kalman filter smoother
//'
//' @param A n x r observation system matrix
//' @param C r x r state transition matrix
//' @param R n x n noise covariance matrix
//' @param Q q x q error covariance matrix
//' @param data_wide n x T data matrix
//' @param P1 r x r initial state variance matrix
//' @param dim_state value r
//' @param only_ll boolean, return only log-likelihood calculated with KFss
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List kfks_cpp(const arma::mat& A, const arma::mat& C,
                    const arma::mat& R, const arma::mat& Q,
                    const arma::mat& data_wide, const arma::mat& P1,
                    double dim_state, bool only_ll) {

  double n_obs = data_wide.n_cols;
  double dim_out = data_wide.n_rows;
  double loglik = n_obs * dim_out * log(2 * arma::datum::pi);
  double lndet;
  double sign;
  double ll2;

  arma::mat s_t = arma::zeros<arma::mat>(dim_state, n_obs);
  arma::mat s_t_T = arma::zeros<arma::mat>(dim_state, n_obs);
  arma::mat e_t = arma::zeros<arma::mat>(dim_out, n_obs);
  arma::mat e_t_T = arma::zeros<arma::mat>(dim_out, n_obs);
  arma::mat v_t_T = arma::zeros<arma::mat>(dim_state, n_obs);
  arma::cube P_t = arma::zeros<arma::cube>(dim_state, dim_state, n_obs);
  arma::cube P_t_T = arma::zeros<arma::cube>(dim_state, dim_state, n_obs);
  arma::cube C_t_T = arma::zeros<arma::cube>(dim_state, dim_state, n_obs);
  arma::cube K_t = arma::zeros<arma::cube>(dim_state, dim_out, n_obs);
  arma::cube L_t = arma::zeros<arma::cube>(dim_state, dim_state, n_obs);
  arma::mat S_t = arma::zeros<arma::mat>(dim_out, dim_out);
  arma::cube S_t_inv = arma::zeros<arma::cube>(dim_out, dim_out, n_obs);
  arma::mat r_t = arma::zeros<arma::mat>(dim_state, n_obs+1);
  arma::cube N_t = arma::zeros<arma::cube>(dim_state, dim_state, n_obs+1);


  P_t.slice(0) = P1;

  for (int ix_t = 0; ix_t < n_obs; ix_t++) {

    S_t = C*P_t.slice(ix_t)*C.t() + R;
    S_t_inv.slice(ix_t) = inv(S_t);
    log_det(lndet, sign, S_t);
    ll2 = as_scalar(e_t.col(ix_t).t() * S_t_inv.slice(ix_t) * e_t.col(ix_t));

    e_t.col(ix_t) = data_wide.col(ix_t) - C * s_t.col(ix_t);

    K_t.slice(ix_t) = A*P_t.slice(ix_t)*C.t()*S_t_inv.slice(ix_t);
    L_t.slice(ix_t) = A - K_t.slice(ix_t)*C;
    loglik += lndet + ll2;

    if(ix_t<(n_obs-1)){
      s_t.col(ix_t+1) = A*s_t.col(ix_t) + K_t.slice(ix_t)*e_t.col(ix_t);
      P_t.slice(ix_t+1) = A*P_t.slice(ix_t)*A.t() + Q - K_t.slice(ix_t) * S_t * K_t.slice(ix_t).t();
    }
  }
  // return early if only log-likelihood is of interest
  if(only_ll) return Rcpp::List::create(Rcpp::Named("loglik") = -loglik/(2*n_obs));

  // P_{t-1|t-1}, note the indexing
  const arma::mat P_tt_hlp = P_t.slice(n_obs-2) - P_t.slice(n_obs-2)*C.t()*S_t_inv.slice(n_obs-2)*C*P_t.slice(n_obs-2);
  C_t_T.slice(n_obs-1) = P_tt_hlp*A.t()*(arma::eye<arma::mat>(dim_state, dim_state) - C.t() * K_t.slice(n_obs-1).t());

  for (int ix_t = 1; ix_t < n_obs+1; ix_t++){
    N_t.slice(n_obs-ix_t) = C.t()*S_t_inv.slice(n_obs-ix_t)*C + L_t.slice(n_obs-ix_t).t()*N_t.slice(n_obs-ix_t+1)*L_t.slice(n_obs-ix_t);
    r_t.col(n_obs-ix_t) = C.t()*S_t_inv.slice(n_obs-ix_t)*e_t.col(n_obs-ix_t) + L_t.slice(n_obs-ix_t).t() * r_t.col(n_obs-ix_t+1);

    s_t_T.col(n_obs-ix_t) = s_t.col(n_obs-ix_t) + P_t.slice(n_obs-ix_t)*r_t.col(n_obs-ix_t);
    P_t_T.slice(n_obs-ix_t) = P_t.slice(n_obs-ix_t) - P_t.slice(n_obs-ix_t)*N_t.slice(n_obs-ix_t)*P_t.slice(n_obs-ix_t).t();

    e_t_T.col(n_obs-ix_t) = R*(S_t_inv.slice(n_obs-ix_t)*e_t.col(n_obs-ix_t)-K_t.slice(n_obs-ix_t).t()*r_t.col(n_obs-ix_t+1));
    v_t_T.col(n_obs-ix_t) = Q*r_t.col(n_obs-ix_t+1);

    if(ix_t>1 && ix_t<n_obs){
      C_t_T.slice(n_obs-ix_t) = P_t.slice(n_obs-ix_t-1) * L_t.slice(n_obs-ix_t-1).t()*(arma::eye<arma::mat>(dim_state, dim_state) - N_t.slice(n_obs-ix_t)*P_t.slice(n_obs-ix_t));
    }
  }

  return Rcpp::List::create(Rcpp::Named("P_t_T") = P_t_T,
                            Rcpp::Named("C_t_T") = C_t_T,
                            Rcpp::Named("s_t_T") = s_t_T,
                            Rcpp::Named("e_t_T") = e_t_T,
                            Rcpp::Named("v_t_T") = v_t_T,
                            Rcpp::Named("N_t") = N_t,
                            Rcpp::Named("L_t") = L_t,
                            Rcpp::Named("loglik") = -loglik/(2*n_obs));

}
