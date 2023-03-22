#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' @keywords internal
// [[Rcpp::export]]
arma::mat l_esreg_covariance(
    arma::mat xq, arma::mat xe, arma::colvec xbq, arma::colvec xbe,
    arma::colvec G1_prime_xq, arma::colvec G2_xe, arma::colvec G2_prime_xe,
    arma::colvec density, arma::colvec conditional_variance, double alpha) {
  try {
    int n = xq.n_rows;
    int kq = xq.n_cols;
    int ke = xe.n_cols;

    // Define some 0-matrices
    arma::mat lambda_11 = arma::zeros<arma::mat>(kq, kq);
    arma::mat lambda_22 = arma::zeros<arma::mat>(ke, ke);

    arma::mat C_11 = arma::zeros<arma::mat>(kq, kq);
    arma::mat C_12 = arma::zeros<arma::mat>(ke, kq);
    arma::mat C_22 = arma::zeros<arma::mat>(ke, ke);

    arma::mat lambda = arma::zeros<arma::mat>(kq+ke, kq+ke);
    arma::mat C = arma::zeros<arma::mat>(kq+ke, kq+ke);

    // Compute the matrix elements
    for (int i = 0; i < n; i++) {
      arma::mat xqi = xq.row(i);
      arma::mat xei = xe.row(i);
      arma::mat xxq = xqi.t() * xqi;
      arma::mat xxe = xei.t() * xei;
      arma::mat xxeq = xei.t() * xqi;

      lambda_11 += 1/alpha * xxq * density(i) * (alpha*G1_prime_xq(i) + G2_xe(i));
      lambda_22 += xxe * G2_prime_xe(i);

      C_11 += (1-alpha)/alpha * xxq * pow(alpha*G1_prime_xq(i) + G2_xe(i), 2);
      C_12 += (1-alpha)/alpha * xxeq * (xbq(i) - xbe(i)) *
        (alpha*G1_prime_xq(i) + G2_xe(i)) * G2_prime_xe(i);
      C_22 += xxe * pow(G2_prime_xe(i), 2) * (conditional_variance(i)/alpha +
        (1-alpha)/alpha * pow(xbq(i) - xbe(i), 2));
    }

    // Fill the matrices
    C.submat(0, 0, kq-1, kq-1) = C_11 / n;
    C.submat(kq, 0, kq+ke-1, kq-1) = C_12 / n;
    C.submat(0, kq, kq-1, kq+ke-1) = C_12.t() / n;
    C.submat(kq, kq, kq+ke-1, kq+ke-1) = C_22 / n;

    lambda.submat(0, 0, kq-1, kq-1) = lambda_11 / n;
    lambda.submat(kq, kq, kq+ke-1, kq+ke-1) = lambda_22 / n;

    // Compute the covariance
    arma::mat cov = inv(lambda) * C * inv(lambda) / n;

    return cov;
  } catch(...) {
    Rcpp::stop("Cannot compute the covariance!");
  }
}
