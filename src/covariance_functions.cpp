#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat l_esreg_covariance(arma::mat x, arma::colvec xq, arma::colvec xe,
                             arma::colvec G1_prime_xq, arma::colvec G2_xe, arma::colvec G2_prime_xe,
                             arma::colvec density, arma::colvec conditional_variance, double alpha) {
  try {
    int n = x.n_rows;
    int k = x.n_cols;

    // Define some 0-matrices
    arma::mat lambda_11 = arma::zeros<arma::mat>(k, k);
    arma::mat lambda_22 = arma::zeros<arma::mat>(k, k);

    arma::mat C_11 = arma::zeros<arma::mat>(k, k);
    arma::mat C_22 = arma::zeros<arma::mat>(k, k);
    arma::mat C_12 = arma::zeros<arma::mat>(k, k);

    arma::mat lambda = arma::zeros<arma::mat>(2*k, 2*k);
    arma::mat C = arma::zeros<arma::mat>(2*k, 2*k);

    arma::mat xi, xx;

    // Compute the matrix elements
    for (int i = 0; i < n; i++) {
      xi = x.row(i);
      xx = xi.t() * xi;
      lambda_11 += xx * (alpha * G1_prime_xq(i) + G2_xe(i)) / alpha * density(i);
      lambda_22 += xx * G2_prime_xe(i);
      C_11 += (1-alpha) / alpha * xx * pow((alpha * G1_prime_xq(i) + G2_xe(i)), 2);
      C_12 += (1-alpha) / alpha * xx * (alpha * G1_prime_xq(i) + G2_xe(i)) * G2_prime_xe(i) * (xq(i) - xe(i));
      C_22 += xx * pow(G2_prime_xe(i), 2) * (conditional_variance(i) / alpha + (1-alpha) / alpha * pow(xq(i) - xe(i), 2));
    }

    // Fill the big matrices
    C.submat(0, 0, k-1, k-1) = C_11 / n;
    C.submat(0, k, k-1, 2*k-1) = C_12 / n;
    C.submat(k, 0, 2*k-1, k-1) = C_12 / n;
    C.submat(k, k, 2*k-1, 2*k-1) = C_22 / n;

    lambda.submat(0, 0, k-1, k-1) = lambda_11 / n;
    lambda.submat(k, k, 2*k-1, 2*k-1) = lambda_22 / n;

    // Compute the covariance
    arma::mat cov = inv(lambda) * C * inv(lambda) / n;

    return cov;
  } catch(...) {
    Rcpp::stop("Cannot compute the covariance!");
  }
}


// [[Rcpp::export]]
arma::mat l_esreg_twostep_covariance(arma::mat x, arma::colvec xq, arma::colvec xe,
                                     arma::colvec density, arma::colvec conditional_variance, double alpha) {
  try {
    int n = x.n_rows;
    int k = x.n_cols;

    // Define some 0-matrices
    arma::mat C_11 = arma::zeros<arma::mat>(k, k);
    arma::mat C_22 = arma::zeros<arma::mat>(k, k);

    arma::mat G_11 = arma::zeros<arma::mat>(k, k);
    arma::mat G_12 = arma::zeros<arma::mat>(k, k);
    arma::mat G_22 = arma::zeros<arma::mat>(k, k);

    arma::mat C = arma::zeros<arma::mat>(2*k, 2*k);
    arma::mat B = arma::zeros<arma::mat>(2*k, 2*k);

    arma::mat xi, xx;

    // Compute the matrix elements
    for (int i = 0; i < n; i++) {
      xi = x.row(i);
      xx = xi.t() * xi;

      C_11 += alpha * xx * conditional_variance(i);
      C_22 += alpha * (1 - alpha) * xx;

      G_11 += -alpha * xx;
      G_12 += xx * density(i) * (xe(i) - xq(i));
      G_22 += xx * density(i);
    }

    // Fill the big matrices
    C.submat(0, 0, k-1, k-1) = C_11 / n;
    C.submat(k, k, 2*k-1, 2*k-1) = C_22 / n;

    B.submat(0, 0, k-1, k-1) = -inv(G_11 / n);
    B.submat(0, k, k-1, 2*k-1) = -inv(G_11 / n) * (G_12 / n) * inv(G_22 / n);
    B.submat(k, k, 2*k-1, 2*k-1) = -inv(G_22 / n);

    // Compute the covariance
    arma::mat cov = B * C * B.t();

    return cov / n;
  } catch(...) {
    Rcpp::stop("Cannot compute the two-step covariance!");
  }
}
