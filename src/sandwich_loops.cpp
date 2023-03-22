#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' @keywords internal
// [[Rcpp::export]]
arma::mat lambda_matrix_loop(
    arma::mat xq, arma::mat xe, arma::vec xbq, arma::vec xbe,
    arma::vec G1_prime_xq, arma::vec G1_prime_prime_xq,
    arma::vec G2_xe, arma::vec G2_prime_xe, arma::vec G2_prime_prime_xe,
    arma::vec density, arma::vec cdf,
    double alpha, bool include_misspecification_terms) {
  int n = xq.n_rows;
  int kq = xq.n_cols;
  int ke = xe.n_cols;

  // Define some 0-matrices
  arma::mat lambda_11 = arma::zeros<arma::mat>(kq, kq);
  arma::mat lambda_12 = arma::zeros<arma::mat>(kq, ke);
  arma::mat lambda_22 = arma::zeros<arma::mat>(ke, ke);
  arma::mat lambda = arma::zeros<arma::mat>(kq+ke, kq+ke);

  // Compute the matrix elements
  for (int i = 0; i < n; i++) {
    arma::mat xqi = xq.row(i);
    arma::mat xei = xe.row(i);
    arma::mat xxq = xqi.t() * xqi;
    arma::mat xxe = xei.t() * xei;
    arma::mat xxqe = xqi.t() * xei;
    double xbqi = xbq(i);
    double cdfi = cdf(i);

    if (include_misspecification_terms) {
      lambda_11 += xxq * (G1_prime_xq(i) + G2_xe(i) / alpha) * density(i) + xxq * G1_prime_prime_xq(i) * (cdfi - alpha);
      lambda_12 += xxqe * G2_prime_xe(i) * (cdfi - alpha) / alpha;
      lambda_22 += xxe * G2_prime_xe(i) + xxe * G2_prime_prime_xe(i) * xbqi * (cdfi - alpha) / alpha;
    } else {
      lambda_11 += xxq * (G1_prime_xq(i) + G2_xe(i) / alpha) * density(i);
      lambda_22 += xxe * G2_prime_xe(i);
    }
  }

  // Fill the matrices
  lambda.submat(0, 0, kq-1, kq-1) = lambda_11 / n;
  lambda.submat(0, kq, kq-1, kq+ke-1) = lambda_12 / n;
  lambda.submat(kq, 0, kq+ke-1, kq-1) = lambda_12.t() / n;
  lambda.submat(kq, kq, kq+ke-1, kq+ke-1) = lambda_22 / n;

  return lambda;
}

//' @keywords internal
// [[Rcpp::export]]
arma::mat sigma_matrix_loop(
    arma::mat xq, arma::mat xe, arma::colvec xbq, arma::colvec xbe,
    arma::colvec G1_prime_xq, arma::colvec G2_xe, arma::colvec G2_prime_xe,
    arma::colvec conditional_variance, arma::vec cdf,
    double alpha, bool include_misspecification_terms) {
  int n = xq.n_rows;
  int kq = xq.n_cols;
  int ke = xe.n_cols;

  // Define some 0-matrices
  arma::mat sigma_11 = arma::zeros<arma::mat>(kq, kq);
  arma::mat sigma_12 = arma::zeros<arma::mat>(ke, kq);
  arma::mat sigma_22 = arma::zeros<arma::mat>(ke, ke);
  arma::mat sigma = arma::zeros<arma::mat>(kq+ke, kq+ke);

  // Compute the matrix elements
  for (int i = 0; i < n; i++) {
    arma::mat xqi = xq.row(i);
    arma::mat xei = xe.row(i);
    arma::mat xxq = xqi.t() * xqi;
    arma::mat xxe = xei.t() * xei;
    arma::mat xxeq = xei.t() * xqi;
    double xbqi = xbq(i);
    double xbei = xbe(i);
    double cdfi = cdf(i);

    arma::mat _sigma_11 = xxq * pow(alpha*G1_prime_xq(i) + G2_xe(i), 2);
    arma::mat _sigma_12 = xxeq * (alpha*G1_prime_xq(i) + G2_xe(i)) * G2_prime_xe(i);
    arma::mat _sigma_22 = xxe * pow(G2_prime_xe(i), 2);

    if (include_misspecification_terms) {
      sigma_11 +=
        _sigma_11 * (1-alpha)/alpha +
        _sigma_11 * (1 - 2*alpha) * (cdfi - alpha) / pow(alpha, 2);
      sigma_12 +=
        _sigma_12 * (1-alpha)/alpha * (xbqi - xbei) +
        _sigma_12 * ((1-alpha)/alpha * xbqi * (cdfi-alpha)/alpha - (cdfi-alpha)/alpha * (xbqi - xbei));
      sigma_22 +=
        _sigma_22 * (conditional_variance(i)/alpha + (1-alpha)/alpha * pow(xbqi - xbei, 2)) +
        _sigma_22 * 2 * (xbqi - xbei) * xbqi * (alpha - cdfi) / alpha;
    } else {
      sigma_11 += _sigma_11 *
        (1-alpha)/alpha;
      sigma_12 += _sigma_12 *
        (1-alpha)/alpha * (xbqi - xbei);
      sigma_22 += _sigma_22 *
        (conditional_variance(i)/alpha + (1-alpha)/alpha * pow(xbqi - xbei, 2));
    }
  }

  // Fill the matrices
  sigma.submat(0, 0, kq-1, kq-1) = sigma_11 / n;
  sigma.submat(kq, 0, kq+ke-1, kq-1) = sigma_12 / n;
  sigma.submat(0, kq, kq-1, kq+ke-1) = sigma_12.t() / n;
  sigma.submat(kq, kq, kq+ke-1, kq+ke-1) = sigma_22 / n;

  return sigma;
}



//' @keywords internal
// [[Rcpp::export]]
arma::mat estimating_function_loop(
    arma::vec y, arma::mat xq, arma::mat xe, arma::colvec xbq, arma::colvec xbe,
    arma::colvec G1_prime_xq, arma::colvec G2_xe, arma::colvec G2_prime_xe, double alpha) {
  int n = xq.n_rows;
  int kq = xq.n_cols;
  int ke = xe.n_cols;

  // Define some 0-matrices
  arma::mat psi = arma::zeros<arma::mat>(n, kq + ke);

  // Compute the matrix elements
  for (int i = 0; i < n; i++) {
    double yi = y(i);
    arma::mat xqi = xq.row(i);
    arma::mat xei = xe.row(i);
    arma::mat xxq = xqi.t() * xqi;
    arma::mat xxe = xei.t() * xei;
    double xbqi = xbq(i);
    double xbei = xbe(i);

    // Hit variable
    bool hit = yi <= xbqi;

    // Fill the matrix
    psi.submat(i, 0, i, kq-1) = xqi * (G1_prime_xq(i) + G2_xe(i)/alpha) * (hit - alpha);
    psi.submat(i, kq, i, kq+ke-1) = xei *  G2_prime_xe(i) * (xbei - xbqi + (xbqi - yi) * hit / alpha);
  }

  return psi;
}
