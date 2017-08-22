#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' @title Specification Function
//' @description G1
//' @param z Data
//' @param type Choice of the G1 function:
//' \itemize{
//'   \item 1: G1(z) = z
//'   \item 2: G1(z) = 0
//' }
//' @keywords internal
//' @export
// [[Rcpp::export]]
double G1_fun(double z, int type) {
  double out;
  if (type == 1) {
    out = z;
  } else if (type == 2) {
    out = 0;
  } else {
    Rcpp::stop("type not in 1, 2!");
  }
  return out;
}


//' @title Specification Function
//' @description G1_prime
//' @param z Data
//' @param type Choice of the G1_prime function:
//' \itemize{
//'   \item 1: G1_prime(z) = 1
//'   \item 2: G1_prime(z) = 0
//' }
//' @keywords internal
//' @export
// [[Rcpp::export]]
double G1_prime_fun(double z, int type) {
  double out;
  if (type == 1) {
    out = 1;
  } else if (type == 2) {
    out = 0;
  } else {
    Rcpp::stop("type not in 1, 2!");
  }
  return out;
}


//' @title Specification Function
//' @description G2_curly
//' @param z Data
//' @param type Choice of the G2_curly function:
//' \itemize{
//'   \item 1: -log(-z), z < 0
//'   \item 2: -sqrt(-z), z < 0
//'   \item 3: -1/z, z < 0
//'   \item 4: log(1 + exp(z))
//'   \item 5: exp(z)
//' }
//' @keywords internal
//' @export
// [[Rcpp::export]]
double G2_curly_fun(double z, int type) {
  double out;
  if ((type == 1) | (type == 2) | (type == 3)) {
    if (z >= 0) {
      Rcpp::warning("z can not be positive for type 1, 2, 3!");
      return NA_REAL;
    }
  }
  if (type == 1) {
    out = -log(-z);
  } else if (type == 2) {
    out = -sqrt(-z);
  } else if (type == 3) {
    out = -1/z;
  } else if (type == 4) {
    out = log1p(exp(z));
  } else if (type == 5) {
    out = exp(z);
  } else {
    Rcpp::stop("type not in 1, 2, 3, 4, 5!");
  }
  return out;
}


//' @title Specification Function
//' @description G2
//' @param z Data
//' @param type Choice of the G2 function:
//' \itemize{
//'   \item 1: -1/z, z < 0
//'   \item 2: 0.5/sqrt(-z), z < 0
//'   \item 3: 1/z^2, z < 0
//'   \item 4: 1 / (1 + exp(-z))
//'   \item 5: exp(z)
//' }
//' @keywords internal
//' @export
// [[Rcpp::export]]
double G2_fun(double z, int type) {
  double out;
  if ((type == 1) | (type == 2) | (type == 3)) {
    if (z >= 0) {
      Rcpp::warning("z can not be positive for type 1, 2, 3!");
      return NA_REAL;
    }
  }
  if (type == 1) {
    out = -1/z;
  } else if (type == 2) {
    out = 0.5/sqrt(-z);
  } else if (type == 3) {
    out = 1/pow(z, 2);
  } else if (type == 4) {
    out = 1/(1 + exp(-z));
  } else if (type == 5) {
    out = exp(z);
  } else {
    Rcpp::stop("type not in 1, 2, 3, 4, 5!");
  }
  return out;
}


//' @title Specification Function
//' @description G2_prime
//' @param z Data
//' @param type Choice of the G2_prime function:
//' \itemize{
//'   \item 1: 1/z^2, z < 0
//'   \item 2: 0.25 / (-z)^(3/2), z < 0
//'   \item 3: -2/z^3, z < 0
//'   \item 4: exp(z) / (1 + exp(z))^2
//'   \item 5: exp(z)
//' }
//' @keywords internal
//' @export
// [[Rcpp::export]]
double G2_prime_fun(double z, int type) {
  double out;
  if ((type == 1) | (type == 2) | (type == 3)) {
    if (z >= 0) {
      Rcpp::warning("z can not be positive for type 1, 2, 3!");
      return NA_REAL;
    }
  }
  if (type == 1) {
    out = 1/pow(z, 2);
  } else if (type == 2) {
    out = 0.25/pow(-z, 1.5);
  } else if (type == 3) {
    out = -2/pow(z, 3);
  } else if (type == 4) {
    out =  exp(z) / pow(1 + exp(z), 2);
  } else if (type == 5) {
    out = exp(z);
  } else {
    Rcpp::stop("type not in 1, 2, 3, 4, 5!");
  }
  return out;
}

//' @title Vectorized call to the G1 / G2 functions
//' @description Vectorized call to the G1 / G2 functions
//' @param z Vector
//' @param g String, either G1, G1_prime, G2_curly, G2 or G2_curly
//' @param type Numeric, for G1: 1-2; G2: 1-5
//' (see \link{G1_fun}, \link{G1_prime_fun}, \link{G2_curly_fun}, \link{G2_fun}, \link{G2_prime_fun})
//' @keywords internal
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector G_vec(Rcpp::NumericVector z, Rcpp::String g, int type) {
  int n = z.size();
  Rcpp::NumericVector out(n);
  if (g == "G1") {
    for (int i = 0; i < n; i++) out[i] = G1_fun(z[i], type);
  } else if (g == "G1_prime") {
    for (int i = 0; i < n; i++) out[i] = G1_prime_fun(z[i], type);
  } else if (g == "G2_curly") {
    for (int i = 0; i < n; i++) out[i] = G2_curly_fun(z[i], type);
  } else if (g == "G2") {
    for (int i = 0; i < n; i++) out[i] = G2_fun(z[i], type);
  } else if (g == "G2_prime") {
    for (int i = 0; i < n; i++) out[i] = G2_prime_fun(z[i], type);
  } else {
    Rcpp::stop("Non supported G-function!");
  }
  return out;
}


//' @title Joint (VaR, ES) loss for a linear predictor
//' @description Returns the loss for the parameter vector b
//' @param b Parameter vector
//' @param y Vector of dependent data
//' @param x Matrix of covariates. Note: intercept needs to be added manually
//' @param alpha Probability level
//' @param g1 1, 2 (see \link{G1_fun})
//' @param g2 1, 2, 3, 4, 5 (see \link{G2_curly_fun}, \link{G2_fun})
//' @param delta Smooth approximation of the indicator function (0 is the indicator function)
//' @importFrom Rcpp sourceCpp
//' @useDynLib esreg
//' @keywords internal
//' @export
// [[Rcpp::export]]
double esr_rho_lp(const arma::colvec& b, const arma::colvec& y, const arma::mat& x,
                  double alpha, int g1=2L, int g2=1L, double delta=0) {

  int n = x.n_rows;                           // Number of observations
  int k = x.n_cols;                           // Number of parameters
  arma::colvec bq = b.subvec(0, k-1);         // Quantile parameters
  arma::colvec be = b.subvec(k, 2*k-1);       // Expected shortfall parameters

  // Initialize variables
  double yi, xq, xe, h, out = 0;
  arma::mat xi;

  // Compute the loss
  for (int i = 0; i < n; i++) {
    yi = y(i);
    xi = x.row(i).t();
    xq = as_scalar(xi.t() * bq);
    xe = as_scalar(xi.t() * be);

    // Check the shortfall
    if (((g2 == 1) | (g2 == 2) | (g2 == 3)) & (xe >= 0)) {
      Rcpp::warning("x'b_e can not be positive for g2 1, 2, 3!");
      return NA_REAL;
    }

    // Hit variable or an approximation
    if (delta > 0) {h = 1 / (1 + exp(delta * (yi - xq)));} else {h = yi <= xq;}

    // Compute the loss
    out += (h - alpha) * G1_fun(xq, g1) - h * G1_fun(yi, g1) +
      G2_fun(xe, g2) * (xe - xq + (xq - yi) * h / alpha) - G2_curly_fun(xe, g2);
  }
  return out / n;
}


//' @title Identification (moment) function for the pair (VaR, ES) for a linear predictor
//' @description Returns the inner product psi' * psi of the moment function
//' for the parameter vector b
//' @param b Parameter vector
//' @param y Vector of dependent data
//' @param x Matrix of covariates. Note: intercept needs to be added manually
//' @param alpha Probability level
//' @param g1 1, 2 (see \link{G1_prime_fun})
//' @param g2 1, 2, 3, 4, 5 (see \link{G2_fun}, \link{G2_prime_fun})
//' @param delta Smooth approximation of the indicator function
//' (the default value 0 is the indicator function)
//' @importFrom Rcpp sourceCpp
//' @useDynLib esreg
//' @keywords internal
//' @export
// [[Rcpp::export]]
double esr_psi_lp(const arma::colvec& b, const arma::colvec& y, const arma::mat& x,
                  double alpha, int g1=2L, int g2=1, double delta=0) {

  int n = x.n_rows;                           // Number of observations
  int k = x.n_cols;                           // Number of parameters
  arma::colvec bq = b.subvec(0, k-1);         // Quantile parameters
  arma::colvec be = b.subvec(k, 2*k-1);       // Expected shortfall parameters

  // Initialize variables
  double yi, ui, xq, xe, h, g1p_xq, g2_xe, g2p_xe;
  arma::mat xi;
  arma::colvec psi = arma::zeros<arma::colvec>(2*k);

  // Loop over the observations
  for (int i = 0; i < n; i++) {
    yi = y(i);
    xi = x.row(i).t();
    xq = as_scalar(xi.t() * bq);
    xe = as_scalar(xi.t() * be);

    // Check the shortfall
    if (((g2 == 1) | (g2 == 2) | (g2 == 3)) & (xe >= 0)) {
      Rcpp::warning("x'b_e can not be positive for g2 1, 2, 3!");
      return NA_REAL;
    }

    // Hit variable or an approximation
    if (delta > 0) {h = 1 / (1 + exp(delta * (yi - xq)));} else {h = yi <= xq;}

    // Fill the vector
    psi.subvec(0, k-1)   += (h - alpha) / alpha * xi * (alpha * G1_prime_fun(xq, g1) + G2_fun(xe, g2));
    psi.subvec(k, 2*k-1) += xi * G2_prime_fun(xe, g2) * (xe - xq + (xq - yi) * h / alpha);
  }
  // Compute the inner product
  double out = as_scalar(psi.t() * psi) / (n*n);

  return out;
}
