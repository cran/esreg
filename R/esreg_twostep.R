#' @title Two Step Quantile and Expected Shortfall Regression
#' @description Estimates the expected shortfall in two steps. First a linear quantile regression,
#' then a weighted least squares regression (see the Oracle estimator in the references for a brief explanation).
#' This estimator is much faster than the joint estimator \code{\link{esreg}}.
#' However, the estimates are often less precise and it is recommended primarily if one needs
#' estimates quickly.
#' @param formula Formula object, e.g.: y ~ x1 + x2 + ...
#' @param data data.frame that holds the variables. Can be missing.
#' @param alpha Probability level
#' @examples
#' # Simulate data (DGP-(2) in the linked paper)
#' set.seed(0)
#' x <- rchisq(2000, df=1)
#' y <- -x + (1 + 0.5 * x) * rnorm(1000)
#'
#' # True quantile and expected shortfall regression parameters (for alpha=0.025)
#' alpha=0.025
#' true_pars <- c(-1.959964, -1.979982, -2.337803, -2.168901)
#'
#' # Joint estimator
#' fit_joint <- esreg(y ~ x, alpha=alpha)
#'
#' # Two-step estimator
#' fit_twostep <- esreg_twostep(y ~ x, alpha=alpha)
#'
#' # Compare the estimates
#' print(rbind(Truth=true_pars, Joint=coef(fit_joint), `Two-Step`=coef(fit_twostep)))
#'
#' # ... and the estimation times
#' print(c(Joint=fit_joint$time, `Two Step`=fit_twostep$time))
#' @seealso \code{\link{vcov.esreg_twostep}} for the covariance estimation and
#' \code{\link{summary.esreg_twostep}} for a summary of the regression results
#' @references \href{https://arxiv.org/abs/1704.02213}{A Joint Quantile and Expected Shortfall Regression Framework}
#' @export
esreg_twostep <- function(formula, data, alpha) {

  # Start the timer
  t0 <- Sys.time()

  # Extract the formula
  if (missing(data))
    data <- environment(formula)
  call <- match.call()
  mf <- stats::model.frame(formula = formula, data = data)
  x <- stats::model.matrix(attr(mf, "terms"), data = mf)
  y <- stats::model.response(mf)

  # Check the data
  if (any(is.na(y)) | any(is.na(x)))
    stop("Data contains NAs!")
  if (!(all(is.finite(y)) & all(is.finite(x))))
    stop("Not all values are finite!")

  # First stage: quantile regression
  fit_rq <- quantreg::rq(y ~ x-1, tau = alpha)
  q <- fit_rq$fitted.values
  b_q <- fit_rq$coefficients

  # Second stage: weighted least squares regression
  fit <- stats::lm(y ~ x-1, weights = (y <= q) * 1)
  b_e <- fit$coefficients

  # Name the coefficents
  names(b_q) <- paste0("bq_", 1:length(b_q) - 1)
  names(b_e) <- paste0("be_", 1:length(b_e) - 1)

  # Return results
  structure(list(call = call, formula = formula,
                 alpha = alpha, y = y, x = x,
                 coefficients = c(b_q, b_e),
                 coefficients_q = b_q,
                 coefficients_e = b_e,
                 time = Sys.time() - t0),
            class = "esreg_twostep")
}

#' @export
print.esreg_twostep <- print.esreg

#' @export
fitted.esreg_twostep <- fitted.esreg

#' @title Covariance Estimation for esreg_twostep
#' @description Estimate the variance-covariance matrix of the joint (VaR, ES) estimator
#' either using the asymptotic formulas or using the bootstrap.
#' @inheritParams vcov.esreg
#' @export
vcov.esreg_twostep <- function(object, sparsity = "iid", cond_var = "ind", bandwidth_type = "Hall-Sheather",
                       bootstrap_method = NULL, B = 1000, block_length = NULL, ...) {
  fit <- object

  if(is.null(bootstrap_method)) {
    if(!(sparsity %in% c("iid", "nid")))
      stop("sparsity can be iid or nid")
    if(!(cond_var %in% c('ind', 'scl_N', 'scl_t')))
      stop('cond_var can be ind, scl_N or scl_t')
    if(!(bandwidth_type %in% c("Bofinger", "Chamberlain", "Hall-Sheather")))
      stop("bandwidth_type can be Bofinger, Chamberlain or Hall-Sheather")

    # Extract some elements
    y <- fit$y
    x <- fit$x
    n <- nrow(x)
    k <- ncol(x)
    coefficients_q <- fit$coefficients_q
    coefficients_e <- fit$coefficients_e

    # Precompute some quantities
    xq <- as.numeric(x %*% coefficients_q)
    xe <- as.numeric(x %*% coefficients_e)
    u <- as.numeric(y - xq)

    # Check the methods in case of sample quantile / es
    if ((k == 1) & sparsity != "iid") {
      warning("Changed sparsity estimation to iid!")
      sparsity <- "iid"
    }
    if ((k == 1) & cond_var != "ind") {
      warning("Changed conditional truncated variance estimation to nid!")
      cond_var <- "ind"
    }

    # Density quantile function
    dens <- density_quantile_function(y = y, x = x, u = u, alpha = fit$alpha,
                                      sparsity = sparsity, bandwidth_type = bandwidth_type)

    # Truncated conditional variance
    cv <- conditional_truncated_variance(y = u, x = x, approach = cond_var)

    # Compute and return the covariance matrix
    cov <- l_esreg_twostep_covariance(x = x, xq = xq, xe = xe, alpha = fit$alpha,
                                      density = dens, conditional_variance = cv)

    # Flip blockwise such that it matches the rest of the esreg package
    cov <- rbind(cbind(cov[(k+1):(2*k), (k+1):(2*k)], cov[1:k, (k+1):(2*k)]),
                 cbind(cov[(k+1):(2*k), 1:k], cov[1:k, 1:k]))
  } else {
    if (!(bootstrap_method %in% c(NULL, "iid", "stationary")))
      stop("bootstrap_method can be NULL, iid or stationary")
    if (B < 1000)
      warning("The number of bootstrap iterations is small!")

    # Draw the bootstrap indices
    n <- length(fit$y)
    set.seed(1)
    if (bootstrap_method == "iid") {
      idx <- matrix(sample(1:n, size = n * B, replace = TRUE), nrow = n)
    } else if (bootstrap_method == "stationary") {
      if (is.null(block_length))
        stop("No average block length provided!")
      idx <- stationary_bootstrap_indices(n = n, avg_block_size = block_length, B = B)
    }

    b <- apply(idx, 2, function(id) {
      fitb <- esreg_twostep(fit$y[id] ~ fit$x[id, -1], alpha = fit$alpha)
      fitb$coefficients
    })

    # Compute the covariance
    cov <- stats::cov(t(b))
  }

  # Set the names
  rownames(cov) <- colnames(cov) <- names(stats::coef(fit))

  # Return the estimated covariance
  cov
}

#' @title esreg_twostep summary
#' @description Summarize details about the regression estimates.
#' @param object An esreg_twostep object
#' @param ... Accepts all parameters you can pass to \code{\link{vcov.esreg_twostep}}.
#' @keywords internal
#' @export
summary.esreg_twostep <- function(object, ...) {
  fit <- object
  x <- fit$x
  n <- nrow(x)
  k <- ncol(x)
  cov <- vcov.esreg_twostep(fit, ...)
  se <- sqrt(diag(cov))
  tval <- stats::coef(fit) / se
  coef_mat <- cbind(Estimate     = stats::coef(fit),
                    `Std. Error` = se,
                    `t value`    = tval,
                    `Pr(>|t|)`   = 2 * stats::pt(abs(tval), n - 2*k, lower.tail = FALSE))

  structure(c(fit, list(cov = cov, coef_mat = coef_mat)), class = "summary.esreg_twostep")
}

#' @export
print.summary.esreg_twostep <- function(x, ...) {
  k <- length(x$coefficients_q)
  cat("Call:\n", paste0(deparse(x$call), sep = "\n", collapse = "\n"))
  cat("\nalpha: ", sprintf("%.3f", x$alpha), "\n")
  cat("\nQuantile Coefficients:\n")
  stats::printCoefmat(x$coef_mat[1:k,], signif.legend = FALSE)
  cat("\nExpected Shortfall Coefficients:\n")
  stats::printCoefmat(x$coef_mat[(k+1):(2*k),])
}

#' @export
residuals.esreg_twostep <- function(object, ...) {
  cbind(object$y - object$x %*% object$coefficients_q,
        object$y - object$x %*% object$coefficients_e)
}
