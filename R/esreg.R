#' @title Joint (VaR, ES) Regression
#' @description  Estimates a joint linear regression model for the pair (VaR, ES):
#' \deqn{Q_\alpha(Y | X) = X'\beta_q}
#' \deqn{ES_\alpha(Y | X) = X'\beta_e}
#' @param formula Formula object, e.g.: y ~ x1 + x2 + ...
#' @param data data.frame that holds the variables.
#' If missing the data is extracted from the environment.
#' @param alpha Probability level
#' @param g1 1, [2] (see \link{G1_fun}, \link{G1_prime_fun})
#' @param g2 [1], 2, 3, 4, 5 (see \link{G2_curly_fun}, \link{G2_fun}, \link{G2_prime_fun})
#' @param target The functions to be optimized: either the loss [rho] or the identification function (psi).
#' Optimization of the rho function is strongly recommended.
#' @param shift_data If g2 is 1, 2 or 3, we can either estimate the model without or with
#' shifting of the Y variable. We either risk biased estimates (no shifting) or slightly different estimates due
#' to the changed loss function (shifting). Defaults to shifting to avoid biased estimates.
#' @param method iterated local search [ils] or simulated annealing (sa)
#' @param control A list with control parameters passed to either the ils or sa:
#' \itemize{
#'   \item terminate_after: Stop the iterated local search if there is no improvement within max_step consecutive steps
#'   \item max.time: Maximum running time of the sa optimizer
#'   \item box: Box around the parameters for the sa optimizer
#' }
#' @return An esreg object
#' @seealso \code{\link{vcov.esreg}} for the covariance estimation and
#' \code{\link{summary.esreg}} for a summary of the regression results
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
#' # Estimate the model using the standard settings
#' fit <- esreg(y ~ x, alpha=alpha)
#'
#' # Compare the different variance-covariance estimators
#' cov1 <- vcov(object=fit, sparsity="iid", cond_var="ind")
#' cov2 <- vcov(object=fit, sparsity="nid", cond_var="scl_N")
#' cov3 <- vcov(object=fit, sparsity="nid", cond_var="scl_sp")
#'
#' print("Comparison of the variance-covariance estimators")
#' print(rbind(Truth=true_pars,
#'             Estimate=coef(fit),
#'             SE_iid_ind=sqrt(diag(cov1)),
#'             SE_ind_N=sqrt(diag(cov2)),
#'             SE_ind_sp=sqrt(diag(cov3))))
#'
#' # Compares estimates using different G2 functions
#' fit1 <- esreg(y ~ x, alpha=alpha, g2=1)
#' fit2 <- esreg(y ~ x, alpha=alpha, g2=2)
#' fit3 <- esreg(y ~ x, alpha=alpha, g2=3)
#' fit4 <- esreg(y ~ x, alpha=alpha, g2=4)
#' fit5 <- esreg(y ~ x, alpha=alpha, g2=5)
#' fits <- sapply(list(fit1, fit2, fit3, fit4, fit5), coef)
#' colnames(fits) <- sapply(1:5, function(i) esreg:::.G_function_names(1, i)[2])
#' print("Comparison of the five G2 functions")
#' print(rbind(Truth=true_pars, t(fits)))
#'
#' # Compare the M- and Z-estimator
#' fit_m <- esreg(y ~ x, alpha=alpha, target="rho")
#' fit_z <- esreg(y ~ x, alpha=alpha, target="psi")
#' print("Comparison of the M- and Z-estimator")
#' print(t(cbind(Truth=true_pars, M=coef(fit_m), Z=coef(fit_z))))
#' @references \href{https://arxiv.org/abs/1704.02213}{A Joint Quantile and Expected Shortfall Regression Framework}
#' @export
esreg <- function(formula, data, alpha, g1 = 2L, g2 = 1L, target = "rho", shift_data = TRUE,
                  method = "ils", control = list(terminate_after = 10, max.time = 10, box = 10)) {

  # Start the timer
  t0 <- Sys.time()

  # Check the inputs parameters
  if (!(g1 %in% c(1, 2))) stop("G1 can be 1 or 2.")
  if (!(g2 %in% c(1, 2, 3, 4, 5))) stop("G2 can be 1, 2, 3, 4 or 5.")
  if ((alpha < 0) | (1 < alpha)) stop("alpha not in (0, 1)")
  if (!(target %in% c("rho", "psi")))
    stop("target can be rho or psi.")
  if (!(method %in% c("ils", "sa")))
    stop("method can be ils or sa.")

  # Extract the formula
  if (missing(data)) {}
    data <- environment(formula)
  call <- match.call()
  mf <- stats::model.frame(formula = formula, data = data)
  x <- stats::model.matrix(attr(mf, "terms"), data = mf)
  y <- stats::model.response(mf)
  k <- ncol(x)

  # Check the input data
  if (any(is.na(y)) | any(is.na(x)))
    stop("Data contains NAs!")
  if (!(all(is.finite(y)) & all(is.finite(x))))
    stop("Not all values are finite!")
  if (length(y) != nrow(x))
    stop("Dimensions of y and x do not match!")

  # Store old random state and set a new seed for reproducibility
  if (!exists(".Random.seed", mode = "numeric", envir = globalenv()))
    sample(NA)
  oldSeed <- get(".Random.seed", mode = "numeric", envir = globalenv());
  set.seed(1)

  # Transform the data
  if ((shift_data == TRUE) & (g2 %in% c(1, 2, 3))) {
    max_y <- max(y)
    y <- y - max_y
  }

  # Set the target function
  if (target == "rho") {
    fun <- function(b) suppressWarnings(esr_rho_lp(b = b, y = y, x = x, alpha = alpha, g1 = g1, g2 = g2))
  } else if (target == "psi") {
    fun <- function(b) suppressWarnings(esr_psi_lp(b = b, y = y, x = x, alpha = alpha, g1 = g1, g2 = g2))
  }

  # Find starting values
  # Match quantile and expected shortfall levels
  e <- -stats::dnorm(stats::qnorm(alpha))/alpha
  alpha_tilde <- stats::uniroot(function(x) stats::qnorm(x) - e, c(0, alpha))$root
  # Fit the quantile regression
  fit_qr <- suppressWarnings(quantreg::rq(y ~ x - 1, tau = c(alpha, alpha_tilde)))
  # Get standard errors for the perturbation
  qr_sum <- suppressWarnings(summary(fit_qr, se="iid"))
  se <- c(qr_sum[[1]]$coefficients[,2], qr_sum[[2]]$coefficients[,2])
  # Store the coefficient estimates
  b0 <- as.vector(fit_qr$coef)

  # Optimize the model
  if (tolower(method) == "ils") { ## Iterated local search
    # Initial optimization
    fit <- try(stats::optim(par = b0, fn = fun, method = "Nelder-Mead"), silent = TRUE)

    # Counts the iterations without decrease of the loss
    counter <- 0
    while (counter < control$terminate_after) {
      # Perturbe b
      bt <- fit$par + stats::rnorm(2*k, sd=se)

      # If G2 is 1, 2, 3 and we do not shift the data, ensure that x'be < 0 by moving the es intercept
      if ((!shift_data) & (g2 %in% c(1, 2, 3))) {
        max_xe <- max(x %*% bt[(k+1):(2*k)])
        bt[k+1] <- bt[k+1] - (max_xe + 0.1) * (max_xe > 0)
      }

      # Fit the model with the perturbed parameters
      tmp_fit <- try(stats::optim(par = bt, fn = fun, method = "Nelder-Mead"), silent = TRUE)

      # Replace the fit if the new loss is smaller than the old. Otherwise increase the counter.
      if (!inherits(tmp_fit, "try-error")) {
        if (tmp_fit$value < fit$value) {
          fit <- tmp_fit
          counter <- 0
        } else {
          counter <- counter + 1
        }
      }
    }
  } else if (tolower(method) == "sa") {  ## Simmulated annealing
    if (!requireNamespace("GenSA", quietly = TRUE)) {
      stop("GenSA needed for this function to work. Please install it.", call. = FALSE)
    }
    fit <- GenSA::GenSA(par = b0, fn = fun,
                        lower = b0 - rep(control$box, 2 * k),
                        upper = b0 + rep(control$box, 2 * k),
                        control = list(max.time = control$max.time))
  }

  # Set names of the parameters
  b <- fit$par
  names(b) <- c(paste0("bq_", 1:k - 1), paste0("be_", 1:k - 1))

  # Undo the transformation
  if ((shift_data == TRUE) & (g2 %in% c(1, 2, 3))) {
    y <- y + max_y
    b[1] <- b[1] + max_y
    b[k + 1] <- b[k + 1] + max_y
    b0[1] <- b0[1] + max_y
    b0[k + 1] <- b0[k + 1] + max_y
  }

  # Reset the random seed to the old state
  assign(".Random.seed", oldSeed, envir=globalenv());

  # Compute the fitted values
  fitted.values <- data.frame(q = x %*% b[1:k],
                              e = x %*% b[(k + 1):(2 * k)])

  # Return results
  structure(list(call = call, formula = formula,
                 target = target, method = method, g1 = g1, g2 = g2, shift_data = shift_data,
                 alpha = alpha, y = y, x = x, b0 = b0,
                 coefficients = b,
                 coefficients_q = b[1:k],
                 coefficients_e = b[(k + 1):(2 * k)],
                 fitted.values = fitted.values,
                 value = fit$value,
                 time = Sys.time() - t0),
            class = "esreg")
}

#' @export
print.esreg <- function(x, digits = 4, ...) {
  cat("Call:\n")
  cat(deparse(x$call), "\n")
  cat("\nQuantile Coefficients:\n")
  print(format(x$coefficients_q, digits = digits), quote = FALSE)
  cat("\nExpected Shortfall Coefficients:\n")
  print(format(x$coefficients_e, digits = digits), quote = FALSE)
}

#' @export
fitted.esreg <- function(object, ...) {
  cbind(object$x %*% object$coefficients_q,
        object$x %*% object$coefficients_e)
}

#' @title Covariance of the joint (VaR, ES) estimator
#' @description Estimate the variance-covariance matrix of the joint (VaR, ES) estimator either using the asymptotic formulas or using the bootstrap.
#' @param object An esreg object
#' @param sparsity Sparsity estimator
#' \itemize{
#'   \item [iid] - Piecewise linear interpolation of the distribution
#'   \item nid - Hendricks and Koenker sandwich
#' }
#' @param cond_var Conditional truncated variance estimator
#' \itemize{
#'   \item [ind] Variance over all negative residuals
#'   \item scl_N Scaling with the normal distribution
#'   \item scl_sp Scaling with the kernel density function
#' }
#' @param bandwidth_type Bofinger, Chamberlain or Hall-Sheather
#' @param bootstrap_method
#' \itemize{
#'   \item [NULL] Use the asymptotic estimator
#'   \item iid bootstrap
#'   \item stationary bootstrap (Politis & Romano, 1994)
#' }
#' @param B Number of bootstrap iterations
#' @param block_length Average block length for the stationary bootstrap
#' @param ... additional arguments
#' @export
vcov.esreg <- function(object, sparsity = "iid", cond_var = "ind", bandwidth_type = "Hall-Sheather",
                       bootstrap_method = NULL, B = 1000, block_length = NULL, ...) {
  fit <- object

  if(is.null(bootstrap_method)) {
    if(!(sparsity %in% c("iid", "nid")))
      stop("sparsity can be iid or nid")
    if(!(cond_var %in% c("ind", "scl_N", "scl_sp")))
      stop("cond_var can be ind, scl_N or scl_sp")
    if(!(bandwidth_type %in% c("Bofinger", "Chamberlain", "Hall-Sheather")))
      stop("bandwidth_type can be Bofinger, Chamberlain or Hall-Sheather")

    # Extract some elements
    y <- fit$y
    x <- fit$x
    n <- nrow(x)
    k <- ncol(x)
    coefficients_q <- fit$coefficients_q
    coefficients_e <- fit$coefficients_e

    # Transform the data and coefficients
    if ((fit$shift_data == TRUE) & (fit$g2 %in% c(1, 2, 3))) {
      max_y <- max(y)
      y <- y - max_y
      coefficients_q[1] <- coefficients_q[1] - max_y
      coefficients_e[1] <- coefficients_e[1] - max_y
    }

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

    # Evaluate G1 / G2 functions
    G1_prime_xq <- G_vec(z = xq, g = "G1_prime", type = fit$g1)
    G2_xe <- G_vec(z = xe, g = "G2", type = fit$g2)
    G2_prime_xe <- G_vec(z = xe, g = "G2_prime", type = fit$g2)

    # Compute the covariance matrix
    cov <- l_esreg_covariance(x = x, xq = xq, xe = xe, alpha = fit$alpha,
                              G1_prime_xq = G1_prime_xq,
                              G2_xe = G2_xe, G2_prime_xe = G2_prime_xe,
                              density = dens, conditional_variance = cv)
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

    # Use the one_shot estimation approach for speed
    b <- apply(idx, 2, function(id) {
      fitb <- esreg(fit$y[id] ~ fit$x[id, -1],
                    alpha = fit$alpha, g1 = fit$g1, g2 = fit$g2,
                    method = "ils", control=list(terminate_after = 0))
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

#' @title esreg summary
#' @description Summarize details about the regression estimates.
#' @param object An esreg object
#' @param ... Accepts all parameters you can pass to \code{\link{vcov.esreg}}.
#' @export
summary.esreg <- function(object, ...) {
  fit <- object
  x <- fit$x
  n <- nrow(x)
  k <- ncol(x)
  cov <- vcov.esreg(fit, ...)
  se <- sqrt(diag(cov))
  tval <- stats::coef(fit) / se
  coef_mat <- cbind(Estimate     = stats::coef(fit),
                    `Std. Error` = se,
                    `t value`    = tval,
                    `Pr(>|t|)`   = 2 * stats::pt(abs(tval), n - 2*k, lower.tail = FALSE))

  structure(c(fit, list(cov = cov, coef_mat = coef_mat)), class = "summary.esreg")
}

#' @export
print.summary.esreg <- function(x, ...) {
  k <- length(x$coefficients_q)
  cat("Call:\n", paste0(deparse(x$call), sep = "\n", collapse = "\n"))
  cat("\nalpha: ", sprintf("%.3f", x$alpha), "\n")
  cat(.G_function_names(x$g1, x$g2)[1], "\n")
  cat(.G_function_names(x$g1, x$g2)[2], "\n")
  cat("Value: ", sprintf("%.9f", x$value), "\n")
  cat("\nQuantile Coefficients:\n")
  stats::printCoefmat(x$coef_mat[1:k,], signif.legend = FALSE)
  cat("\nExpected Shortfall Coefficients:\n")
  stats::printCoefmat(x$coef_mat[(k+1):(2*k),])
}

#' @export
residuals.esreg <- function(object, ...) {
  cbind(object$y - object$x %*% object$coefficients_q,
        object$y - object$x %*% object$coefficients_e)
}
