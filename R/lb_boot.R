#' @title Generic Bootstrap Interface
#'
#' @description
#' A unified interface for performing bootstrap resampling on different types of objects.
#' It automatically dispatches the appropriate method based on the input data type
#' (e.g., numeric vectors for simple statistics, or linear models for residual bootstrap).
#'
#' @param data The object to perform bootstrap on. Can be a numeric vector or an object of class \code{lm}.
#' @param ... Additional arguments passed to specific methods.
#'
#' @return An object of class \code{"lb_boot"} containing the original statistics and the bootstrap distribution.
#' @export
lb_boot <- function(data, ...) {
  UseMethod("lb_boot")
}

#' @describeIn lb_boot Bootstrap method for linear models (Residual Bootstrap)
#'
#' @param R Integer. The number of bootstrap replicates. Default is 1000.
#' @param type Character. "residual"(default) or "pairs".
#' @export
lb_boot.lm <- function(data, R = 1000, type = "residual", ...) {
  model <- data
  orig_coef <- coef(model)
  boot_coefs <- matrix(NA, nrow = R, ncol = length(orig_coef))
  colnames(boot_coefs) <- names(orig_coef)
  response_var <- all.vars(formula(model))[1]
  boot_data <- model$model

  if (type == "residual") {
    fit_vals <- fitted(model)
    res_vals <- resid(model)
    n <- length(res_vals)

    for (i in 1:R) {
      boot_res <- sample(res_vals, n, replace = TRUE)
      new_y <- fit_vals + boot_res
      boot_data[[response_var]] <- new_y
      boot_model <- update(model, data = boot_data)
      boot_coefs[i, ] <- coef(boot_model)
    }

  } else if (type == "pairs") {
    n <- nrow(boot_data)

    for (i in 1:R) {
      indices <- sample(1:n, n, replace = TRUE)
      boot_data_subset <- boot_data[indices, ]
      boot_model <- update(model, data = boot_data_subset)
      boot_coefs[i, ] <- coef(boot_model)
    }

  } else {
    stop("Invalid bootstrap type. Choose 'residual' or 'pairs'.")
  }

  result <- list(
    t0 = orig_coef,
    t = boot_coefs,
    R = R,
    type = type,
    call = match.call()
  )

  class(result) <- "lb_boot"
  return(result)
}

#' @describeIn lb_boot Bootstrap method for numeric vectors (Simple Sampling)
#'
#' @param stat A function to compute the statistic (e.g., mean, median, sd). Default is \code{mean}.
#' @export
lb_boot.numeric <- function(data, R = 1000, stat = mean, ...) {
  x <- data
  n <- length(x)

  # Calculate original statistic
  t0 <- stat(x)

  # Perform Bootstrap Resampling
  t_values <- replicate(R, {
    x_star <- sample(x, n, replace = TRUE)
    stat(x_star)
  })

  # Convert to matrix format to maintain consistency with other methods
  t_matrix <- matrix(t_values, ncol = 1)
  colnames(t_matrix) <- "statistic"

  # Package the results
  result <- list(
    t0 = t0,
    t = t_matrix,
    R = R,
    call = match.call()
  )

  class(result) <- "lb_boot"
  return(result)
}

#' Print method for lb_boot objects
#'
#' @description
#' Nicely formats the output of the lb_boot function, showing bias and standard error.
#'
#' @param x An object of class \code{lb_boot}.
#' @param ... Additional arguments (not used).
#' @export
print.lb_boot <- function(x, ...) {
  cat("\n=== LearnBootstrap Result ===\n")
  cat("Class:     ", class(x)[1], "\n")
  if (!is.null(x$type)) cat("Method:    ", x$type, "bootstrap\n")
  cat("Replicates:", x$R, "\n\n")

  cat("Statistics & Bootstrap Estimates:\n")
  se <- apply(x$t, 2, sd)
  bias <- colMeans(x$t) - x$t0
  tbl <- data.frame(Original = x$t0, Bias = bias, Std.Error = se)
  print(round(tbl, 4))
  cat("\n* Use confint() to calculate confidence intervals.\n")
  invisible(x)
}

#' Calculate Confidence Intervals for lb_boot Objects
#'
#' @description
#' Computes confidence intervals for bootstrap estimates using either the Percentile method or the Normal approximation.
#'
#' @param object An object of class \code{lb_boot}.
#' @param parm A specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level The confidence level required. Default is 0.95.
#' @param type Character string specifying the type of interval. Either "percentile" (default) or "normal".
#' @param ... Additional arguments (not used).
#'
#' @return A matrix containing the confidence intervals.
#' @export
confint.lb_boot <- function(object, parm, level = 0.95, type = "percentile", ...) {
  # 1. Prepare parameters
  # Extract original estimates and bootstrap distribution
  t0 <- object$t0
  t_dist <- object$t

  # Handle 'parm' argument (select specific parameters if requested)
  if (missing(parm)) {
    parm <- names(t0)
  } else if (is.numeric(parm)) {
    parm <- names(t0)[parm]
  }

  # Filter data for selected parameters
  # Check if t0 is a scalar (named vector of length 1) or vector
  if(length(t0) == 1) {
    # Special handling for single statistic (like mean) to keep matrix structure
    t_sub <- t_dist
    t0_sub <- t0
  } else {
    t_sub <- t_dist[, parm, drop = FALSE]
    t0_sub <- t0[parm]
  }

  alpha <- 1 - level
  probs <- c(alpha / 2, 1 - alpha / 2)

  # 2. Calculate Intervals based on Type
  if (type == "percentile") {
    # Method A: Percentile Method (Quantiles of the bootstrap distribution)
    ci <- t(apply(t_sub, 2, quantile, probs = probs, na.rm = TRUE))

  } else if (type == "normal") {
    # Method B: Normal Approximation (t0 +/- Z * SE)
    se <- apply(t_sub, 2, sd, na.rm = TRUE)
    z_score <- qnorm(1 - alpha / 2)

    lower <- t0_sub - z_score * se
    upper <- t0_sub + z_score * se
    ci <- cbind(lower, upper)
    colnames(ci) <- paste(format(probs * 100, trim = TRUE, scientific = FALSE, digits = 3), "%")

  } else {
    stop("Type must be 'percentile' or 'normal'.")
  }

  return(ci)
}


