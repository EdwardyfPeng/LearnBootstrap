#' Unified Hypothesis Testing Interface
#'
#' @description
#' A unified interface for performing various bootstrap and permutation tests.
#' Currently supports Bootstrap ANOVA (equal/unequal variances) and Permutation Tests.
#'
#' @param x A formula specifying the model (e.g., \code{y ~ group} or \code{y ~ A * B}).
#' @param data A data frame containing the variables in the formula.
#' @param method The type of test to perform:
#' \itemize{
#'   \item \code{"anova"}: Bootstrap F-test (assuming equal variances, residual resampling).
#'   \item \code{"anova_unequal"}: Bootstrap F-test (assuming unequal variances, within-group resampling).
#'   \item \code{"perm"}: Permutation test for factor effects.
#' }
#' @param R Integer. Number of replicates. Default is 1000.
#' @param ... Additional arguments.
#'
#' @return An object of class \code{"lb_test"} containing the test statistic and p-value.
#' @export
lb_test <- function(x, ...) {
  UseMethod("lb_test")
}

#' @describeIn lb_test Formula method for hypothesis testing
#' @export
lb_test.formula <- function(x, data, method = "anova", R = 1000, ...) {
  formula <- x
  fit <- lm(formula, data = data)
  obs_anova <- anova(fit)
  f_obs <- obs_anova[1, "F value"]

  response_var <- all.vars(formula)[1] # Y
  group_var <- all.vars(formula)[2]    # Group
  boot_stats <- numeric(R)

  if (method == "anova") {
    y_grand_mean <- mean(data[[response_var]])
    res_vals <- resid(fit)
    n <- length(res_vals)

    for (i in 1:R) {
      boot_res <- sample(res_vals, n, replace = TRUE)
      new_y <- y_grand_mean + boot_res
      boot_data <- data
      boot_data[[response_var]] <- new_y
      boot_fit <- lm(formula, data = boot_data)
      boot_stats[i] <- anova(boot_fit)[1, "F value"]
    }

  } else if (method == "anova_unequal") {
    y_grand_mean <- mean(data[[response_var]])
    split_data <- split(data[[response_var]], data[[group_var]])
    centered_data_list <- lapply(split_data, function(g) g - mean(g))
    group_lengths <- sapply(split_data, length)
    group_names <- names(split_data)

    for (i in 1:R) {
      new_y_list <- lapply(centered_data_list, function(g_centered) {
        boot_resid <- sample(g_centered, length(g_centered), replace = TRUE)
        return(y_grand_mean + boot_resid)
      })

      new_df <- data.frame(
        y = unlist(new_y_list),
        grp = rep(group_names, group_lengths)
      )
      names(new_df) <- c(response_var, group_var)

      boot_fit <- lm(formula, data = new_df)
      boot_stats[i] <- anova(boot_fit)[1, "F value"]
    }

  } else if (method == "perm") {

    for (i in 1:R) {
      perm_data <- data
      perm_data[[group_var]] <- sample(perm_data[[group_var]])

      boot_fit <- lm(formula, data = perm_data)
      boot_stats[i] <- anova(boot_fit)[1, "F value"]
    }

  } else {
    stop("Unknown method. Choose 'anova', 'anova_unequal', or 'perm'.")
  }

  p_value <- mean(boot_stats >= f_obs)

  result <- list(
    statistic = f_obs,
    p_value = p_value,
    distribution = boot_stats,
    method = method,
    R = R,
    formula = formula
  )

  class(result) <- "lb_test"
  return(result)
}

#' Print method for lb_test
#' @export
print.lb_test <- function(x, ...) {
  cat("\n=== LearnBootstrap Hypothesis Test ===\n")

  method_name <- switch(x$method,
                        "anova" = "Bootstrap ANOVA (Equal Variance)",
                        "anova_unequal" = "Bootstrap ANOVA (Unequal Variance)",
                        "perm" = "Permutation Test")

  cat("Method:     ", method_name, "\n")
  cat("Replicates: ", x$R, "\n\n")

  cat("Observed F statistic: ", format(x$statistic, digits = 4), "\n")
  cat("P-value:              ", format(x$p_value, digits = 4), "\n")

  if (x$p_value < 0.05) {
    cat("Result:               Significant (*)\n")
  } else {
    cat("Result:               Not Significant\n")
  }
  invisible(x)
}
