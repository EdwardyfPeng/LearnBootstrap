#' @title Bootstrap Sampling from Several Populations
#' @description This is a function that can help you do bootstrap test of several populations using F statistics under same error variance.
#' @param data_list a list of data, each entry is data of one population
#' @param n_bootstrap the number of bootstrap samples that you want to construct
#'
#' @return F statistics based on the observed data, p value and F statistic based on the bootstrap samples
#' @export
#'
#' @examples
#' bootstrap_F_test(list(c1 = c(37.8, 37.1), c2 = c(37, 38.7,39.1), c3 = c(35.6, 32.4)),5000)
bootstrap_F_test <- function(data_list, n_bootstrap = 5000) {
  # Concatenate all groups to compute overall mean
  all_data <- unlist(data_list)
  overall_mean <- mean(all_data)

  # Calculate group means and errors
  group_means <- sapply(data_list, mean)
  errors <- lapply(data_list, function(group, group_mean) group - group_mean, group_mean = group_means)

  # Calculate the observed F statistic
  k <- length(data_list)
  n <- length(all_data)
  ni <- sapply(data_list, length)
  numerator_obs <- sum(ni * (group_means - overall_mean)^2)
  denominator_obs <- sum(unlist(lapply(errors, function(e) sum(e^2))))
  F_obs <- (numerator_obs / (k - 1)) / (denominator_obs / (n - k))

  # Perform bootstrap to generate F statistics
  bootstrap_F <- numeric(n_bootstrap)
  for (i in 1:n_bootstrap) {
    # Generate bootstrap samples by resampling errors within each group
    bootstrap_errors <- lapply(ni, function(size) sample(unlist(errors), size, replace = TRUE))
    bootstrap_group_means <- sapply(bootstrap_errors, mean)
    bootstrap_numerator <- sum(ni * (bootstrap_group_means - mean(bootstrap_group_means))^2)
    bootstrap_denominator <- sum(unlist(lapply(bootstrap_errors, function(e) sum(e^2))))
    bootstrap_F[i] <- (bootstrap_numerator / (k - 1)) / (bootstrap_denominator / (n - k))
  }

  # Calculate p-value
  p_value <- mean(bootstrap_F >= F_obs)

  # Return the results
  return(list(observed_F = F_obs, p_value = p_value, bootstrap_F_values = bootstrap_F))
}
