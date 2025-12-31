#' Aligned Rank Transform (ART) for Split-Plot Designs
#'
#' @description
#' Performs the Aligned Rank Transform test for interaction effects in split-plot designs.
#' This method aligns the data by removing main effects and then applying a rank-based ANOVA
#' on the residuals, as described in Higgins (Section 9.2).
#'
#' @param data A data frame containing the experimental data.
#' @param response Character. The name of the response variable column.
#' @param whole_plot Character. The name of the whole-plot factor.
#' @param sub_plot Character. The name of the sub-plot (within-block) factor.
#' @param block Character. The name of the block/subject identifier (e.g., "farm" or "plotID").
#'
#' @return An object of class \code{"anova"} containing the F-test results for the interaction.
#' @export
lb_art <- function(data, response, whole_plot, sub_plot, block) {
  y <- data[[response]]
  A <- data[[whole_plot]] # Whole plot factor
  B <- data[[sub_plot]]   # Sub plot factor
  Blk <- data[[block]]    # Block/Subject

  A <- as.factor(A)
  B <- as.factor(B)
  Blk <- as.factor(Blk)

  grand_mean <- mean(y, na.rm = TRUE)
  mu_ij <- aggregate(y, by = list(A, Blk), FUN = mean, na.rm = TRUE)
  colnames(mu_ij) <- c(whole_plot, block, "mean_ij")

  mu_k <- aggregate(y, by = list(B), FUN = mean, na.rm = TRUE)
  colnames(mu_k) <- c(sub_plot, "mean_k")

  merged_data <- merge(data, mu_ij, by = c(whole_plot, block))
  merged_data <- merge(merged_data, mu_k, by = sub_plot)

  merged_data$aligned_value <- merged_data[[response]] -
    merged_data$mean_ij -
    merged_data$mean_k +
    grand_mean

  merged_data$aligned_rank <- rank(merged_data$aligned_value)

  f_str <- paste("aligned_rank ~", whole_plot, "*", sub_plot,
                 "+ Error(", block, "/", whole_plot, ")")
  fit_art <- aov(as.formula(f_str), data = merged_data)

  return(summary(fit_art))
}
