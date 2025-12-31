#' Lattice-Ordered Alternatives Test
#'
#' @description
#' Performs a nonparametric test for lattice-ordered alternatives using the Jonckheere-Terpstra type statistic.
#' It tests the null hypothesis of no difference against an alternative where the treatment effects follow a specific matrix ordering.
#'
#' @param matrix_data A numeric matrix where rows represent levels of Factor A and columns represent levels of Factor B.
#'        The data should be arranged such that the hypothesized order increases with row and column indices.
#' @param n_permutations Integer. Number of permutations for the P-value calculation. Default is 1000.
#'
#' @return A list containing the observed test statistic (TL) and the permutation p-value.
#' @export
lb_lattice <- function(matrix_data, n_permutations = 1000) {

  data_mat <- as.matrix(matrix_data)
  rows <- nrow(data_mat)
  cols <- ncol(data_mat)

  calc_TL <- function(mat) {
    TL <- 0
    for (i in 1:rows) {
      for (j in 1:cols) {
        for (s in i:rows) {
          for (t in j:cols) {
            if (i == s && j == t) next
            if (mat[i, j] < mat[s, t]) {
              TL <- TL + 1
            }
          }
        }
      }
    }
    return(TL)
  }

  observed_TL <- calc_TL(data_mat)

  perm_stats <- numeric(n_permutations)
  flat_data <- as.vector(data_mat)
  total_len <- length(flat_data)

  for (k in 1:n_permutations) {
    shuffled_vec <- sample(flat_data, total_len)
    perm_mat <- matrix(shuffled_vec, nrow = rows, ncol = cols)
    perm_stats[k] <- calc_TL(perm_mat)
  }

  p_value <- mean(perm_stats >= observed_TL)

  result <- list(
    statistic = observed_TL,
    p_value = p_value,
    distribution = perm_stats
  )

  class(result) <- "lb_lattice"
  return(result)
}

#' Print method for lb_lattice
#' @export
print.lb_lattice <- function(x, ...) {
  cat("\n=== Lattice-Ordered Alternatives Test ===\n")
  cat("Statistic (TL): ", x$statistic, "\n")
  cat("P-value:        ", format(x$p_value, digits = 4), "\n")

  if (x$p_value < 0.05) {
    cat("Result:         Significant Trend (*)\n")
  } else {
    cat("Result:         No Significant Trend\n")
  }
  invisible(x)
}
