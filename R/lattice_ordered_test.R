#' @title Function to perform the lattice-ordered alternatives test
#'
#' @param data the data you want to perform the test
#' @param A_levels A level
#' @param B_levels B level
#' @param n_permutations permutation time
#'
#' @return observed_TL, p value
#' @export
#'
#' @examples
#' lattice_ordered_test(matrix(c(1,2,3,4,5,6,7,8,9),nrow = 3,byrow = TRUE,dimnames=list(c("P1", "P2", "P3"), c("Week 1", "Week 2", "Week 4"))), c("P1", "P2", "P3"), c("Week 1", "Week 2", "Week 4"))
lattice_ordered_test <- function(data, A_levels, B_levels, n_permutations = 1000) {
  # Ensure the data is in matrix form
  data <- as.matrix(data)

  # Number of rows (levels of A) and columns (levels of B)
  a <- length(A_levels)
  b <- length(B_levels)

  # Calculate the observed test statistic T_L
  observed_TL <- 0
  for (i in 1:(a - 1)) {
    for (j in 1:(b - 1)) {
      for (s in (i + 1):a) {
        for (t in (j + 1):b) {
          observed_TL <- observed_TL + (data[i, j] < data[s, t])
        }
      }
    }
  }

  # Permutation test to obtain the distribution of T_L under the null hypothesis
  permuted_TL <- numeric(n_permutations)
  for (k in 1:n_permutations) {
    permuted_data <- data[sample(a * b)]
    dim(permuted_data) <- c(a, b)
    permuted_TL[k] <- 0
    for (i in 1:(a - 1)) {
      for (j in 1:(b - 1)) {
        for (s in (i + 1):a) {
          for (t in (j + 1):b) {
            permuted_TL[k] <- permuted_TL[k] + (permuted_data[i, j] < permuted_data[s, t])
          }
        }
      }
    }
  }

  # Calculate p-value
  p_value <- mean(permuted_TL >= observed_TL)

  # Return the test statistic and p-value
  return(list(observed_TL = observed_TL, p_value = p_value))
}
