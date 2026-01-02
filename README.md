# LearnBootstrap: Nonparametric Inference with Resampling Methods

## 📦 Installation

You can install the development version of `LearnBootstrap` from GitHub with:
```{r}
# install.packages("devtools")
devtools::install_github("EdwardyfPeng/LearnBootstrap")
```

## 🚀 Key features

* **Bootstrapped Estimation**:
    * Estimates standard errors and confidence intervals (Normal approximation & Percentile) for:
        * Sample statistics (e.g., standard deviation).
        * Regression coefficients (via Residual Bootstrap and Pairwise Bootstrap).

* **Bootstrap ANOVA**:
    * **Residual Bootstrap**: For one-way ANOVA under the assumption of equal variances (Homoscedasticity).
    * **Stratified Bootstrap**: For one-way ANOVA when variances are unequal (Heteroscedasticity).

* **Permutation Tests**:
    * Implements Exact or Monte Carlo permutation tests for assessing factor effects (e.g., two-sample comparisons or one-way designs).

* **Aligned Rank Transform (ART)**:
    * A nonparametric rank test specifically designed for detecting **interaction effects** in Split-Plot designs.

* **Lattice-Ordered Tests**:
    * A specialized test for detecting specific monotonic trends (e.g., dose-response surfaces) in two-factor experiments.
 
## 📚 Examples
1. Bootstrapping Statistics (Confidence Intervals)
Standard formulas for the standard error of the Median or CV (Coefficient of Variation) are complex or non-existent. `lb_boot()` makes this easy.

```{R}
library(LearnBootstrap)
set.seed(123)

# Generate skewed data (e.g., income distribution)
income_data <- rlnorm(50, meanlog = 2, sdlog = 1)

# 1. Estimate the Median and its 95% CI
# Note: Traditional t-tests cannot do this for the median!
boot_median <- lb_boot(income_data, stat = median, R = 2000)
print(boot_median)

# 2. Estimate a custom statistic: Coefficient of Variation (CV)
# Define the function: CV = SD / Mean
calc_cv <- function(x) { sd(x) / mean(x) }

boot_cv <- lb_boot(income_data, stat = calc_cv, R = 2000)
print(boot_cv)
```

2. Robust ANOVA (Handling Unequal Variances)
When groups have different variances (Heteroscedasticity), standard ANOVA fails. Use Stratified Bootstrap to fix this.

```{r}
# Create data with unequal variances
# Group A: small variance, Group C: huge variance
df_anova <- data.frame(
  group = factor(rep(c("A", "B", "C"), each = 15)),
  value = c(rnorm(15, 10, 1),   # SD = 1
            rnorm(15, 12, 2),   # SD = 2
            rnorm(15, 15, 5))   # SD = 5 (Heteroscedasticity!)
)

# A. Standard Residual Bootstrap (Assumes Equal Variance)
# This might yield incorrect Type I error rates here.
res_boot <- lb_boot(value ~ group, data = df_anova, test = "residual")
print(res_boot)

# B. Stratified Bootstrap (Robust)
# Resamples within each group separately, preserving the variance structure.
strat_boot <- lb_boot(value ~ group, data = df_anova, test = "stratified")
print(strat_boot)
```
3. Bootstrap for Linear Regression
In linear regression, standard p-values rely on the assumption that residuals are normally distributed. Residual Bootstrap allows you to estimate standard errors and confidence intervals for regression coefficients without this assumption.

```{R}
# Generate synthetic regression data with non-normal noise
set.seed(42)
n <- 30
x <- runif(n, 0, 10)
# True relation: y = 2 + 3x + noise (exponential noise, right-skewed)
y <- 2 + 3 * x + rexp(n) 

# Fit a linear model
model <- lm(y ~ x)

# Run Residual Bootstrap to estimate coefficients
# R = 1000 repetitions
boot_reg <- lb_boot(model, R = 1000)

print(boot_reg)

# Optional: Calculate 95% Confidence Intervals
# Compare these with confint(model) which assumes normality
# confint(boot_reg) # If you have implemented a confint method

4. Aligned Rank Transform (ART) for Split-Plot
Designed for Interaction Effects in split-plot designs where normality assumptions are violated.

Scenario: An experiment on dogs (Block) testing the effect of Diabetes (Whole Plot) and Drug Delivery Method (Sub Plot) on insulin turnover.

```{R}
# Construct the dataset
dogs_data <- data.frame(
  dog_id   = factor(rep(1:10, each=2)), 
  diabetes = factor(rep(c("Diabetic", "Healthy"), each=10)), 
  method   = factor(rep(c("Inject", "Infuse"), 10)),
  turnover = c(44,28, 33,23, 38,34, 59,19, 46,26,  # Diabetic group
               54,42, 43,23, 55,23, 71,27, 57,35)  # Healthy group
)

# Run ART
# We are specifically looking for the 'diabetes:method' interaction
art_result <- lb_art(
  data       = dogs_data, 
  response   = "turnover", 
  whole_plot = "diabetes", 
  sub_plot   = "method", 
  block      = "dog_id"
)

print(art_result)
```

5. Lattice-Ordered Test (Detecting Matrix Trends)
Detects if values increase across both rows and columns simultaneously (e.g., higher dose combination = higher toxicity).

```{R}
# A 3x3 dose-response matrix
# Hypothesized trend: Values increase from Top-Left to Bottom-Right
toxicity_matrix <- matrix(
  c(5,  10, 15,
    12, 20, 28,
    25, 35, 50), 
  nrow = 3, 
  byrow = TRUE
)

# Run the permutation test for lattice ordering
lattice_res <- lb_lattice(toxicity_matrix, n_permutations = 5000)

print(lattice_res)
```
