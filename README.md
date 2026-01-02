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
