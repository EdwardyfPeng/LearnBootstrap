#' Bootstrap法估计均方误差（MSE）
#'
#' 这个函数使用Bootstrap法来估计一个统计函数应用于给定数据集时的均方误差（MSE）。
#'
#' @param x 一个数值向量，包含原始数据。
#' @param fun 一个函数，其MSE将被估计。该函数应接受一个数值向量作为输入，并返回一个数值。
#' @param rep 一个整数，指定自助法过程的重复次数。
#'
#' @return 一个数值，代表函数\code{fun}应用于数据\code{x}时的估计均方误差（MSE）。
#'
#' @export
#'
#' @examples
#' # 使用均值函数的例子
#' x <- rnorm(100)
#' mean_mse <- bootest(x, mean, 1000)
#' mean_mse
#'
bootest_mse <- function(x, fun, rep) {
  theta <- fun(x)                        # 从原始数据中估计感兴趣的统计量。
  stats <- replicate(rep, fun(sample(x, replace = TRUE))) # 执行自助法过程，通过有放回地重采样并估计统计量。
  mse <- mean((stats - theta)^2)         # 计算自助法估计与原始估计之间的均方误差（MSE）。
  return(mse)                            # 返回估计的MSE。
}



#' 计算变异系数
#'
#' 这个函数计算输入向量的变异系数，它是标准差与均值的比例，乘以100转换为百分比形式。
#'
#' @param x 一个数值向量，其变异系数将被计算。
#'
#' @return 变异系数的值，以百分比形式表示。
#'
#' @export
#'
#' @examples
#' x <- c(7, 11, 15, 16, 20, 22, 24, 25, 29, 33)
#' cv_result <- cv(x)
#' print(cv_result)
cv <- function(x) {
  100 * sqrt(var(x)) / mean(x)
}



#' Bootstrap方法估计偏差与方差
#'
#' 这个函数使用Bootstrap方法来估计一个统计函数应用于给定数据集时的偏差和方差。
#'
#' @param x 一个数值向量，包含原始数据。
#' @param fun 一个函数，其偏差和方差将被估计。该函数应接受一个数值向量作为输入，并返回一个数值。
#' @param rep 一个整数，指定Bootstrap过程的重复次数。
#'
#' @return 一个数值向量，包含估计的偏差和方差。
#'
#' @export
#'
#' @examples
#' x <- c(7, 11, 15, 16, 20, 22, 24, 25, 29, 33)
#' bootest_result <- bootest_bias_var(x, mean, 1000)
#' print(bootest_result)
bootest_bias_var <- function(x, fun, rep) {
  theta <- fun(x)                        # 计算原始数据集的统计量。
  stats <- replicate(rep, fun(sample(x, replace = TRUE))) # 执行Bootstrap过程，重复计算统计量。
  mse <- mean((stats - theta)^2)         # 计算均方误差（MSE）。
  bias <- mean(stats) - theta            # 计算偏差，即Bootstrap统计量的平均值与原始统计量的差。
  var <- var(stats)                      # 计算方差，即Bootstrap统计量的方差。
  c(bias, var)                           # 返回一个包含偏差和方差的数值向量。
}



#' 计算样本统计量与已知参数之间的差异的标准化值
#'
#' 这个函数用于计算一个给定样本的统计量与一个已知参数之间的差异的标准化值。
#' 这个标准化值可以帮助我们理解样本统计量与已知参数值之间的差异，并将其标准化到与样本大小无关的量度。
#' 在假设检验中，这个标准化值可以用于计算统计检验的统计量，如t统计量或z统计量。
#'
#' @param x 一个数值向量，包含样本数据。
#' @param theta 一个数值，代表已知的参数值。
#'
#' @return 一个数值，表示样本统计量与已知参数之间的差异的标准化值。
#'
#' @export
#'
#' @examples
#' x <- c(7, 11, 15, 16, 20, 22, 24, 25, 29, 33)
#' theta <- 18
#' stt_result <- stt(x, theta)
#' print(stt_result)
stt <- function(x, theta) {
  mx <- mean(x)                      # 计算样本均值
  n <- length(x)                     # 计算样本长度
  s <- sqrt(var(x))                  # 计算样本方差，然后开平方得到标准差
  sqrt(n) * (mx - theta) / s         # 计算样本统计量与已知参数之间的差异的标准化值
}

#' 计算卡方统计量
#'
#' 这个函数用于计算给定样本的卡方统计量。
#' 卡方统计量是一种用于比较样本方差与假设方差（s2）的统计量。
#'
#' @param x 一个数值向量，包含样本数据。
#' @param s2 一个数值，代表假设方差。
#'
#' @return 一个数值，表示卡方统计量。
#'
#' @export
#'
#' @examples
#' x <- c(7, 11, 15, 16, 20, 22, 24, 25, 29, 33)
#' s2 <- 100
#' chi_result <- chi(x, s2)
#' print(chi_result)
chi <- function(x, s2) {
  n <- length(x)                     # 计算样本长度
  chi <- (n - 1) * var(x)            # 计算样本方差的卡方统计量
  chi / s2                           # 将卡方统计量除以假设方差，得到卡方统计量的值
}


#' 计算样本的Aest统计量
#'
#' 这个函数用于计算给定样本的Aest统计量，这是一种用于估计样本数据中潜在异常值的统计量。
#' Aest统计量可以用于识别和处理数据集中的异常值。
#'
#' @param x 一个数值向量，包含样本数据。
#' @param func 一个函数，用于计算每个子集的统计量。这个函数应该接受一个数值向量作为输入，并返回一个数值。
#'
#' @return 一个数值，表示Aest统计量。
#'
#' @export
#'
#' @examples
#' x <- c(7, 11, 15, 16, 20, 22, 24, 25, 29, 33)
#' func <- function(subset) mean(subset)
#' aest_result <- aest(x, func)
#' print(aest_result)
aest <- function(x, func) {
  n <- length(x)                     # 计算样本长度
  tmp <- rep(0, n)                   # 创建一个与x长度相同的向量，初始化为0
  for(i in 1:n) {                   # 遍历样本的每个元素
    tmp[i] <- func(x[-i])            # 调用func函数，去除当前元素后的子集
  }
  err <- mean(tmp) - tmp             # 计算子集统计量与原始统计量之间的差异
  beta2 <- sum(err * err)            # 计算差异平方和
  beta3 <- sum(err * err * err)      # 计算差异立方和
  a <- beta3 / beta2^(1.5)          # 计算Aest统计量
  a / 6                              # 返回Aest统计量的值
}


#' 计算BCA置信区间
#'
#' 这个函数用于计算基于bootstrap方法的BCA（Bootstrap Percentile Confidence Interval）置信区间。
#' BCA置信区间是一种基于bootstrap方法的置信区间，它考虑了数据的整体分布。
#'
#' @param a 参数a
#' @param boot 一个向量
#' @param alpha 一个概率向量，通常包含一个或多个置信水平。
#' @param theta 参数theta，用于计算booted样本中不大于theta的比例
#'
#' @return 一个矩阵，包含BCA置信区间的下界和上界。
#'
#' @export
#'
#' @examples
#' a <- 0.5
#' boot <- c(7, 11, 15, 16, 20, 22, 24, 25, 29, 33)
#' alpha <- 0.05
#' theta <- 20
#' bca_interval_result <- bca.interval(a, boot, alpha, theta)
#' print(bca_interval_result)
bca.interval <- function(a, boot, alpha, theta) {
  zp <- qnorm(1 - alpha / 2)  # 计算正态分布的分位数，用于计算置信区间
  tmp <- length(boot[boot <= theta]) / length(boot)  # 计算booted样本中不大于theta的比例
  c <- qnorm(tmp)                                   # 计算对应于比例tmp的z分数
  L <- c + (c - zp) / (1 - a * (c - zp))           # 计算BCA置信区间的下界
  U <- c + (c + zp) / (1 - a * (c + zp))           # 计算BCA置信区间的上界
  lower <- quantile(boot, pnorm(L))                 # 计算booted样本的置信区间的下界
  upper <- quantile(boot, pnorm(U))                 # 计算booted样本的置信区间的上界
  result <- rbind(lower, upper)                     # 创建一个包含上下界的矩阵
  colnames(result) <- as.character(1 - alpha)       # 设置列名
  rownames(result) <- c('BCA lower limit: ', 'BCA upper limit: ') # 设置行名
  result                                            # 返回结果矩阵
}



#' 多元线性回归中的residuals bootstrap
#' @description 多元线性回归中的residuals bootstrap是一种半参数方法，所依赖的假设为误差独立且同方差。
#' @param data 数据类型为dataframe的数据框
#' @param formula 类似于y ~ x1 + x2的表达式
#' @param rep 重抽样次数
#' @param conf_level 置信水平
#'
#' @return 返回系数估计值和对应的置信区间
#' @export
#'
#' @examples \donttest{
#' resboot(data = bodyfat,
#'         formula = bodyfat ~ Weight + Chest + Abdomen,
#'         rep = 1000,
#'         conf_level = 0.95)
#' }
#'
resboot <- function(data, formula, rep, conf_level){
  model <- lm(formula, data = data) #明确回归模型
  boot_population <- resid(model)   #确定bootstrap总体
  coef_estimates <- coef(model)     #得到参数估计
  fitted_values <- fitted(model)    #得到y_hat
  boot_coefs <- matrix(NA, rep, length(coef_estimates))
  terms <- all.vars(formula)[-1]   # 提取特征变量名称

  for (i in 1:rep) {
    boot_residuals <- sample(boot_population, nrow(data), replace = TRUE) # 从残差中有放回抽样得到bootstrap sample
    data$Y_boot <- fitted_values + boot_residuals          # compute bootstrapped y values
    boot_formula <- as.formula(paste("Y_boot ~",           # 更新公式，使Y_boot作为响应变量
                                     paste(terms, collapse = " + ")))
    fit_boot <- lm(boot_formula, data = data)              # 拟合新的模型，特征变量保持不变
    boot_coefs[i, ] <- coef(fit_boot)                      # 记录bootstrap样本的系数估计值
  }                                                        # 得到rep组bootstrapped coefficients

  #计算置信区间
  alpha <- (1 - conf_level) / 2
  lower_bounds <- apply(boot_coefs, 2, quantile, probs = alpha)
  upper_bounds <- apply(boot_coefs, 2, quantile, probs = 1 - alpha)
  conf_intervals <- cbind(lower_bounds, upper_bounds)

  #返回结果（返回系数估计值和置信区间）
  result <- list(
    estimates = coef_estimates,
    bootstrapped_cf = conf_intervals
  )
  return(result)
}


#' 使用Bootstrap方法构造两样本均值之差的置信区间
#'
#' @param x 数值向量，第一个样本
#' @param y 数值向量，第二个样本
#' @param conf_level 数值，置信水平（默认0.95）
#' @param n_bootstrap 整数，自助法的抽样次数（默认1000）
#'
#' @return 一个列表，包括以下内容：
#' \item{observed_diff}{原始样本的均值之差}
#' \item{conf_int}{置信区间的上下界}
#' \item{conf_level}{置信水平}
#' @export
#'
#' @examples
#' set.seed(123)
#' x <- rnorm(30, mean = 5, sd = 2)
#' y <- rnorm(30, mean = 3, sd = 3)
#' result <- bootstrap_mean_diff_ci(x, y, conf_level = 0.95, n_bootstrap = 1000)
#' print(result)

bootstrap_mean_diff_ci <- function(x, y, conf_level = 0.95, n_bootstrap = 1000) {
  # 检查输入
  if (!is.numeric(x) || !is.numeric(y)) {
    stop("x 和 y 必须是数值向量")
  }

  # 样本大小
  n_x <- length(x)
  n_y <- length(y)

  # 样本均值之差
  observed_diff <- mean(x) - mean(y)

  # 初始化一个向量来存储每次bootstrap的结果
  bootstrap_diffs <- numeric(n_bootstrap)

  # 进行bootstrap抽样
  for (i in 1:n_bootstrap) {
    sample_x <- sample(x, n_x, replace = TRUE)
    sample_y <- sample(y, n_y, replace = TRUE)
    bootstrap_diffs[i] <- mean(sample_x) - mean(sample_y)
  }

  # 计算置信区间
  alpha <- 1 - conf_level
  lower_bound <- quantile(bootstrap_diffs, alpha / 2)
  upper_bound <- quantile(bootstrap_diffs, 1 - alpha / 2)

  # 返回结果
  list(
    observed_diff = observed_diff,
    conf_int = c(lower_bound, upper_bound),
    conf_level = conf_level
  )
}

#' 多因子试验设计中的bootstrap方法
#'
#' @param data 数据类型为dataframe的数据框
#' @param response 目标变量，如"lifetime"
#' @param factor 感兴趣的因子，如交互效应"temp:material"
#' @param formula 用于anova的一个式子
#' @param rep 抽样次数
#'
#' @return F统计量和bootstrapped p-value
#' @export
#'
#' @examples \donttest{
#' FactDesignBoot(data = battery,
#' response = "lifetime",
#' factor = "temp:material",
#' formula = lifetime ~ temp * material,
#' rep = 1000)
#' }
FactDesignBoot <- function(data, response, factor, formula, rep) {
  # 计算原始数据的F统计量
  aov <- anova(lm(formula, data = data))
  obs_stat <- aov[factor, "F value"]

  # 初始化bootstrap统计量向量
  boot_F <- numeric(rep)
  #进行bootstrap
  for (i in 1:rep) {
    # 对响应变量进行重采样
    new_Y <- sample(data[[response]], length(data[[response]]), replace = TRUE)
    # 计算新的F统计量
    boot_anova <- anova(lm(as.formula(paste("new_Y", as.character(formula)[3],
                                            sep = "~")), data = data))
    boot_F[i] <- boot_anova[factor, "F value"]
  }
  # 计算p值
  p_value <- mean(boot_F >= obs_stat)
  return(list(F_value = obs_stat, p_value = p_value))
}

#' 使用Aligned-Rank Test检验裂区试验中的交互效应
#' @description 通过剥离wholeplot因子的效应和withinblock因子的效应并转化为相应的秩，再使用anova
#'    可以探究wholeplot因子和withinblock因子之间的交互效应
#' @param data 数据类型为dataframe的数据框
#' @param wholeplot wholeplot因子
#' @param withinblock withinblock因子
#' @param random 作为随机效应的因子
#' @param response 响应变量
#'
#' @return 方差分析表格
#' @export
#' @import tidyverse
#' @examples \donttest{
#' AlignedRankSP(data = dog,
#'               wholeplot = "diabetes",
#'               withinblock = "method",
#'               random = "dogs",
#'               response = "turnover")
#' }
AlignedRankSP <- function(data, wholeplot, withinblock, random, response){
  library(tidyverse)
  baryij <- data %>%
    group_by(across(all_of(c(wholeplot, random)))) %>%
    summarise(mean_ij = mean(get(response)))

  baryk <- data %>% group_by(across(all_of(withinblock))) %>%
    summarise(mean_k = mean(get(response)))

  bary <- mean(data[[response]])

  result <- data %>%
    left_join(baryij, by = c(wholeplot, random)) %>%
    left_join(baryk, by = withinblock) %>%
    mutate(aligned_value = get(response) - mean_ij - mean_k + bary)
  aov_formula <- as.formula(paste("aligned_value ~", wholeplot, "+", random, "+",
                                  withinblock, "+",
                                  paste(wholeplot, withinblock, sep = ":")))
  aov_result <- summary(aov(aov_formula, data = result))
  return(aov_result)
}

# 置换检验函数
#' Title: 因子整体效应的置换检验
#'
#' @param data 数据框，包含因子A、因子B和响应变量
#' @param factor_A 因子A的列名（字符串）
#' @param factor_B 因子B的列名（字符串）
#' @param response 响应变量的列名（字符串）
#' @param num_permutations 置换次数，默认值为5000
#'
#' @return 返回一个列表，包括观测到的F统计量、p值和置换分布的F统计量
#' @export
#'
#' @examples
#' # 示例数据
#' data <- data.frame(
#'   A = rep(c("A1", "A2", "A3"), each = 9),
#'   B = rep(rep(c("B1", "B2", "B3"), each = 3), times = 3),
#'   response = c(1.38, 2.54, 2.14, 2.32, 3.77, 2.59, 4.57, 3.69, 4.41,
#'                1.92, 2.48, 2.10, 3.64, 2.72, 2.58, 3.78, 4.17, 4.31,
#'                1.84, 2.34, 2.12, 3.56, 2.70, 2.60, 3.72, 4.08, 4.29)
#' )
#'
#' # 执行置换检验
#' result <- permutation_test(data, "A", "B", "response", num_permutations = 5000)
#'
#' # 输出结果
#' print(result)

permutation_test <- function(data, factor_A, factor_B, response, num_permutations = 5000) {
  # data: 数据框，包含因子A、因子B和响应变量
  # factor_A: 因子A的列名（字符串）
  # factor_B: 因子B的列名（字符串）
  # response: 响应变量的列名（字符串）
  # num_permutations: 置换次数

  library(car)

  # 完整模型的平方误差和与自由度
  full_model <- lm(as.formula(paste(response, "~", factor_A, "*", factor_B)), data = data)
  #使用lm函数拟合包含因子A、因子B及其交互作用的完整模型。

  SSE_full <- sum(residuals(full_model)^2)
  #计算完整模型的平方误差和（SSE_full）。

  df_full <- df.residual(full_model)
  #获取完整模型的自由度（df_full）。


  # 简化模型的平方误差和与自由度
  reduced_model <- lm(as.formula(paste(response, "~", factor_B)), data = data)
  #使用lm函数拟合仅包含因子B的简化模型。
  SSE_reduced <- sum(residuals(reduced_model)^2)
  #简化模型的平方误差和
  df_reduced <- df.residual(reduced_model)
  #简化模型的自由度

  # 观测到的F统计量
  F_obs <- ((SSE_reduced - SSE_full) / (df_reduced - df_full)) / (SSE_full / df_full)

  # 置换分布
  F_perm <- numeric(num_permutations)

  set.seed(123)  # 设置种子以确保结果可重复
  for (i in 1:num_permutations) {
    # 对因子B的每个水平进行置换
    perm_data <- data
    for (level in unique(data[[factor_B]])) {
      subset_indices <- which(data[[factor_B]] == level)
      perm_data[subset_indices, factor_A] <- sample(data[subset_indices, factor_A])
    }

    # 计算置换样本的F统计量
    perm_full_model <- lm(as.formula(paste(response, "~", factor_A, "*", factor_B)), data = perm_data)
    SSE_full_perm <- sum(residuals(perm_full_model)^2)
    F_perm[i] <- ((SSE_reduced - SSE_full_perm) / (df_reduced - df_full)) / (SSE_full_perm / df_full)
  }

  # 计算p值
  p_value <- mean(F_perm >= F_obs)

  return(list(F_obs = F_obs, p_value = p_value, F_perm = F_perm))
}



#' 相关系数Bootstrap区间估计
#'
#' @param x 向量，数据集的第一个变量
#' @param y 向量，数据集的第二个变量
#' @param R 整数，Bootstrap抽样次数
#' @param method 字符串，区间估计的方法，可选"bca"或"perc"
#'
#' @return 一个列表，包含原始相关系数和Bootstrap置信区间
#' @import boot
#' @export
#'
#' @examples \donttest{
#' x <- rnorm(100)
#' y <- 0.5 * x + rnorm(100)
#' bootstrap_cor_ci(x, y, R = 1000, method = "bca")
#' }
bootstrap_cor_ci <- function(x, y, R = 1000, method = "bca") {
  library(boot)
  # 检查输入数据
  if (length(x) != length(y)) {
    stop("x和y的长度必须相同")
  }
  if (!method %in% c("bca", "perc")) {
    stop("method参数必须为'bca'或'perc'")
  }

  # 将数据组合成一个数据框
  data <- data.frame(x, y)

  # 自定义函数，计算相关系数
  correlation <- function(data, indices) {
    d <- data[indices, ]  # 根据indices取样
    return(cor(d$x, d$y)) # 计算相关系数
  }

  # 计算原始相关系数
  original_cor <- cor(x, y)

  # 进行Bootstrap抽样并计算相关系数
  bootstrap_results <- boot(data, statistic = correlation, R = R)

  # 计算置信区间
  bootstrap_ci <- boot.ci(bootstrap_results, type = method)

  # 返回结果
  return(list(
    original_cor = original_cor,
    bootstrap_ci = bootstrap_ci
  ))
}


#' 固定X的Bootstrap抽样
#'
#' @param x 向量，数据集的第一个变量（自变量）
#' @param y 向量，数据集的第二个变量（因变量）
#' @param R 整数，Bootstrap抽样次数
#' @param method 字符串，区间估计的方法，可选"bca"或"perc"
#' @param return_samples 逻辑值，是否返回Bootstrap样本，默认为FALSE
#'
#' @return 一个列表，包含原始相关系数、Bootstrap置信区间和（可选的）Bootstrap样本
#' @import boot
#' @export
#'
#' @examples \donttest{
#' x <- rnorm(100)
#' y <- 0.5 * x + rnorm(100)
#' result <- fixed_x_bootstrap_cor_ci(x, y, R = 1000, method = "bca", return_samples = FALSE)
#' str(result)
#' }
fixed_x_bootstrap_cor_ci <- function(x, y, R = 1000, method = "bca", return_samples = FALSE) {
  library(boot)
  # 检查输入数据
  if (length(x) != length(y)) {
    stop("x和y的长度必须相同")
  }
  if (!method %in% c("bca", "perc")) {
    stop("method参数必须为'bca'或'perc'")
  }

  # 拟合回归模型，估计h(X)
  model <- lm(y ~ x)
  h_x_hat <- fitted(model)

  # 计算残差
  residuals <- residuals(model)

  # 初始化存储Bootstrap样本的列表
  bootstrap_samples <- if (return_samples) vector("list", R) else NULL

  # 自定义函数，计算相关系数并存储Bootstrap样本
  correlation <- function(data, indices) {
    x_boot <- data$x
    e_boot <- sample(data$residuals, length(data$residuals), replace = TRUE)
    y_boot <- data$h_x_hat + e_boot

    # 存储当前的Bootstrap样本
    if (return_samples) {
      bootstrap_samples[[indices[1]]] <<- data.frame(x = x_boot, y = y_boot)
    }

    return(cor(x_boot, y_boot))
  }

  # 准备数据框
  data <- data.frame(x = x, h_x_hat = h_x_hat, residuals = residuals)

  # 计算原始相关系数
  original_cor <- cor(x, y)

  # 进行Bootstrap抽样并计算相关系数
  bootstrap_results <- boot(data, statistic = correlation, R = R)

  # 计算置信区间
  bootstrap_ci <- boot.ci(bootstrap_results, type = method)

  # 返回结果
  result <- list(
    original_cor = original_cor,
    bootstrap_ci = bootstrap_ci
  )

  if (return_samples) {
    result$bootstrap_samples <- bootstrap_samples
  }

  return(result)
}


#' 回归斜率的Bootstrap区间估计
#'
#' @param x 向量，数据集的第一个变量（自变量）
#' @param y 向量，数据集的第二个变量（因变量）
#' @param R 整数，Bootstrap抽样次数，默认为1000
#' @param alpha 置信水平，默认为0.05
#'
#' @return 一个列表，包含原始斜率估计、标准误差、t_枢轴量的经验分布和Bootstrap置信区间
#' @import boot
#' @export
#'
#' @examples \donttest{
#' x <- rnorm(100)
#' y <- 0.5 * x + rnorm(100)
#' result <- bootstrap_slope_ci(x, y, R = 1000, alpha = 0.05)
#' str(result)
#' }
bootstrap_slope_ci <- function(x, y, R = 1000, alpha = 0.05) {
  library(boot)
  # 检查输入数据
  if (length(x) != length(y)) {
    stop("x和y的长度必须相同")
  }

  # 计算原始样本的最小二乘估计
  model <- lm(y ~ x)
  beta_0_hat <- coef(model)[1]
  beta_1_hat <- coef(model)[2]
  residuals <- residuals(model)
  n <- length(x)

  # 自定义函数，计算新的斜率估计和t_枢轴量
  calc_slope_t <- function(data, indices) {
    x_boot <- data$x
    e_boot <- sample(data$residuals, n, replace = TRUE)
    y_boot <- data$h_x_hat + e_boot

    # 拟合新的回归模型
    model_boot <- lm(y_boot ~ x_boot)
    beta_1_hat_boot <- coef(model_boot)[2]
    se_beta_1_hat_boot <- summary(model_boot)$coefficients[2, 2]

    # 计算t_枢轴量
    t_epsilon <- beta_1_hat_boot / se_beta_1_hat_boot
    return(t_epsilon)
  }

  # 准备数据框
  data <- data.frame(x = x, h_x_hat = beta_0_hat + beta_1_hat * x, residuals = residuals)

  # 进行Bootstrap抽样并计算t_枢轴量
  t_epsilon_values <- boot(data, statistic = calc_slope_t, R = R)$t

  # 计算t_枢轴量的分位数
  t_epsilon_alpha <- quantile(t_epsilon_values, probs = c(alpha / 2, 1 - alpha / 2))

  # 计算标准误差
  se_beta_1_hat <- summary(model)$coefficients[2, 2]

  # 计算Bootstrap置信区间
  ci_lower <- beta_1_hat - se_beta_1_hat * t_epsilon_alpha[2]
  ci_upper <- beta_1_hat + se_beta_1_hat * t_epsilon_alpha[1]

  # 返回结果
  result <- list(
    beta_1_hat = beta_1_hat,
    se_beta_1_hat = se_beta_1_hat,
    ci_lower = ci_lower,
    ci_upper = ci_upper
  )

  return(result)
}


#' Aligning Data in Completely Randomized Design
#' This function aligns the data for a completely randomized design by
#' adjusting each data point based on its group median and the overall median.
#'
#' @param data A data frame containing the data
#' @param group_col The name of the column representing the group variable
#' @param value_col The name of the column representing the value variable
#'
#' @return A data frame with an additional column aligned_value containing the aligned data.
#' @export
#'
#' @examples \donttest{
#' data <- data.frame(
#'   group = rep(1:3, each = 5),
#'   value = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)
#' )
#' align_data(data, "group", "value")
#' }
align_data <- function(data, group_col, value_col) {
  if (!is.data.frame(data)) {
    stop("Input data must be a data frame")
  }

  if (!(group_col %in% names(data))) {
    stop(paste("Column", group_col, "does not exist in the data"))
  }

  if (!(value_col %in% names(data))) {
    stop(paste("Column", value_col, "does not exist in the data"))
  }

  overall_median <- median(data[[value_col]])

  medians <- aggregate(data[[value_col]], by = list(data[[group_col]]), FUN = median)
  names(medians) <- c(group_col, "group_median")

  aligned_data <- merge(data, medians, by = group_col)
  aligned_data$aligned_value <- aligned_data[[value_col]] - aligned_data$group_median + overall_median

  return(aligned_data)
}
