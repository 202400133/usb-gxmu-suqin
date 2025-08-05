# 加载所需包
library(dplyr)
library(broom)
library(readr)

# 读取数据
data <- read_csv("output/2025.8.4/filtered_postbaseline_CVD_with_isCVD.csv")

# 定义 CVD 变量（结局）
cvd_vars <- c(
  "is_participant.p131382", "is_participant.p131348", "is_participant.p131350", "is_participant.p131352",
  "is_participant.p131342", "is_participant.p131344", "is_participant.p131354", "is_participant.p131296",
  "is_participant.p131298", "is_participant.p131304", "is_participant.p131306", "is_participant.p131314",
  "is_participant.p131316", "is_participant.p131386", "is_participant.p131320", "is_participant.p131366",
  "is_participant.p131362", "is_participant.p131056", "is_participant.p131322", "is_participant.p131324",
  "is_participant.p131328", "is_participant.p131400", "is_participant.p131308"
)

# 协变量（用于 Adjusted Model 1）
covariates <- c(
  "participant.p21022", "participant.p21001_i0", "participant.p4080_i0_a1", "participant.p2966_i0",
  "participant.p1239_i0", "participant.p30690_i0", "participant.p30760_i0", "participant.p30870_i0",
  "participant.p30630_i0", "participant.p30750_i0", "participant.p30700_i0", "participant.p30620_i0",
  "participant.p30730_i0", "participant.p30710_i0", "participant.p2814_i0", "participant.p2834_i0",
  "participant.p3581_i0", "participant.p2714_i0", "participant.p20116_i0"
)

# 创建空结果表
results <- data.frame(
  CVD = character(),
  Adjusted_OR_95CI = character(),
  P_value = character(),
  stringsAsFactors = FALSE
)

# 循环每个CVD变量进行逻辑回归
for (cvd in cvd_vars) {
  # 构建模型公式
  formula_str <- paste(cvd, "~ is_RPL +", paste(covariates, collapse = " + "))
  formula <- as.formula(formula_str)
  
  # 拟合模型，使用完整案例分析避免NA
  model_data <- data %>%
    select(all_of(c(cvd, "is_RPL", covariates))) %>%
    na.omit()
  
  if (nrow(model_data) > 50) {
    fit <- try(glm(formula, data = model_data, family = binomial), silent = TRUE)
    
    if (!inherits(fit, "try-error")) {
      # 提取OR、CI、P值
      tidy_fit <- tidy(fit, conf.int = TRUE, exponentiate = TRUE)
      rpl_row <- tidy_fit %>% filter(term == "is_RPL")
      
      if (nrow(rpl_row) == 1) {
        OR_CI <- paste0(
          round(rpl_row$estimate, 2), "(", 
          round(rpl_row$conf.low, 2), "-", 
          round(rpl_row$conf.high, 2), ")"
        )
        pval <- formatC(rpl_row$p.value, format = "f", digits = 3)
      } else {
        OR_CI <- NA
        pval <- NA
      }
    } else {
      OR_CI <- NA
      pval <- NA
    }
  } else {
    OR_CI <- NA
    pval <- NA
  }
  
  # 添加到结果表
  results <- rbind(results, data.frame(
    CVD = cvd,
    Adjusted_OR_95CI = OR_CI,
    P_value = pval,
    stringsAsFactors = FALSE
  ))
}

# 输出结果
write_csv(results, "output/2025.8.5/CVD_RPL_AdjustedModel1_results.csv")
