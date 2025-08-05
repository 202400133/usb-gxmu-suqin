# 加载必要包
library(readr)
library(dplyr)
library(broom)

# 读取数据
data <- read_csv("output/2025.8.4/filtered_postbaseline_CVD_with_isCVD.csv")

# 结局变量：CVD列表
cvd_vars <- c(
  "is_participant.p131382", "is_participant.p131348", "is_participant.p131350", "is_participant.p131352",
  "is_participant.p131342", "is_participant.p131344", "is_participant.p131354", "is_participant.p131296",
  "is_participant.p131298", "is_participant.p131304", "is_participant.p131306", "is_participant.p131314",
  "is_participant.p131316", "is_participant.p131386", "is_participant.p131320", "is_participant.p131366",
  "is_participant.p131362", "is_participant.p131056", "is_participant.p131322", "is_participant.p131324",
  "is_participant.p131328", "is_participant.p131400", "is_participant.p131308"
)

# 协变量：Model 2 中使用的协变量
covariates <- c(
  "participant.p26223", "participant.p22009_a1", "participant.p22009_a2", "participant.p22009_a3",
  "participant.p22009_a4", "participant.p22009_a5", "participant.p22009_a6", "participant.p22009_a7",
  "participant.p22009_a8", "participant.p22009_a9", "participant.p22009_a10", "participant.p22000"
)

# 初始化结果表
results <- data.frame(
  CVD = character(),
  Adjusted_OR_95CI = character(),
  P_value = character(),
  stringsAsFactors = FALSE
)

# 遍历每一个CVD变量进行回归分析
for (cvd in cvd_vars) {
  formula_str <- paste0(cvd, " ~ is_RPL + ", paste(covariates, collapse = " + "))
  formula <- as.formula(formula_str)
  
  model_data <- data %>%
    select(all_of(c(cvd, "is_RPL", covariates))) %>%
    na.omit()
  
  if (nrow(model_data) > 50) {
    fit <- try(glm(formula, data = model_data, family = binomial), silent = TRUE)
    
    if (!inherits(fit, "try-error")) {
      fit_summary <- tidy(fit, conf.int = TRUE, exponentiate = TRUE)
      rpl_row <- fit_summary %>% filter(term == "is_RPL")
      
      if (nrow(rpl_row) == 1) {
        or_ci <- paste0(
          round(rpl_row$estimate, 2), "(", 
          round(rpl_row$conf.low, 2), "-", 
          round(rpl_row$conf.high, 2), ")"
        )
        p_val <- formatC(rpl_row$p.value, format = "f", digits = 3)
      } else {
        or_ci <- NA
        p_val <- NA
      }
    } else {
      or_ci <- NA
      p_val <- NA
    }
  } else {
    or_ci <- NA
    p_val <- NA
  }
  
  results <- rbind(results, data.frame(
    CVD = cvd,
    Adjusted_OR_95CI = or_ci,
    P_value = p_val,
    stringsAsFactors = FALSE
  ))
}

# 创建输出目录并保存
output_dir <- "output/2025.8.5"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
write_csv(results, file.path(output_dir, "CVD_RPL_AdjustedModel2_results.csv"))
