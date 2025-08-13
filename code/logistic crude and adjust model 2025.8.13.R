# logistics回归分析粗模型、调整1、2、3、4模型
library(dplyr)
library(broom)
library(readr)
library(ggplot2)
library(forcats)

# 读取数据
data <- read_csv("output/2025.8.11/participant_with_htn_med_flags.csv")

# 定义结局变量（22种CVD）
cvd_vars <- c(
  "is_participant.p131382", "is_participant.p131348", "is_participant.p131350", "is_participant.p131352",
  "is_participant.p131342", "is_participant.p131344", "is_participant.p131354", "is_participant.p131296",
  "is_participant.p131298", "is_participant.p131304", "is_participant.p131306", "is_participant.p131314",
  "is_participant.p131316", "is_participant.p131386", "is_participant.p131366",
  "is_participant.p131362", "is_participant.p131056", "is_participant.p131322", "is_participant.p131324",
  "is_participant.p131328", "is_participant.p131400", "is_participant.p131308"
)

# 模型协变量设置
model1_cov <- c(
  "participant.p21022",
  "participant.p21001_i0",
  "participant.p20116_i0",
  "participant.p1558_i0...3",
  "participant.p1488_i0",
  "participant.p1160_i0",
  "participant.p1498_i0",
  "participant.p884_i0"
)
model2_cov <- c(model1_cov, "participant.p4080_i0_a1","participant.p4079_i0_a1","participant.p30690_i0","participant.p30760_i0","participant.p30710_i0")
model3_cov <- c(model2_cov, "hypertension_status")
model4_cov <- c(model3_cov, "participant.p3581_i0","participant.p2814_i0","participant.p26223",
                "participant.p22009_a1","participant.p22009_a2","participant.p22009_a3","participant.p22009_a4",
                "participant.p22009_a5","participant.p22009_a6","participant.p22009_a7","participant.p22009_a8",
                "participant.p22009_a9","participant.p22009_a10")

# 用于保存所有模型结果
all_results <- data.frame()

# 循环每个 CVD 变量
for (cvd in cvd_vars) {
  for (model_type in c("Crude", "Model 1", "Model 2", "Model 3", "Model 4")) {
    covariates <- switch(model_type,
                         "Crude"   = c(),
                         "Model 1" = model1_cov,
                         "Model 2" = model2_cov,
                         "Model 3" = model3_cov,
                         "Model 4" = model4_cov)
    
    formula_str <- paste(cvd, "~ is_RPL", if (length(covariates) > 0) paste("+", paste(covariates, collapse = " + ")) else "")
    formula <- as.formula(formula_str)
    
    model_data <- data %>%
      select(any_of(c(cvd, "is_RPL", covariates))) %>%  # <- 用 any_of 更稳
      na.omit()
    
    if (nrow(model_data) > 50) {
      fit <- try(glm(formula, data = model_data, family = binomial), silent = TRUE)
      
      if (!inherits(fit, "try-error")) {
        tidy_fit <- tidy(fit, conf.int = TRUE, exponentiate = TRUE)
        rpl_row <- tidy_fit %>% filter(term == "is_RPL")
        
        if (nrow(rpl_row) == 1) {
          OR    <- round(rpl_row$estimate, 2)
          lower <- round(rpl_row$conf.low, 2)
          upper <- round(rpl_row$conf.high, 2)
          pval  <- formatC(rpl_row$p.value, format = "e", digits = 2)
          or_ci_str <- paste0(OR, " (", lower, "-", upper, ")")
        } else {
          OR <- lower <- upper <- NA
          or_ci_str <- pval <- NA
        }
      } else {
        OR <- lower <- upper <- or_ci_str <- pval <- NA
      }
    } else {
      OR <- lower <- upper <- or_ci_str <- pval <- NA
    }
    
    all_results <- rbind(all_results, data.frame(
      CVD = cvd,
      Model = model_type,
      OR = OR,
      Lower_CI = lower,
      Upper_CI = upper,
      OR_95CI = or_ci_str,
      P_value = pval,
      stringsAsFactors = FALSE
    ))
  }
}

# 保存结果（保留你原来的文件名）
write_csv(all_results, "output/2025.8.12/CVD_RPL_Crude_Model1_results111.csv")
