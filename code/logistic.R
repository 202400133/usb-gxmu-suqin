#logistic crude model
# 加载必要的包
library(dplyr)
library(broom)
library(readr)

# 设置文件路径
file_path <- "output/2025.8.4/filtered_postbaseline_CVD_with_isCVD.csv"

# 读取数据
data <- read_csv(file_path)

# CVD变量列名列表
cvd_vars <- c(
  "is_participant.p131382", "is_participant.p131348", "is_participant.p131350", "is_participant.p131352",
  "is_participant.p131342", "is_participant.p131344", "is_participant.p131354", "is_participant.p131296",
  "is_participant.p131298", "is_participant.p131304", "is_participant.p131306", "is_participant.p131314",
  "is_participant.p131316", "is_participant.p131386", "is_participant.p131320", "is_participant.p131366",
  "is_participant.p131362", "is_participant.p131056", "is_participant.p131322", "is_participant.p131324",
  "is_participant.p131328", "is_participant.p131400", "is_participant.p131308"
)

# 初始化结果表
results <- data.frame(
  CVD = character(),
  RPL_case_prop = character(),
  nonRPL_case_prop = character(),
  OR_CI = character(),
  P_value = character(),
  stringsAsFactors = FALSE
)

# 遍历每个CVD变量
for (cvd in cvd_vars) {
  # 计算病例数和比例
  rpl_total <- sum(data$is_RPL == 1, na.rm = TRUE)
  nonrpl_total <- sum(data$is_RPL == 0, na.rm = TRUE)
  
  rpl_case <- sum(data[[cvd]] == 1 & data$is_RPL == 1, na.rm = TRUE)
  nonrpl_case <- sum(data[[cvd]] == 1 & data$is_RPL == 0, na.rm = TRUE)
  
  rpl_prop <- paste0(rpl_case, " (", round(rpl_case / rpl_total, 3), ")")
  nonrpl_prop <- paste0(nonrpl_case, " (", round(nonrpl_case / nonrpl_total, 3), ")")
  
  # 逻辑回归 crude model
  formula <- as.formula(paste(cvd, "~ is_RPL"))
  model <- glm(formula, data = data, family = binomial)
  tidy_model <- tidy(model, conf.int = TRUE, exponentiate = TRUE)
  
  # 提取OR和置信区间
  if ("is_RPL" %in% tidy_model$term) {
    or_row <- tidy_model %>% filter(term == "is_RPL")
    or_ci <- paste0(round(or_row$estimate, 2), " (", 
                    round(or_row$conf.low, 2), "-", 
                    round(or_row$conf.high, 2), ")")
    p_value <- formatC(or_row$p.value, format = "f", digits = 3)
  } else {
    or_ci <- NA
    p_value <- NA
  }
  
  # 添加到结果表
  results <- rbind(results, data.frame(
    CVD = cvd,
    RPL_case_prop = rpl_prop,
    nonRPL_case_prop = nonrpl_prop,
    OR_CI = or_ci,
    P_value = p_value,
    stringsAsFactors = FALSE
  ))
}

# 保存结果为CSV文件
write_csv(results, "output/2025.8.5/CVD_RPL_association_results.csv")
