#死亡原因统计
#第一步先根据字段40001找出全因死亡和CVD原因死亡
library(readr)
library(dplyr)
library(stringr)

# 读取数据
file_path <- "output/2025.8.4/filtered_postbaseline_CVD_with_isCVD.csv"
df <- read_csv(file_path)

# 生成 CVD_death 列（ICD-10代码以"I"开头代表CVD死亡）
df <- df %>%
  mutate(
    CVD_death = if_else(
      !is.na(participant.p40001_i0.y) & str_starts(participant.p40001_i0.y, "I"),
      1, 0
    ),
    # 生成 All_cause_death 列（只要非空就代表死亡）
    All_cause_death = if_else(
      is.na(participant.p40001_i0.y),
      0, 1
    )
  )

# 写入原路径（覆盖保存）
write_csv(df, file_path)
# 加载必要库
library(dplyr)
library(readr)

# 读取数据
file_path <- "output/2025.8.4/filtered_postbaseline_CVD_with_isCVD.csv"
df <- read_csv(file_path)

# 统计0/1频数
death_counts <- df %>%
  summarise(
    CVD_death_0 = sum(CVD_death == 0, na.rm = TRUE),
    CVD_death_1 = sum(CVD_death == 1, na.rm = TRUE),
    All_cause_death_0 = sum(All_cause_death == 0, na.rm = TRUE),
    All_cause_death_1 = sum(All_cause_death == 1, na.rm = TRUE)
  )

# 打印结果
print(death_counts)

#第二步统计两组病人在全因死亡和cvd死亡中的例数

library(dplyr)
library(readr)

# 读取数据
file_path <- "output/2025.8.4/filtered_postbaseline_CVD_with_isCVD.csv"
df <- read_csv(file_path)

# 按 is_RPL 和 CVD_death 分组统计
cvd_death_counts <- df %>%
  group_by(is_RPL, CVD_death) %>%
  summarise(n = n(), .groups = "drop")

# 按 is_RPL 和 All_cause_death 分组统计
allcause_death_counts <- df %>%
  group_by(is_RPL, All_cause_death) %>%
  summarise(n = n(), .groups = "drop")

# 输出查看
print(cvd_death_counts)
print(allcause_death_counts)

#第三步用卡方检验比较两组病人在两个死亡原因的分布上是否有统计学意义
install.packages("epiR")  # 只需运行一次
library(epiR)

library(readr)
library(dplyr)
library(epitools)  # 用于计算OR和置信区间

# 读取数据
file_path <- "output/2025.8.4/filtered_postbaseline_CVD_with_isCVD.csv"
df <- read_csv(file_path)

# 函数：卡方检验 + OR 计算
calculate_or_chisq <- function(data, outcome_col, group_col = "is_RPL") {
  # 创建2x2列联表
  tbl <- table(data[[group_col]], data[[outcome_col]])
  
  # 卡方检验
  chisq <- chisq.test(tbl, correct = FALSE)
  
  # OR计算（使用fisher.test可以获得更精确CI）
  or_result <- oddsratio(tbl)
  
  result <- data.frame(
    Outcome = outcome_col,
    OR = round(or_result$measure[2,1], 3),
    Lower_CI = round(or_result$measure[2,2], 3),
    Upper_CI = round(or_result$measure[2,3], 3),
    P_value = chisq$p.value
  )
  
  return(result)
}

# 对两个死亡结局分别计算
res_cvd <- calculate_or_chisq(df, "CVD_death")
res_allcause <- calculate_or_chisq(df, "All_cause_death")

# 合并结果
results <- bind_rows(res_cvd, res_allcause)

# 打印
print(results)


