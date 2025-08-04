# 第一步提取所有29种cvd疾病发病日期，提取最早发病日期。
library(readr)
library(dplyr)

# 1. 读取文件
df <- read_csv("output/filtered_merged_data_p3839_filtered.csv")

# 2. 指定 CVD 发病日期相关列（共 29 个字段 ID）
cvd_cols <- c(
  "participant.p131382", "participant.p131350", "participant.p131352", "participant.p131342",
  "participant.p131344", "participant.p131354", "participant.p131296", "participant.p131298",
  "participant.p131300", "participant.p131302", "participant.p131304", "participant.p131306",
  "participant.p131334", "participant.p131314", "participant.p131316", "participant.p131318",
  "participant.p131386", "participant.p131320", "participant.p131366", "participant.p131362",
  "participant.p131056", "participant.p131322", "participant.p131324", "participant.p131326",
  "participant.p131328", "participant.p131330", "participant.p131400", "participant.p131308",
  "participant.p131348"
)

# 3. 处理日期列：确保它们是 Date 类型（readr 已自动识别）
# 4. 计算每行的最早发病日期和是否有病
df <- df %>%
  mutate(
    CVD_date = do.call(pmin, c(across(all_of(cvd_cols)), na.rm = TRUE)),  # 最早日期
    is_CVD = if_else(!is.na(CVD_date), 1, 0)  # 是否患过任何一种 CVD
  )

# 5. 保存结果
write_csv(df, "output/filtered_merged_data_with_CVD_info.csv")

# 可视化查看
View(df)

# 加载 tidyverse（如果尚未加载）
library(dplyr)
library(readr)


#第二步与回答RPL时间即p53字段比较，排除入组前患CVD病例

library(readr)
library(dplyr)

# 1. 读取数据
df <- read_csv("output/filtered_merged_data_with_CVD_info.csv")

# 2. 确保日期列为 Date 类型
df$CVD_date <- as.Date(df$CVD_date)
df$participant.p53_i0 <- as.Date(df$participant.p53_i0)

# 3. 原始行数
original_n <- nrow(df)

# 4. 保留：CVD_date > 入组时间，或 CVD_date 为 NA（未患病）
df_filtered <- df %>%
  filter(is.na(CVD_date) | CVD_date > participant.p53_i0)

# 5. 统计
filtered_n <- nrow(df_filtered)
removed_n <- original_n - filtered_n

# 6. 输出统计信息
cat("📊 原始样本数：", original_n, "\n")
cat("✅ 保留的样本数（CVD_date > 入组 或 未患病）：", filtered_n, "\n")
cat("❌ 剔除的样本数（入组前已发病）：", removed_n, "\n")

# 7. 保存新文件
write_csv(df_filtered, "output/filtered_postbaseline_CVD.csv")

#第三步流产次数在字段 3839>=3为病例组，<3为非病例组
library(readr)
library(dplyr)

# 1. 读取原文件
df <- read_csv("output/filtered_postbaseline_CVD.csv")

# 2. 添加 is_RPL 列：>=3 为 1，<3 或 NA 为 0
df <- df %>%
  mutate(
    is_RPL = if_else(!is.na(participant.p3839_i0) & participant.p3839_i0 >= 3, 1, 0)
  )

# 3. 统计结果（可选）
table_RPL <- table(df$is_RPL)
cat("📊 RPL 分组统计：\n")
cat("非RPL组（is_RPL = 0）：", table_RPL["0"], "人\n")
cat("RPL组（is_RPL = 1）：", table_RPL["1"], "人\n")

# 4. 保存回原文件（覆盖写入）
write_csv(df, "output/filtered_postbaseline_CVD.csv")

#第五步横断面研究计算RPL组与非RPL组总的CVD患病率比较

# 如未安装，先安装 epitools 包（只需运行一次）
if (!require(epitools)) {
  install.packages("epitools")
  library(epitools)
} else {
  library(epitools)
}

# 加载数据
data <- read.csv("output/filtered_postbaseline_CVD.csv")

# 创建二维列联表
table_RPL_CVD <- table(data$is_RPL, data$is_CVD)
colnames(table_RPL_CVD) <- c("未患CVD", "患CVD")
rownames(table_RPL_CVD) <- c("非RPL组", "RPL组")
print(table_RPL_CVD)

# 卡方检验
chisq_result <- chisq.test(table_RPL_CVD)
cat("✅ 卡方检验 p 值为：", chisq_result$p.value, "\n")

# 计算比值比（OR）和95%置信区间
or_result <- oddsratio(table_RPL_CVD)
cat("✅ OR 值：\n")
print(or_result$measure)

cat("✅ 置信区间：\n")
print(or_result$conf.int)

cat("✅ p 值（Fisher 精确检验）：\n")
print(or_result$p.value)


# 计算比值比（OR）和95%置信区间
or_result <- oddsratio(table_RPL_CVD)
print(or_result$measure)         # OR值
print(or_result$p.value)         # p值
print(or_result$conf.int)        # 置信区间

#第六步将每个cvd疾病变成二分类
# 读取数据
data <- read.csv("output/filtered_postbaseline_CVD.csv")

# CVD 日期列向量
cvd_date_cols <- c(
  "participant.p131382", "participant.p131350", "participant.p131352", "participant.p131342",
  "participant.p131344", "participant.p131354", "participant.p131296", "participant.p131298",
  "participant.p131300", "participant.p131302", "participant.p131304", "participant.p131306",
  "participant.p131334", "participant.p131314", "participant.p131316", "participant.p131318",
  "participant.p131386", "participant.p131320", "participant.p131366", "participant.p131362",
  "participant.p131056", "participant.p131322", "participant.p131324", "participant.p131326",
  "participant.p131328", "participant.p131330", "participant.p131400", "participant.p131308",
  "participant.p131348"
)

# 初始化统计结果列表
cvd_counts <- list()

# 遍历每列生成新列 is_xxx（是否患病），并统计患病人数
for (col in cvd_date_cols) {
  new_col <- paste0("is_", col)
  data[[new_col]] <- ifelse(is.na(data[[col]]), 0, 1)
  cvd_counts[[col]] <- sum(data[[new_col]])
}

# 输出每种 CVD 的患病人数
cat("📊 各种 CVD 疾病患病人数统计：\n")
for (col in names(cvd_counts)) {
  cat(col, ": ", cvd_counts[[col]], "\n")
}

# 保存更新后的数据表
write.csv(data, "output/filtered_postbaseline_CVD_with_isCVD.csv", row.names = FALSE)

#第七步除去发病率为0的cvd疾病，分别统计每一个cvd疾病在两组病人中的发病率比较
# 读取更新后的数据
data <- read.csv("output/filtered_postbaseline_CVD_with_isCVD.csv")

# 所有 CVD 标志列
cvd_flags <- grep("^is_participant\\.p131", names(data), value = TRUE)

# 剔除病例数为0的 CVD 类型
cvd_nonzero <- cvd_flags[sapply(data[cvd_flags], function(x) sum(x, na.rm = TRUE) > 0)]

# 初始化结果存储
result_list <- data.frame(
  CVD_Type = character(),
  RPL_Cases = integer(),
  NonRPL_Cases = integer(),
  OR = numeric(),
  `95% CI Lower` = numeric(),
  `95% CI Upper` = numeric(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)

# 遍历每个 CVD 做横断面分析
for (cvd in cvd_nonzero) {
  # 构造二维列联表
  tab <- table(data[[cvd]], data$is_RPL)
  if (all(dim(tab) == c(2, 2))) {
    test <- chisq.test(tab, correct = FALSE)
    or <- (tab[2,2] * tab[1,1]) / (tab[1,2] * tab[2,1])
    se_log_or <- sqrt(1/tab[2,2] + 1/tab[1,1] + 1/tab[1,2] + 1/tab[2,1])
    ci_lower <- exp(log(or) - 1.96 * se_log_or)
    ci_upper <- exp(log(or) + 1.96 * se_log_or)
    
    result_list <- rbind(result_list, data.frame(
      CVD_Type = cvd,
      RPL_Cases = tab[2,2],
      NonRPL_Cases = tab[2,1],
      OR = round(or, 2),
      `95% CI Lower` = round(ci_lower, 2),
      `95% CI Upper` = round(ci_upper, 2),
      P_value = signif(test$p.value, 3)
    ))
  }
}

# 打印三线表格式（示意，实际应用 xtable 或 kable 输出为论文格式）
print(result_list)

# 可选择保存为 CSV
write.csv(result_list, "output/CVD_cross_sectional_RPL_vs_nonRPL.csv", row.names = FALSE)

# 读取数据
data <- read.csv("output/CVD_cross_sectional_RPL_vs_nonRPL.csv")

# 确保 P_value 列为数值型
data$P_value <- as.numeric(data$P_value)

# 统计 P 值小于 0.05 的行数
sig_count <- sum(data$P_value < 0.05, na.rm = TRUE)

# 输出结果
cat("P值小于0.05的CVD种类数量为：", sig_count, "\n")
