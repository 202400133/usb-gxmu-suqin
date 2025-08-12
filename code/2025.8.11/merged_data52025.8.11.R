#第一步合并新下载的Date5.csv 
library(dplyr)
library(readr)

#读取两个文件
data1 <- read_csv("data/Date5.csv")
data2 <- read_csv("output/2025.8.4/filtered_postbaseline_CVD_with_isCVD.csv")

# 合并（只保留两个表中都有的 participant.eid）
merged_data <- inner_join(data1, data2, by = "participant.eid")

# 创建输出文件夹（如果不存在）
if(!dir.exists("output/2025.8.11")){
  dir.create("output/2025.8.11", recursive = TRUE)
}

# 保存到新的 CSV 文件
write_csv(merged_data, "output/2025.8.11/merged_Date5_with_filtered_CVD.csv")

cat("合并完成，结果已保存到 output/2025.8.11/merged_Date5_with_filtered_CVD.csv\n")

#查看字段6153用药情况
library(readr)
library(dplyr)

# 文件路径
input_file <- "output/2025.8.11/merged_Date5_with_filtered_CVD.csv"
output_file <- "output/2025.8.11/participant_p6153_i0_distribution.csv"

# 读取数据
df <- read_csv(input_file)

# 统计数值构成
result <- df %>%
  count(participant.p6153_i0) %>%
  mutate(percent = round(100 * n / sum(n), 2)) %>%
  arrange(desc(n))

# 保存结果为 CSV
write_csv(result, output_file)

cat("统计完成，结果已保存到：", output_file, "\n")

#第二步，提取字段2966高血压的病人，以及是否服用降压药6153中回答为2的
# 稳健生成 hypertension_status 和 antihypertensive_med
# 读取与处理
library(readr)
library(dplyr)
library(stringr)
library(purrr)

infile  <- "output/2025.8.11/merged_Date5_with_filtered_CVD.csv"
outfile <- "output/2025.8.11/participant_with_htn_med_flags.csv"

# 以字符读入，避免 [] 被转义
df <- read_csv(infile, col_types = cols(.default = col_character()), show_col_types = FALSE)

df2 <- df %>%
  mutate(
    # 高血压状态：患病年龄>0 记 1
    age_num = suppressWarnings(as.numeric(participant.p2966_i0)),
    hypertension_status = if_else(!is.na(age_num) & age_num > 0, 1L, 0L),
    
    # 用药：只要方括号里/任意分隔里出现独立的“2”就算服用降压药
    med_raw = ifelse(is.na(participant.p6153_i0), "", participant.p6153_i0),
    # 抽取所有数字 token，再判断是否包含 "2"
    med_has2 = map_lgl(str_extract_all(med_raw, "\\d+"), ~ any(.x == "2")),
    antihypertensive_med = if_else(med_has2, 1L, 0L)
  ) %>%
  select(-age_num, -med_raw, -med_has2)

write_csv(df2, outfile)
cat("已生成：", outfile, "\n")

#第三步读取字段2966字段中为高血压的病人，比较两组病人用药和未用药，cvd的患病情况
library(readr)
library(dplyr)

# 输入文件路径
infile <- "output/2025.8.11/participant_with_htn_med_flags.csv"
outfile <- "output/2025.8.11/CVD_by_RPL_HTNmed.csv"

# 读取数据
df <- read_csv(infile, show_col_types = FALSE) %>%
  mutate(
    is_RPL = suppressWarnings(as.integer(is_RPL)),
    is_CVD = suppressWarnings(as.integer(is_CVD)),
    hypertension_status = suppressWarnings(as.integer(hypertension_status)),
    antihypertensive_med = suppressWarnings(as.integer(antihypertensive_med))
  )

# 只看有高血压的人
summary_tbl <- df %>%
  filter(hypertension_status == 1, !is.na(is_RPL), !is.na(antihypertensive_med)) %>%
  group_by(is_RPL, antihypertensive_med) %>%
  summarise(
    n = n(),
    CVD_cases = sum(is_CVD == 1, na.rm = TRUE),
    CVD_rate = CVD_cases / n,
    .groups = "drop"
  ) %>%
  arrange(is_RPL, desc(antihypertensive_med))

# 写入 CSV
write_csv(summary_tbl, outfile)

cat("结果已保存到：", outfile, "\n")

