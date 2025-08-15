#logistics回归分析第一步合并新下载的字段
library(dplyr)
library(readr)

# 文件路径
file_main <- "output/2025.8.4/filtered_merged_data_with_CVD_info.csv"
file_new  <- "data/Date5.csv"

# 读取数据
df_main <- read_csv(file_main)
df_new  <- read_csv(file_new)

# 合并数据（按 participant.eid）
df_merged <- df_main %>%
  left_join(df_new, by = "participant.eid")

# 创建保存目录（如果不存在）
output_dir <- "output/2025.8.14"
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 保存到新文件
output_file <- file.path(output_dir, "merged_filtered_with_Date5.csv")
write_csv(df_merged, output_file)

cat("合并完成，文件已保存到：", output_file, "\n")

#第二步提取高血压患病状态
library(dplyr)
library(readr)

# 读取刚才生成的合并表
file_path <- "output/2025.8.14/merged_filtered_with_Date5.csv"
df <- read_csv(file_path)

# 添加 hypertension_status 列
df <- df %>%
  mutate(
    hypertension_status = ifelse(!is.na(participant.p2966_i0), 1, 0)
  )

# 保存回原文件（覆盖）
write_csv(df, file_path)

cat("已在表格中加入 hypertension_status 列，并保存到原文件：", file_path, "\n")

#第三步RPL分组
library(dplyr)
library(readr)

# 路径
file_path <- "output/2025.8.14/merged_filtered_with_Date5.csv"

# 读取
df <- read_csv(file_path, show_col_types = FALSE)

# 防御：确认字段存在
stopifnot("participant.p3839_i0" %in% names(df))

# 新增 is_RPL 列：p3839 >= 3 记为 1，否则（含 NA）记为 0
df <- df %>%
  mutate(
    .p3839_num = suppressWarnings(as.numeric(`participant.p3839_i0`)),
    is_RPL = dplyr::if_else(!is.na(.p3839_num) & .p3839_num >= 3, 1L, 0L)
  ) %>%
  select(-.p3839_num)

# 覆盖保存到原文件
write_csv(df, file_path)

# 可选：快速查看分组计数
cat("is_RPL 分组计数：\n")
print(df %>% count(is_RPL))
cat("\n已写回：", file_path, "\n")

#第四步把22个cvd亚组疾病变为二分类
# 加载必要库
library(dplyr)
library(readr)

# 读取原始数据
df <- read_csv("output/2025.8.14/merged_filtered_with_Date5.csv")

# 要处理的 CVD 日期列（22个）
cvd_date_vars <- c(
  "participant.p131382", "participant.p131350", "participant.p131352", "participant.p131342",
  "participant.p131344", "participant.p131354", "participant.p131296", "participant.p131298",
  "participant.p131304", "participant.p131306", "participant.p131314", "participant.p131316",
  "participant.p131386", "participant.p131366", "participant.p131362", "participant.p131056",
  "participant.p131322", "participant.p131324", "participant.p131328", "participant.p131400",
  "participant.p131308", "participant.p131348"
)

# 遍历每个列，生成对应 s_ 开头的新列
for (col in cvd_date_vars) {
  new_col <- paste0("is_", col)  # 新列名
  df[[new_col]] <- ifelse(!is.na(df[[col]]), 1, 0)
}

# 可选：覆盖原文件（如需保存）
write_csv(df, "output/2025.8.14/merged_filtered_with_Date5.csv")

