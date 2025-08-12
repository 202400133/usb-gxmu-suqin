# 加载必要的包
library(readr)
library(dplyr)

# 读取CSV文件
file_path <- "output/2025.8.4/filtered_postbaseline_CVD_with_isCVD.csv"
data <- read_csv(file_path)

# 检查字段是否存在
if (!"participant.p3839_i0" %in% colnames(data)) {
  stop("字段 'participant.p3839_i0' 不存在于数据中。请检查列名是否正确。")
}

# 统计流产次数的频数
abortion_counts <- data %>%
  count(participant.p3839_i0, name = "频数") %>%
  rename("流产次数" = participant.p3839_i0) %>%
  arrange(流产次数)  # 按流产次数排序

# 打印结果
print(abortion_counts)

# 可选：保存结果到新文件
write_csv(abortion_counts, "output/2025.8.10/abortion_counts_summary.csv") 