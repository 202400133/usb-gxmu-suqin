# 加载所需包
library(dplyr)

# 设置文件路径
file_path <- "output/2025.8.4/filtered_postbaseline_CVD_with_isCVD.csv"

# 读取数据
df <- read.csv(file_path)

# 指定要进行 t 检验的变量列
target_vars <- c(
  "participant.p21022", "participant.p4080_i0_a1", "participant.p4079_i0_a1",
  "participant.p2966_i0", "participant.p21001_i0", "participant.p1239_i0",
  "participant.p30690_i0", "participant.p30760_i0", "participant.p30870_i0",
  "participant.p30780_i0", "participant.p30790_i0", "participant.p30630_i0",
  "participant.p30640_i0", "participant.p30750_i0", "participant.p30700_i0",
  "participant.p30600_i0", "participant.p30620_i0", "participant.p30650_i0",
  "participant.p30730_i0", "participant.p30710_i0", "participant.p2784_i0",
  "participant.p2814_i0", "participant.p2834_i0", "participant.p3581_i0",
  "participant.p2714_i0"
)

# 检查是否存在 is_RPL 分组列
if (!"is_RPL" %in% colnames(df)) stop("⚠️ 缺少 is_RPL 分组列")

# 创建结果数据框
results <- data.frame(
  Variable = character(),
  Mean_SD_RPL = character(),
  Mean_SD_nonRPL = character(),
  P_value = character(),
  stringsAsFactors = FALSE
)

# 自定义 p 值格式化函数
format_pval <- function(p) {
  if (is.na(p)) return(NA)
  if (p < 0.001) return("<0.001")
  return(sprintf("%.3f", p))
}

# 循环处理每个变量
for (var in target_vars) {
  if (var %in% colnames(df)) {
    sub_df <- df[, c("is_RPL", var)] %>% na.omit()
    group0 <- sub_df[sub_df$is_RPL == 0, var]
    group1 <- sub_df[sub_df$is_RPL == 1, var]
    
    mean_sd_0 <- sprintf("%.2f ± %.2f", mean(group0), sd(group0))
    mean_sd_1 <- sprintf("%.2f ± %.2f", mean(group1), sd(group1))
    
    p_val <- tryCatch(t.test(group1, group0)$p.value, error = function(e) NA)
    formatted_p <- format_pval(p_val)
    
    results <- rbind(results, data.frame(
      Variable = var,
      Mean_SD_RPL = mean_sd_1,
      Mean_SD_nonRPL = mean_sd_0,
      P_value = formatted_p,
      stringsAsFactors = FALSE
    ))
  } else {
    results <- rbind(results, data.frame(
      Variable = var,
      Mean_SD_RPL = "列不存在",
      Mean_SD_nonRPL = "列不存在",
      P_value = NA,
      stringsAsFactors = FALSE
    ))
  }
}

# 保存为 CSV 文件
write.csv(results, "output/2025.8.4/t_test_RPL_vs_nonRPL_results.csv", row.names = FALSE)