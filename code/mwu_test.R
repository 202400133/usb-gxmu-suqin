# 加载必要包
library(dplyr)

# 设置文件路径
file_path <- "output/2025.8.4/filtered_postbaseline_CVD_with_isCVD.csv"

# 读取数据
df <- read.csv(file_path)

# 检查分组变量
if (!"is_RPL" %in% colnames(df)) stop("❌ 缺少 is_RPL 分组列")

# 目标列
target_vars <- c(
  "participant.p30870_i0", "participant.p30790_i0", "participant.p30700_i0",
  "participant.p30620_i0", "participant.p30650_i0", 
  "participant.p30730_i0", "participant.p30710_i0"
)

# 创建结果表
results <- data.frame(
  Variable = character(),
  RPL = character(),
  nonRPL = character(),
  P_value = character(),
  stringsAsFactors = FALSE
)

# 定义中位数 + 四分位数格式函数
median_q_format <- function(x) {
  q <- quantile(x, probs = c(0.25, 0.5, 0.75))
  sprintf("%.2f[%.2f,%.2f]", q[2], q[1], q[3])
}

# 循环每个变量
for (var in target_vars) {
  if (var %in% colnames(df)) {
    sub_df <- df[, c("is_RPL", var)] %>% na.omit()
    group0 <- sub_df[sub_df$is_RPL == 0, var]
    group1 <- sub_df[sub_df$is_RPL == 1, var]
    
    # 中位数格式化
    med0 <- median_q_format(group0)
    med1 <- median_q_format(group1)
    
    # Mann-Whitney U 检验
    p_val <- tryCatch(wilcox.test(group0, group1)$p.value, error = function(e) NA)
    p_fmt <- ifelse(is.na(p_val), NA, sprintf("%.3f", p_val))
    
    # 添加结果
    results <- rbind(results, data.frame(
      Variable = var,
      RPL = med1,
      nonRPL = med0,
      P_value = p_fmt,
      stringsAsFactors = FALSE
    ))
  } else {
    results <- rbind(results, data.frame(
      Variable = var,
      RPL = "列不存在",
      nonRPL = "列不存在",
      P_value = NA,
      stringsAsFactors = FALSE
    ))
  }
}

# 查看结果
print(results)

#保存结果为 CSV
write.csv(results, "output/2025.8.4/mwu_test_RPL_vs_nonRPL.csv", row.names = FALSE)
