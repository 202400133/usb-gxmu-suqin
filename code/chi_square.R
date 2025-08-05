# 加载必要包
library(dplyr)

# 读取数据
file_path <- "output/2025.8.4/filtered_postbaseline_CVD_with_isCVD.csv"
df <- read.csv(file_path)

# 检查 is_RPL 是否存在
if (!"is_RPL" %in% colnames(df)) stop("❌ 缺少 is_RPL 分组列")

# 定义变量分组
binary_vars <- c("participant.p2784_i0", "participant.p2814_i0", "participant.p2834_i0")
smoking_var <- "participant.p20116_i0"
hypertension_var <- "participant.p2966_i0"

# 统一结果表
results <- data.frame(
  Variable = character(),
  RPL = character(),
  nonRPL = character(),
  P_value = character(),
  stringsAsFactors = FALSE
)

# --- 1. 二分类变量处理 ---
for (var in binary_vars) {
  if (var %in% colnames(df)) {
    temp <- df %>% filter(!is.na(.data[[var]]) & .data[[var]] %in% c(0,1))
    tab <- table(temp[[var]], temp$is_RPL)
    
    # 计算频率
    count_table <- as.data.frame.matrix(tab)
    if (all(c("0", "1") %in% rownames(count_table))) {
      rpl_1 <- count_table["1", "1"]
      rpl_total <- sum(count_table[, "1"])
      nonrpl_1 <- count_table["1", "0"]
      nonrpl_total <- sum(count_table[, "0"])
      
      rpl_str <- sprintf("%d(%.2f%%)", rpl_1, rpl_1 / rpl_total * 100)
      nonrpl_str <- sprintf("%d(%.2f%%)", nonrpl_1, nonrpl_1 / nonrpl_total * 100)
      
      # 卡方检验
      p <- chisq.test(tab)$p.value
      p_fmt <- sprintf("%.3f", p)
      
      results <- rbind(results, data.frame(
        Variable = var,
        RPL = rpl_str,
        nonRPL = nonrpl_str,
        P_value = p_fmt,
        stringsAsFactors = FALSE
      ))
    }
  }
}

# --- 2. 吸烟状况（0=never, 1=previous, 2=current）按分类分三行显示 ---
if (smoking_var %in% colnames(df)) {
  temp <- df %>% filter(.data[[smoking_var]] %in% c(0,1,2))
  tab <- table(temp[[smoking_var]], temp$is_RPL)
  
  # 获取总数
  rpl_total <- sum(tab[, "1"])
  nonrpl_total <- sum(tab[, "0"])
  
  # 获取分类标签
  smoke_labels <- c("Never", "Previous", "Current")
  
  # 卡方检验整体 p 值
  p <- chisq.test(tab)$p.value
  p_fmt <- sprintf("%.3f", p)
  
  for (code in 0:2) {
    rpl_n <- tab[as.character(code), "1"]
    nonrpl_n <- tab[as.character(code), "0"]
    
    rpl_str <- sprintf("%d(%.2f%%)", rpl_n, rpl_n / rpl_total * 100)
    nonrpl_str <- sprintf("%d(%.2f%%)", nonrpl_n, nonrpl_n / nonrpl_total * 100)
    
    results <- rbind(results, data.frame(
      Variable = paste0(smoking_var, "_", smoke_labels[code + 1]),
      RPL = rpl_str,
      nonRPL = nonrpl_str,
      P_value = p_fmt,
      stringsAsFactors = FALSE
    ))
  }
}
# --- 3. 高血压处理：NA/-1/-3 → 0，无高血压；其他为 1，有高血压 ---
if (hypertension_var %in% colnames(df)) {
  df$hypertension_status <- NA
  df$hypertension_status[is.na(df[[hypertension_var]]) | df[[hypertension_var]] %in% c(-1, -3)] <- 0
  df$hypertension_status[!is.na(df[[hypertension_var]]) & !(df[[hypertension_var]] %in% c(-1, -3))] <- 1
  
  temp <- df %>% filter(!is.na(hypertension_status))
  tab <- table(temp$hypertension_status, temp$is_RPL)
  
  count_table <- as.data.frame.matrix(tab)
  if (all(c("1", "0") %in% rownames(count_table))) {
    rpl_1 <- count_table["1", "1"]
    rpl_total <- sum(count_table[, "1"])
    nonrpl_1 <- count_table["1", "0"]
    nonrpl_total <- sum(count_table[, "0"])
    
    rpl_str <- sprintf("%d(%.2f%%)", rpl_1, rpl_1 / rpl_total * 100)
    nonrpl_str <- sprintf("%d(%.2f%%)", nonrpl_1, nonrpl_1 / nonrpl_total * 100)
    
    p <- chisq.test(tab)$p.value
    p_fmt <- sprintf("%.3f", p)
    
    results <- rbind(results, data.frame(
      Variable = hypertension_var,
      RPL = rpl_str,
      nonRPL = nonrpl_str,
      P_value = p_fmt,
      stringsAsFactors = FALSE
    ))
  }
}

# 显示结果
print(results)

#导出结果
write.csv(results, "output/2025.8.4/chi_square_RPL_vs_nonRPL_results.csv", row.names = FALSE)
