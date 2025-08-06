#第一步计算随访时间，以年为单位，保留小数点后2位
library(readr)
library(dplyr)
library(lubridate)

# 读取数据
file_path <- "output/2025.8.4/filtered_postbaseline_CVD_with_isCVD.csv"
df <- read_csv(file_path)

# 转换日期列格式
df <- df %>%
  mutate(
    CVD_date = as.Date(CVD_date),
    participant.p53_i0 = as.Date(participant.p53_i0)
  )

# 设置随访截止日期
end_date <- as.Date("2022-10-31")

# 计算随访时间
df <- df %>%
  mutate(
    follow_up = round(
      as.numeric(
        if_else(
          !is.na(CVD_date),
          CVD_date - participant.p53_i0,
          end_date - participant.p53_i0
        )
      ) / 365.25, 2  # 保留两位小数
    )
  )

# 保存回原文件
write_csv(df, file_path)

#第二步查看随访时间中位数
library(readr)
library(dplyr)

# 读取数据
file_path <- "output/2025.8.4/filtered_postbaseline_CVD_with_isCVD.csv"
df <- read_csv(file_path)

# 计算 follow_up 中位数（忽略缺失值）
median_follow_up <- median(df$follow_up, na.rm = TRUE)

# 输出结果
print(paste("Median follow-up time (years):", round(median_follow_up, 2)))

#第三步
#计算cox回归分析
library(readr)
library(survival)
library(dplyr)

# 读取数据
file_path <- "output/2025.8.4/filtered_postbaseline_CVD_with_isCVD.csv"
df <- read_csv(file_path)

# 拟合 Cox 模型（仅含 is_RPL）
cox_model <- coxph(Surv(follow_up, is_CVD) ~ is_RPL, data = df)

# 提取 HR、95% CI、P 值
hr <- exp(coef(cox_model))
ci <- exp(confint(cox_model))
p_val <- summary(cox_model)$coefficients[,"Pr(>|z|)"]

# 整理成表格
cox_result <- data.frame(
  Variable = "is_RPL",
  HR = round(hr, 3),
  CI_lower = round(ci[1], 3),
  CI_upper = round(ci[2], 3),
  P_value = round(p_val, 4)
)

# 保存结果为 CSV
write_csv(cox_result, "output/2025.8.5/cox_RPL_vs_CVD.csv")

#用cuminc函数做可视化图
# 如尚未安装，请先安装一次
install.packages("cmprsk")

library(cmprsk)
library(readr)

# 读取数据
file_path <- "output/2025.8.4/filtered_postbaseline_CVD_with_isCVD.csv"
df <- read_csv(file_path)

# 提取所需列
ftime <- df$follow_up         # 随访时间
fstatus <- df$is_CVD          # 是否发生CVD（事件=1，删失=0）
group <- df$is_RPL            # 分组变量（RPL vs 非RPL）

# 第四步使用 cuminc 函数进行分组累计发生概率分析
ci <- cuminc(ftime = ftime, fstatus = fstatus, group = group)

# 绘图
plot(ci,
     col = c("blue", "red"),
     lty = 1:2,
     xlab = "Follow-up time (years)",
     ylab = "Cumulative Incidence of CVD",
     main = "Cumulative Incidence of CVD by RPL Status")

legend("bottomright",
       legend = c("Non-RPL", "RPL"),
       col = c("blue", "red"),
       lty = 1:2)

