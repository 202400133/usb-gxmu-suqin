#第一步计算随访时间，以年为单位，保留小数点后2位
library(readr)
library(dplyr)
library(lubridate)


# liaoming：E:\nas\ukbA\data\suqin\RPL and CVD原始数据\output2025.8.4
# liaoming：file_path <- "E:/nas/ukbA/data/suqin/RPL and CVD原始数据/output2025.8.4/filtered_postbaseline_CVD_with_isCVD.csv"


# 读取数据
file_path <- "output/2025.8.4/filtered_postbaseline_CVD_with_isCVD.csv"
df <- read_csv(file_path)

# 转换日期列格式
df <- df %>%
  mutate(
    CVD_date = as.Date(CVD_date),
    participant.p53_i0 = as.Date(participant.p53_i0),
    # 设置随访截止日期
    end_date = as.Date("2022-10-31")
  )


## liaoming 
df$CVD_date_new<-ifelse(df$CVD_date > df$end_date,
                        df$end_date,
                        df$CVD_date) %>% as.Date()
df[,c("CVD_date","CVD_date_new","participant.p53_i0","end_date","is_CVD", "is_RPL", "follow_up")] %>%View()

df <- df %>%
  mutate(
    follow_up_new = round(
      as.numeric(
        if_else(
          !is.na(CVD_date_new),
          CVD_date_new - participant.p53_i0,
          end_date - participant.p53_i0
        )
      ) / 365.25, 2  # 保留两位小数
    )
  )



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

colnames(df)
#df[,c("CVD_date","participant.p53_i0","end_date","is_CVD", "is_RPL", "follow_up")] %>%View()
df %>% mutate(check=df$participant.p53_i0>df$CVD_date) %>% select(check)%>% sum(na.rm = TRUE)

# 保存回原文件
# write_csv(df, file_path)


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
cox_model <- coxph(Surv(follow_up_new, is_CVD) ~ is_RPL, data = df)

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

# 加载所需包
library(cmprsk)
library(readr)

# 读取数据
file_path <- "output/2025.8.4/filtered_postbaseline_CVD_with_isCVD.csv"
df <- read_csv(file_path)

# 设置变量
ftime <- df$follow_up_new
fstatus <- df$is_CVD
group <- df$is_RPL

# 计算累积发病率
ci <- cuminc(ftime = ftime, fstatus = fstatus, group = group)

# 绘图（设置 y 轴范围为 0 到 0.25）
plot(ci,
     col = c("blue", "red"),
     lty = 1:2,
     xlab = "time (years)",
     ylab = "Cumulative Incidence of CVD",
     main = "Cumulative Incidence of CVD by RPL Status",
     ylim = c(0, 0.25))   # 👈 设置 y 轴范围

# 图例放在左上角
legend("topleft",
       legend = c("Non-RPL", "RPL"),
       col = c("blue", "red"),
       lty = 1:2)


<<<<<<< Updated upstream
##### plot and save #####
df[,c("CVD_date","participant.p53_i0","end_date","is_CVD", "is_RPL", "follow_up")] %>%str()
df$is_CVD<-as.factor(df$is_CVD);df$is_RPL<-as.factor(df$is_RPL)
table(df$is_RPL,df$is_CVD)
#df[,c("CVD_date","participant.p53_i0","end_date","is_CVD", "is_RPL", "follow_up")] %>%View()

median(df$follow_up, na.rm = TRUE)
median(df$follow_up_new, na.rm = TRUE)

library(tidycmprsk)
library(survminer)
library(ggsurvfit)
library(survival)
png("C:/Users/Administrator/Documents/GitHub/usb-gxmu-suqin/output/cvd_incidence_RPL.png",width = 960, height = 480)
cuminc(Surv(follow_up_new, is_CVD) ~ is_RPL, df) %>%
  ggcuminc(outcome = "1") +
  add_confidence_interval() +
  #add_risktable() +
  scale_ggsurvfit()
dev.off()


# 单因素Cox回归分析，检验RPL是否为HCC发病的预测因素
cox_fit <- coxph(Surv(follow_up_new, as.numeric(is_CVD)) ~ is_RPL, data = df) #可以校正年龄等多因素分析
summary(cox_fit)

=======
>>>>>>> Stashed changes
