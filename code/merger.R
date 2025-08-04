# 第一步：读取原始数据，路径从 ../data/ 中读取
Date1 <- read.csv("../data/Date1.csv")
Date2 <- read.csv("../data/Date2.csv")
Date3 <- read.csv("../data/Date3.csv")
Date4 <- read.csv("../data/Date4.csv")

# 确保列名一致（统一为 participant.eid）
names(Date1)[1] <- "participant.eid"
names(Date2)[1] <- "participant.eid"
names(Date3)[1] <- "participant.eid"
names(Date4)[1] <- "participant.eid"

# 合并数据
merged_data <- Reduce(function(x, y) merge(x, y, by = "participant.eid", all = TRUE),
                      list(Date1, Date2, Date3, Date4))

# 查看维度
dim(merged_data)

# 保存合并结果到 data 文件夹
write.csv(merged_data, "../data/merged_data.csv", row.names = FALSE)

# 第二步：剔除 participant.p22001 为 1 和 NA 的行
merged_data_filtered <- merged_data
num_1 <- sum(merged_data_filtered$participant.p22001 == 1, na.rm = TRUE)
num_na <- sum(is.na(merged_data_filtered$participant.p22001))
merged_data_filtered <- merged_data_filtered[!is.na(merged_data_filtered$participant.p22001) &
                                               merged_data_filtered$participant.p22001 != 1, ]
cat("删除值为 1 的行数：", num_1, "\n")
cat("删除 NA 的行数：", num_na, "\n")
cat("筛选后的总行数：", nrow(merged_data_filtered), "\n")
write.csv(merged_data_filtered, "../data/filtered_merged_data.csv", row.names = FALSE)

# 第三步：剔除 participant.p41267 为 0 和 99 的行
filtered_data <- read.csv("../data/filtered_merged_data.csv")
num_zero <- sum(filtered_data$participant.p41267 == 0, na.rm = TRUE)
num_99 <- sum(filtered_data$participant.p41267 == 99, na.rm = TRUE)
filtered_data_new <- filtered_data[!(filtered_data$participant.p41267 %in% c(0, 99)), ]
cat("剔除 participant.p41267 == 0 的行数：", num_zero, "\n")
cat("剔除 participant.p41267 == 99 的行数：", num_99, "\n")
cat("剔除后剩余的总行数：", nrow(filtered_data_new), "\n")
write.csv(filtered_data_new, "../data/filtered_merged_data_p41267_filtered.csv", row.names = FALSE)

# 第四步：剔除 participant.p3839 为 -3 和 -1 的行
filtered_data <- read.csv("../data/filtered_merged_data_p41267_filtered.csv")
num_minus3 <- sum(filtered_data$participant.p3839 == -3, na.rm = TRUE)
num_minus1 <- sum(filtered_data$participant.p3839 == -1, na.rm = TRUE)
filtered_data_new <- filtered_data[!(filtered_data$participant.p3839 %in% c(-3, -1)), ]
cat("剔除 participant.p3839 == -3 的行数：", num_minus3, "\n")
cat("剔除 participant.p3839 == -1 的行数：", num_minus1, "\n")
cat("剔除后剩余的总行数：", nrow(filtered_data_new), "\n")
write.csv(filtered_data_new, "../data/filtered_merged_data_p3839_filtered.csv", row.names = FALSE)
