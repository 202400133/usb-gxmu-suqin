setwd("E:\\01rsa\I31 Other diseases of pericardium")

library(ggplot2)
library(dplyr)
library(data.table)
library(devtools)
library(MendelianRandomization)
library(TwoSampleMR) 
library(MRPRESSO) 

library(ieugwasr)

EXP1<- data.table::fread("E:/01rsa/GWAS3239672/N96_meta_finnest_22062022.tsv")
EXP1<- data.frame(EXP1)
head(EXP1)

colnames(EXP1)[c(3,4,9,8,5,6,7)] <-
  c("effect_allele.exposure", "other_allele.exposure", "rsid", "pval.exposure", "beta.exposure", "se.exposure", "eaf.exposure")

EXP1$id.exposure <- "Recurrent_Spontaneous_Abortion"
EXP1$exposure <- "Recurrent_Spontaneous_Abortion"
EXP1$samplesize.exposure <- 150965
head(EXP1, 6)

EXP2 <- subset(EXP1, pval.exposure < 5e-06)


# 删除 rsid 为 NA 或空字符串的行
EXP2 <- EXP2 %>%
  filter(rsid != "")# 删除空字符串 ""
#EXP2 <- EXP2 %>%
  #filter(nearest_genes != "")

# 测试 API 连通性
#httr::GET("https://api.opengwas.io/api/ld/clump")$status_code

#EXP2 <- clump_data(EXP2, clump_kb = 10000,clump_r2 = 0.001, pop = "EUR")

colnames(EXP2)[c(8)] <-
  c("pval")

# 方法2：使用ieugwasr包的本地模式
library(ieugwasr)

clumped_EXP2 <- ld_clump(
  dat = EXP2,
  bfile = "E:/01rsa/EUR/1000G_EUR",
  plink = "E:/plink_win64_20250615/plink.exe",
  clump_kb = 10000,
  clump_r2 = 0.001
)

#OUT <- data.table::fread("E:/01rsa/I71 Aortic aneurysm and dissection/I71_gwas_with_snp_info.tsv")

# 读取 TSV 文件（自动处理 NA）
OUT <- data.table::fread("E:/01rsa/I44 Atrioventricular and left bundle-branch block/I44_gwas_with_snp_info.tsv", na.strings = c("", "NA"))
# 删除所有包含 NA 的行
OUT <- na.omit(OUT)
#OUT <- data.frame(OUT)
head(OUT)


#data.table::fwrite(OUT, "E:/01rsa/I71 Aortic aneurysm and dissection/I71_processed.tsv", 
                   #sep = "\t", na = "NA", row.names = FALSE)


colnames(OUT)[c(3,4,5,7,8)] <-
  c("other_allele.outcome", "effect_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome")

#colnames(OUT)[c(9)] <-
  #c("neglog10_pval_meta")
OUT <- OUT %>%
   mutate(pval.outcome = 10^(-neglog10_pval_meta))

OUT$id.outcome <- "Atrioventricular and left bundle-branch block"
OUT$outcome <- "Atrioventricular and left bundle-branch block"
OUT$samplesize.outcome <- 424099
head(OUT)

colnames(clumped_EXP2)[c(8)] <-
  c("pval.exposure")

total <- merge(OUT, clumped_EXP2, by.x = "rsid", by.y = "rsid", all = F)

total <- subset(total, pval.outcome > 5e-08)
total <- total[!duplicated(total$rsid), ]



EXP3 <- total[, c("rsid", "chr", "pos", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure", "beta.exposure", "se.exposure", "pval.exposure", "id.exposure", "exposure", "samplesize.exposure")]
OUT3 <- total[, c("rsid", "chr", "pos", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome", "pval.outcome", "id.outcome", "outcome", "samplesize.outcome")]

colnames(EXP3)[c(1)] <-
  c("SNP")
colnames(OUT3)[c(1)] <-
  c("SNP")

data_h <-harmonise_data(exposure_dat = EXP3, outcome_dat = OUT3, action = 2)

res <- mr(data_h)

OR <- generate_odds_ratios(res)
OR

#通过调参进行限定方法的MR检测
mr(data_h, method_list = c("mr_egger_regression", "mr_ivw"))

RE <-mr(data_h,method_list=c('mr_ivw_mre'))
REOR <-generate_odds_ratios(RE)

FE <-mr(data_h,method_list=c('mr_ivw_fe'))
FEOR <-generate_odds_ratios(FE)

res_single <- mr_singlesnp(data_h)
ORR <-generate_odds_ratios(res_single)
#异质性检验
het <-mr_heterogeneity(data_h)
het

#离群值检验
mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure",
          OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = data_h, NbDistribution = 10000, 
          SignifThreshold = 0.05)
#水平多效性检验（独立性）
pleio <- mr_pleiotropy_test(data_h)
pleio

#逐个剔除检验留一法图
single <- mr_leaveoneout(data_h)
mr_leaveoneout_plot(single)

#散点图
mr_scatter_plot(res,data_h)

#森林图
res_single <- mr_singlesnp(data_h)
mr_forest_plot(res_single)

#漏斗图
mr_funnel_plot(res_single)


# 计算每个SNP的F统计量
data_h$f_stat <- (data_h$beta.exposure / data_h$se.exposure)^2

# 计算平均F统计量
mean_f_stat <- mean(data_h$f_stat)
cat("Mean F-statistic:", mean_f_stat, "\n")

# 合并F统计量到工具变量数据
iv_data <- data_h[, c(
  "SNP", "chr.x", "pos.x", 
  "effect_allele.exposure", "other_allele.exposure",
  "beta.exposure", "se.exposure", "pval.exposure", 
  "eaf.exposure", "samplesize.exposure", "f_stat"
)]

# 添加平均F统计量作为备注
iv_data$note <- ifelse(
  row.names(iv_data) == 1, 
  paste0("Mean_F=", round(mean_f_stat, 2)), 
  ""
)

# 保存工具变量数据 + F统计量
write.csv(
  iv_data,
  file = "I44_3239672.csv",
  row.names = FALSE,
  quote = FALSE  # 避免字符串带引号
)