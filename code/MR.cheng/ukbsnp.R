library(data.table)

# 直接读取 .bgz 文件
#gwas_data <- fread("E:\\01rsa\\I37 Pulmonary valve disorders\\icd10-I37-both_sexes.tsv.bgz")

# 查看前几行
#head(gwas_data)

setwd("E:\\01rsa\\I31 Other diseases of pericardium")
library(data.table)
library(dplyr)

# 1. 加载 SNP 信息
snp_info <- as.data.table(readRDS("E:\\01rsa\\I71 Aortic aneurysm and dissection\\ukb_snp.rds"))

# 2. 加载 GWAS 数据
gwas_data <- fread("E:/01rsa/I31 Other diseases of pericardium/icd10-I31-both_sexes.tsv.bgz")

# 3. 从 variant 列提取染色体和位置（如果需要）
#gwas_data[, c("chr", "pos", "alleles") := tstrsplit(variant, "[:_]")]
#gwas_data[, pos := as.numeric(pos)]"E:/01rsa/I31 Other diseases of pericardium/icd10-I31-both_sexes.tsv.bgz"

# 4. 合并数据
merged_data <- left_join(gwas_data, snp_info, by = c("chr" = "chrom", "pos" = "pos"))

# 5. 保存结果
fwrite(merged_data, "I31_gwas_with_snp_info.tsv", sep = "\t")
 